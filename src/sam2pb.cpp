#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <regex>
#include <vector>
#include <numeric>
#include <random>
#include <iomanip>

#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>
#include <tbb/parallel_for.h>

#include "dataset.hpp"
#include "read.hpp"
#include "timer.hpp"

#include "sam2pb.hpp"
#include "sam.pb.h"

#include "config.hpp"  

// freyja - includes indels, litreally every single possible mutation of a chunk along with frequency
// freyja-depth - how many reads/chunks were there starting at a given nucleotide
// vcf - for all non indels with a mutation, look at every single read and tell which mutation it corresponds to
// also for pairs of mutation frequencies
static double frequency_read_cutoff ;
static int phred_score_cutoff;

std::string CHROM;


static void 
dump_sub_table(const sub_table& sub, std::string const& name) {
    std::ofstream fout(name);
    fout << "Position\tAllele\tFrequency\tDepth\n";
    for (size_t i = 0; i < sub.size(); ++i) {
        int sum = std::accumulate(sub[i].begin(), sub[i].end(), 0);
        
        std::vector<std::pair<int, int>> res;
        for (size_t j = 0; j < 6; ++j) {
            if (sub[i][j]) {
                res.emplace_back(sub[i][j], j);
            }
        }
        
        std::sort(res.begin(), res.end(), std::greater<>());
        for (auto [_, j] : res) {
            fout << i + 1 << '\t'
                 << GENOME_STRING[j] << '\t'
                 << std::fixed << std::setprecision(10) << (double) sub[i][j] / sum << '\t'
                 << sum << '\n';
        }
    }
}

void sam2PB(const dataset& d) {
    std::string proto_filename = d.pb_path();
    std::string sam_file = d.sam_path();
    phred_score_cutoff = d.min_phred();
    frequency_read_cutoff = d.min_af();

    //Reading reference genome
    timer timer;
    timer.start();
    
    std::string ref_seq = d.reference();
    CHROM = d.reference_name();

    //std::ofstream outfile_freyja_vcf(freyja_vcf_file, std::ios::out | std::ios::binary);
    //std::ofstream outfile_freyja_depth(freyja_depth_file, std::ios::out | std::ios::binary);
    //boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_freyja_vcf, outbuf_freyja_depth;
    
    //if (freyja_vcf_file.find(".gz\0") != std::string::npos)
    //    outbuf_freyja_vcf.push(boost::iostreams::gzip_compressor());
    //if (freyja_depth_file.find(".gz\0") != std::string::npos)
    //    outbuf_freyja_depth.push(boost::iostreams::gzip_compressor());
    
    //outbuf_freyja_vcf.push(outfile_freyja_vcf);
    //outbuf_freyja_depth.push(outfile_freyja_depth);
    //std::ostream freyja_vcf(&outbuf_freyja_vcf);
    //std::ostream freyja_depth(&outbuf_freyja_depth);

    sam sam{ref_seq, (int) d.max_reads()};
    boost::filesystem::ifstream fileHandler(sam_file);
    
    std::vector<std::string> lines;
    std::string s;
    while (getline(fileHandler, s)) {
        lines.emplace_back(std::move(s));
    }

    tbb::queuing_mutex mutex;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, lines.size()), 
        [&](tbb::blocked_range<size_t> range)  {
            sam.add_reads(lines, range, &mutex);
        }
    );

    sam.build();

    sam.dump_proto(proto_filename);
    // sam.dump_freyja(freyja_depth, freyja_vcf);
    fprintf(stderr, "\nFiles generated in %ld sec \n\n", timer.seconds());
}

sam::sam(const std::string& ref, int subsampled_reads) : reference_seq{ref}, subsampled_reads{subsampled_reads} {
    this->frequency_table.resize(ref.size());
    this->collapsed_frequency_table.resize(ref.size());
}


void sam::dump_proto(const std::string& filename) {
    Sam::sam data;

    for (const sam_read& read: aligned_reads) {
        auto dump = data.add_reads();
        dump->set_start_idx(read.start_idx + 1);
        dump->set_content(read.aligned_string);
        dump->set_degree(read.degree);
        dump->set_read(read.degree_name());
    }

    for (const auto& kv: reverse_merge) {
        Sam::column_info *col = data.add_reverse_columns();
        col->set_column_name(kv.first);
        for (const auto& str: kv.second) {
            std::string* dump = col->add_input_columns();
            *dump = str;
        }
    }

    // Boost library used to stream the contents to the output protobuf file in
    // uncompressed or compressed .gz format
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf< boost::iostreams::output> outbuf;

    if (filename.find(".gz\0") != std::string::npos) {
        try {
            outbuf.push(boost::iostreams::gzip_compressor());
            outbuf.push(outfile);
            std::ostream outstream(&outbuf);
            data.SerializeToOstream(&outstream);
            boost::iostreams::close(outbuf);
            outfile.close();
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        data.SerializeToOstream(&outfile);
        outfile.close();
    }
}

void sam::add_reads(const std::vector<std::string> &lines, tbb::blocked_range<size_t> range, tbb::queuing_mutex* mutex)
{
    std::vector<sam_read> reads;

    for (size_t i = range.begin(); i < range.end(); ++i)
    {
        const std::string& line = lines[i];

        std::vector<std::string> tokens;
        std::stringstream ss(line);
        std::string tmp;
        while (ss >> tmp) {
            tokens.push_back(std::move(tmp));
        }

        /* skip if unmapped */
        if (tokens.empty() || tokens[0][0] == '@' || std::stoi(tokens[1]) & 4)
        {
            return;
        }

        /* zero indexed */
        int start_idx = std::stoi(tokens[3]);
        const std::string &read_seq = tokens[9];
        const std::string &phred_seq = tokens[10];

        /* process cigar string, updating frequency tables */
        std::regex cigar_pattern(R"(\d+[A-Za-z])"); // Matches one or more digits followed by a character
        std::sregex_iterator r_iter(tokens[5].begin(), tokens[5].end(), cigar_pattern);
        std::sregex_iterator end;
        std::vector<std::pair<int, char>> cigar_chunks;
        while (r_iter != end)
        {
            std::string sub_cigar = (*r_iter).str();
            std::regex pos_char("([0-9]+)([A-Za-z])");
            std::sregex_iterator regex_iter(sub_cigar.begin(), sub_cigar.end(), pos_char);
            std::smatch pc_match = *regex_iter;
            // Convert the captured substrings to int and char and store in cigar_chunks
            cigar_chunks.emplace_back(std::stoi(pc_match.str(1)), pc_match.str(2)[0]);
            r_iter++;
        }

        int seq_idx = 0, ins_count = 0, del_count = 0;
        std::string ref_nuc, alt_nuc, build;
        for (const auto &[cig_len, cig_val] : cigar_chunks)
        {
            int curr_sub_len = 0;
            while (curr_sub_len < cig_len)
            {
                int nuc_pos = start_idx + seq_idx - ins_count + del_count;
                switch (cig_val)
                {
                case 'I':
                    ref_nuc = "";
                    alt_nuc = "";
                    curr_sub_len += cig_len;
                    seq_idx += cig_len;
                    ins_count += cig_len;
                    /* build stays the same*/
                    break;

                case 'D':
                    alt_nuc = "";
                    ref_nuc += reference_seq.substr(nuc_pos - 1, cig_len);
                    curr_sub_len += cig_len;
                    del_count += cig_len;
                    build += std::string(cig_len, '_');
                    break;

                case 'N':
                    ref_nuc = reference_seq[nuc_pos - 1];
                    alt_nuc = "N";
                    curr_sub_len++;
                    seq_idx++;
                    build += "N";
                    break;

                case 'H':
                    curr_sub_len += cig_len;
                    break;

                case 'S':
                    curr_sub_len += cig_len;
                    seq_idx += cig_len;
                    ins_count += cig_len;
                    break;

                default:
                    ref_nuc = reference_seq[nuc_pos - 1];
                    if ((static_cast<int>(phred_seq[seq_idx]) - 33) < phred_score_cutoff)
                    {
                        alt_nuc = std::string(1, 'N');
                    }
                    else
                        alt_nuc = std::string(1, read_seq[seq_idx]);

                    // replace non ACGTN with N
                    if (GENOME_STRING.find(alt_nuc) == std::string::npos)
                    {
                        alt_nuc = 'N';
                    }

                    build += alt_nuc;
                    curr_sub_len++;
                    seq_idx++;
                }
            }
        }

        reads.emplace_back(sam_read{std::move(tokens[0]), start_idx - 1, 1, std::move(build)});
    }

    tbb::queuing_mutex::scoped_lock lock{*mutex};
    for (sam_read& read : reads) {
        // update frequency table
        for (size_t i = 0; i < read.aligned_string.size(); ++i) {
            size_t pos = i + read.start_idx;
            char current = read.aligned_string[i];
            if (current == 'N') {
                continue;
            }

            int sub = (int)GENOME_STRING.find(current);
            assert(0 <= sub && sub < (int)GENOME_STRING.size());

            frequency_table[pos][sub] += read.degree;
        }

        aligned_reads.push_back(std::move(read));
    }
}

void sam::read_correction() {
    std::vector<int> total_occurences(reference_seq.size());
    std::vector<char> majority(reference_seq.size());
    for (int i = 0; i < (int) reference_seq.size(); ++i) {
        total_occurences[i] = std::accumulate(collapsed_frequency_table[i].begin(), collapsed_frequency_table[i].end(), 0); 
        majority[i] = std::max_element(collapsed_frequency_table[i].begin(), collapsed_frequency_table[i].end()) - collapsed_frequency_table[i].begin();
    }

    /* first correct aligned reads */
    for (int i = 0; i < (int) aligned_reads.size(); ++i) {
        std::string& align = aligned_reads[i].aligned_string;
        int start = aligned_reads[i].start_idx;
        for (int j = 0; j < (int) align.size(); ++j) {
            int indx = start + j;
            int curr = GENOME_STRING.find(align[j]);

            if (frequency_read_cutoff - (double) collapsed_frequency_table[indx][curr] / total_occurences[indx] > SCORE_EPSILON) {
                if (MAP_TO_MAJORITY_INSTEAD_OF_N) {
                    align[j] = GENOME_STRING[majority[indx]];
                }
                else {
                    align[j] = 'N';
                }
            }
        }
    }

    /* update frequency table */
    for (int i = 0; i < (int) reference_seq.size(); ++i) {
        for (int j = 0; j < (int) GENOME_STRING.size(); ++j) {
            if (j == majority[i]) continue;

            if (frequency_read_cutoff - (double) collapsed_frequency_table[i][j] / total_occurences[i] > SCORE_EPSILON) {
                collapsed_frequency_table[i][majority[i]] += std::exchange(collapsed_frequency_table[i][j], 0);
            }
        }
    }
}

void sam::merge_duplicates() {
    assert(!aligned_reads.empty());

    std::vector<sam_read> merged{aligned_reads[0]};

    for (int i = 1; i < (int) aligned_reads.size(); ++i) {
        if (aligned_reads[i] == merged.back()) {
            merged.back().degree++;
        }
        else {
            merged.push_back(aligned_reads[i]);
        }
    }

    /* create reverse table */
    int reverse = 0;
    for (const auto& merge: merged) {
        std::string generated = merge.degree_name();

        for (int j = 0; j < merge.degree; ++j) {
            this->reverse_merge[generated].push_back(aligned_reads[j + reverse].raw_name);
        }
        reverse += merge.degree;
    }

    this->aligned_reads = merged;
}

void sam::subsample() {
    if ((int) this->aligned_reads.size() <= subsampled_reads) {
        return;
    }

    // transform aligned reads to a subset of 
    // aligned reads

    double best_score = std::numeric_limits<double>::max();
    std::vector<sam_read> best_set;

    std::random_device rd;
    std::mt19937 g(rd());

    // reference allele frequency
    auto frequency_vector = [&](const sub_table& t) {
        std::vector<double> p;
        for (size_t i = 0; i < reference_seq.size(); ++i) {
            int sum = std::accumulate(t[i].begin(), t[i].end(), 0);
            for (size_t j = 0; j < GENOME_STRING.size(); ++j) {
                if (!sum) p.push_back(0);
                else p.push_back((double) t[i][j] / sum);
            }
        }

        double sum = std::accumulate(p.begin(), p.end(), 0.0);
        std::for_each(p.begin(), p.end(), 
            [&](double &p) { p /= sum; });

        return p;
    };

    std::vector<double> p = frequency_vector(this->collapsed_frequency_table);
    tbb::queuing_mutex my_mutex;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, SUBSAMPLE_ITERS),
                      [&](tbb::blocked_range<size_t> block)
                      {
                          for (size_t i = block.begin(); i < block.end(); ++i)
                          {
                              std::vector<sam_read> current;

                              std::vector<int> index_set(aligned_reads.size());
                              std::iota(index_set.begin(), index_set.end(), 0);
                              {
                                  tbb::queuing_mutex::scoped_lock my_lock{my_mutex};
                                  std::shuffle(index_set.begin(), index_set.end(), g);
                              }
                              for (int j = 0; j < subsampled_reads; ++j)
                              {
                                  current.push_back(aligned_reads[index_set[j]]);
                              }

                              sub_table af(this->reference_seq.size());

                              for (const auto &read : current)
                              {
                                  for (size_t j = 0; j < read.aligned_string.size(); ++j)
                                  {
                                      int pos = j + read.start_idx;
                                      int ind = GENOME_STRING.find(read.aligned_string[j]);
                                      // this is before merging so degree is 1, but bettter
                                      // to be explicit
                                      af[pos][ind] += read.degree;
                                  }
                              }

                              std::vector<double> q = frequency_vector(af);
                              double divergence = 0;
                              for (size_t j = 0; j < p.size(); ++j)
                              {
                                  if (p[j] == 0)
                                      continue;
                                  double effective_q = std::max(q[j], 1e-10);

                                  divergence += p[j] * (std::log(p[j]) - std::log(effective_q));
                              }

                              {
                                tbb::queuing_mutex::scoped_lock my_lock{my_mutex};
                                if (divergence < best_score)
                                {
                                    best_score = divergence;
                                    best_set = current;
                                }
                              }
                          }
                      });

    this->aligned_reads = best_set;
    sort(this->aligned_reads.begin(), this->aligned_reads.end());
}

void sam::build() {
    collapsed_frequency_table = frequency_table;

    this->subsample();

    if (USE_READ_CORRECTION) {
        this->read_correction();
    }

    std::sort(this->aligned_reads.begin(), this->aligned_reads.end());


    if (USE_COLUMN_MERGING) {
        this->merge_duplicates();
    }
    else {
        for (const auto& col: aligned_reads) {
            std::string const col_name = col.degreeless_name();
            reverse_merge[col_name] = std::vector<std::string>{col.raw_name};
        }
    }
}

void sam::dump_reverse_merge(std::ostream &out) {
    for (const auto& field: reverse_merge) {
        out << field.first << '\t';
        for (const std::string& str: field.second) {
            out << str << '\t';
        }
        out << '\n';
    }
}


std::vector<raw_read> load_reads_from_proto(std::string const& reference, std::string const& filename, std::unordered_map<std::string, std::vector<std::string>> &reverse_merge) {
    timer t;

    Sam::sam data;

    boost::iostreams::filtering_istream instream;
    std::ifstream inpfile(filename, std::ios::in | std::ios::binary);
    if (!inpfile) {
        fprintf(stderr, "ERROR: Could not load the read protobuf from file: %s!\n", filename.c_str());
        exit(1);
    }
    instream.push(inpfile);
    google::protobuf::io::IstreamInputStream stream(&instream);
    google::protobuf::io::CodedInputStream input(&stream);
    
    //input.SetTotalBytesLimit(BIG_SIZE, BIG_SIZE);
    data.ParseFromCodedStream(&input);

    int read_count = data.reads_size();
    std::vector<raw_read> reads(read_count);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, read_count),
                      [&](tbb::blocked_range<size_t> range)
                      {
                          for (size_t i = range.begin(); i < range.end(); ++i)
                          {
                              raw_read &out = reads[i];
                              const auto &curr = data.reads()[i];
                              out.start = curr.start_idx();
                              out.end = curr.start_idx() + curr.content().size() - 1;
                              out.degree = curr.degree();
                              out.read = curr.read();
                              for (size_t i = 0; i < curr.content().size(); ++i)
                              {
                                  if (curr.content()[i] != reference[curr.start_idx() + i - 1])
                                  {
                                      mutation mut; 
                                      mut.ref = nuc_from_char(reference[curr.start_idx() + i - 1]);
                                      mut.mut = nuc_from_char(curr.content()[i]);
                                      mut.pos = out.start + i;

                                      out.mutations.push_back(std::move(mut));
                                  }
                              }
                          }
                      });

    /* finally, reverse merge table */
    int unmerged_count = 0;
    for (int i = 0; i < data.reverse_columns_size(); ++i) {
        Sam::column_info const &inv = data.reverse_columns()[i];
        for (int j = 0; j < inv.input_columns_size(); ++j) {
            reverse_merge[inv.column_name()].push_back(inv.input_columns()[j]);
        }
        unmerged_count += inv.input_columns().size();
    }

    printf("--- parsed %s containing %d merged / %d unmerged reads in %ld sec\n\n", filename.c_str(), read_count, unmerged_count, t.seconds());   
    return reads;
}