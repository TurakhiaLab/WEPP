#include "wbe.hpp"

// freyja - includes indels, litreally every single possible mutation of a chunk along with frequency
// freyja-depth - how many reads/chunks were there starting at a given nucleotide
// vcf - for all non indels with a mutation, look at every single read and tell which mutation it corresponds to
constexpr bool USE_READ_CORRECTION = true;
constexpr bool USE_COLUMN_MERGING  = true;
constexpr double frequency_read_cutoff = 0.005;
const std::string CHROM = "NC_045512v2";

void sam2VCF(po::parsed_options parsed) {
    //main argument for the complex extract command
    po::variables_map vm = parseWBEcommand(parsed);
    std::string dir_prefix = vm["output-directory"].as<std::string>();

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string freyja_vcf_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads_freyja.vcf";
    std::string freyja_depth_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads_freyja.depth";
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    std::string sam_file = dir_prefix + vm["align-sam"].as<std::string>();

    //Reading reference genome
    timer.Start();
    std::ifstream fasta_f(ref_fasta);
    if (!fasta_f.is_open()) {
        std::cerr << "Error: Unable to open file " << ref_fasta << std::endl;
        exit(1); 
    }
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }

    std::ofstream outfile_vcf(vcf_filename_reads, std::ios::out | std::ios::binary);
    std::ofstream outfile_freyja_vcf(freyja_vcf_file, std::ios::out | std::ios::binary);
    std::ofstream outfile_freyja_depth(freyja_depth_file, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_vcf, outbuf_freyja_vcf, outbuf_freyja_depth;
    
    if (vcf_filename_reads.find(".gz\0") != std::string::npos)
        outbuf_vcf.push(boost::iostreams::gzip_compressor());
    if (freyja_vcf_file.find(".gz\0") != std::string::npos)
        outbuf_freyja_vcf.push(boost::iostreams::gzip_compressor());
    if (freyja_depth_file.find(".gz\0") != std::string::npos)
        outbuf_freyja_depth.push(boost::iostreams::gzip_compressor());
    
    outbuf_vcf.push(outfile_vcf);
    outbuf_freyja_vcf.push(outfile_freyja_vcf);
    outbuf_freyja_depth.push(outfile_freyja_depth);
    std::ostream vcf(&outbuf_vcf);
    std::ostream freyja_vcf(&outbuf_freyja_vcf);
    std::ostream freyja_depth(&outbuf_freyja_depth);

    SAM sam{ref_seq};
    boost::filesystem::ifstream fileHandler(sam_file);
    std::string s;
    while (getline(fileHandler, s))
        sam.add_read(s);
    sam.build();
    sam.dump_vcf(vcf);
    sam.dump_freyja(freyja_depth, freyja_vcf);
    fprintf(stderr, "\nFiles generated in %ld sec \n\n", (timer.Stop() / 1000));
}

SAM::SAM(const std::string& ref) : reference_seq{ref} {
    this->frequency_table.resize(ref.size());
    this->indel_frequency_table.resize(ref.size());
    this->collapsed_frequency_table.resize(ref.size());
}

void load_reads_from_proto(std::string const& filename, std::unordered_map<size_t, struct read_info *>& reads, std::unordered_map<std::string, std::vector<std::string>> &reverse_merge) {
    Sam::sam data;

    boost::iostreams::filtering_istream instream;
    std::ifstream inpfile(filename, std::ios::in | std::ios::binary);
    if (filename.find(".gz\0") != std::string::npos) {
        if (!inpfile) {
            fprintf(stderr, "ERROR: Could not load the mutation-annotated tree object from file: %s!\n", filename.c_str());
            exit(1);
        }
        try {
            instream.push(boost::iostreams::gzip_decompressor());
            instream.push(inpfile);
        } catch(const boost::iostreams::gzip_error& e) {
            std::cout << e.what() << '\n';
        }
    } else {
        instream.push(inpfile);
    }
    google::protobuf::io::IstreamInputStream stream(&instream);
    google::protobuf::io::CodedInputStream input(&stream);
//    input.SetTotalBytesLimit(BIG_SIZE, BIG_SIZE);
    data.ParseFromCodedStream(&input);

    int read_count = data.reads_size();
    for (int i = 0; i < read_count; ++i) {
        struct read_info *out = new read_info{};
        const auto& curr = data.reads()[i];
        out->start = curr.start_idx();
        out->end = curr.end_idx();
        out->degree = curr.degree();
        out->read = curr.read();
        for (int j = 0; j < curr.mutations_size(); ++j) {
            const auto& mut = curr.mutations()[j];
            MAT::Mutation mutation;
            mutation.is_missing = mut.is_missing();
            mutation.chrom = mut.chrom();
            mutation.mut_nuc = mut.mut_nuc();
            mutation.ref_nuc = mut.ref_nuc();
            mutation.par_nuc = mut.par_nuc();
            mutation.position = mut.position();
            out->mutations.push_back(std::move(mutation));
        }

        reads[i] = out;
    }

    /* finally, reverse merge table */
    for (int i = 0; i < data.reverse_columns_size(); ++i) {
        Sam::column_info const &inv = data.reverse_columns()[i];
        for (int j = 0; j < inv.input_columns_size(); ++j) {
            reverse_merge[inv.column_name()].push_back(inv.input_columns()[j]);
        }
    }
}

void SAM::add_read(const std::string& line) {
    std::vector<std::string> tokens;
    MAT::string_split(line, tokens);

    /* skip if unmapped */
    if (tokens.empty() || tokens[0][0] == '@' || std::stoi(tokens[1]) & 4)
    {
        return;
    }

    /* zero indexed */
    int start_idx = std::stoi(tokens[3]);
    const std::string& read_seq = tokens[9];

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
    for (const auto& [cig_len, cig_val]: cigar_chunks) {
        int curr_sub_len = 0;
        while (curr_sub_len < cig_len)
        {
            int nuc_pos = start_idx + seq_idx - ins_count + del_count;
            switch (cig_val)
            {
            case 'I':
                // Insertion at pos 1 requires ref_nuc at pos 1 after the alt_nuc
                if (nuc_pos == 1)
                {
                    ref_nuc = reference_seq[0];
                    alt_nuc = read_seq.substr(seq_idx, cig_len) + ref_nuc;
                }
                // Otherwise the ref_nuc needs to be from prev pos
                else
                {
                    alt_nuc = ref_nuc + read_seq.substr(seq_idx, cig_len);
                    nuc_pos -= 1;
                }
                curr_sub_len += cig_len;
                seq_idx += cig_len;
                ins_count += cig_len;
                /* build stays the same*/
                break;

            case 'D':
                // Deletion at pos 1 requires alt_nuc at pos 1 after the ref_nuc
                if (nuc_pos == 1)
                {
                    alt_nuc = reference_seq[cig_len];
                    ref_nuc = reference_seq.substr(0, cig_len) + alt_nuc;
                }
                else
                {
                    alt_nuc = ref_nuc;
                    ref_nuc += reference_seq.substr(nuc_pos - 1, cig_len);
                    nuc_pos -= 1;
                }
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

            default:
                ref_nuc = reference_seq[nuc_pos - 1];
                alt_nuc = std::string(1, read_seq[seq_idx]);
                build += alt_nuc;
                curr_sub_len++;
                seq_idx++;
            }
            /* update frequency tables */
            if (ref_nuc.size() == 1 && alt_nuc.size() == 1) {
                int sub = (int) GENOME_STRING.find(alt_nuc[0]);
                frequency_table[nuc_pos - 1][sub]++;
            }
            else {
                indel_frequency_table[nuc_pos - 1][std::make_pair(ref_nuc.size(), alt_nuc)]++; 
            }
        }
    }

    /* add to column name list */
    SAM_read read{tokens[0], start_idx - 1, 1, build};
    aligned_reads.push_back(std::move(read));
}

void SAM::read_correction() {
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
            char effective = align[j] == '_' ? reference_seq[j] : align[j];
            int curr = GENOME_STRING.find(effective);
            int indx = start + j;

            if (frequency_read_cutoff - (double) collapsed_frequency_table[indx][curr] / total_occurences[indx] > 1e-9) {
                align[j] = GENOME_STRING[majority[indx]];
            }
        }
    }

    /* update frequency table */
    for (int i = 0; i < (int) reference_seq.size(); ++i) {
        for (int j = 0; j < (int) GENOME_STRING.size(); ++j) {
            if (j == majority[i]) continue;

            if (frequency_read_cutoff - (double) collapsed_frequency_table[i][j] / total_occurences[i] > 1e-9) {
                collapsed_frequency_table[i][majority[i]] += std::exchange(collapsed_frequency_table[i][j], 0);
            }
        }
    }
}

void SAM::merge_duplicates() {
    assert(!aligned_reads.empty());

    std::vector<SAM_read> merged{aligned_reads[0]};

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
        std::string generated = merge.raw_name + 
            "_READ_" + std::to_string(merge.start_idx + 1) + 
            "_" + std::to_string(merge.start_idx + 1 + merge.aligned_string.size() - 1) + 
            "_" + std::to_string(merge.degree);

        for (int j = 0; j < merge.degree; ++j) {
            this->reverse_merge[generated].push_back(aligned_reads[j + reverse].raw_name);
        }
        reverse += merge.degree;
    }

    this->aligned_reads = merged;
}

void SAM::build()
{
    /* create collapsed frequencies, treating indels as identity transformations */
    for (size_t i = 0; i < indel_frequency_table.size(); ++i) {
        collapsed_frequency_table[i] = frequency_table[i];
        for (const auto& freq: indel_frequency_table[i]) {
            int ref = (int) GENOME_STRING.find(reference_seq[i]);
            collapsed_frequency_table[i][ref] += freq.second;
        }
    }

    if (USE_READ_CORRECTION) {
        this->read_correction();
    }

    std::sort(this->aligned_reads.begin(), this->aligned_reads.end());

    if (USE_COLUMN_MERGING) {
        this->merge_duplicates();
    }
    else {
        for (const auto& col: aligned_reads) {
            std::string const col_name = col.raw_name + "_READ_" 
                + std::to_string(col.start_idx + 1) + "_"
                + std::to_string(col.start_idx + 1 + col.aligned_string.size() - 1);
            reverse_merge[col_name] = std::vector<std::string>{col.raw_name};
        }
    }
}

void SAM::dump_reverse_merge(std::ostream &out) {
    for (const auto& field: reverse_merge) {
        out << field.first << '\t';
        for (const std::string& str: field.second) {
            out << str << '\t';
        }
        out << '\n';
    }
}

void SAM::dump_vcf(std::ostream& out) {
    out << "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    for (const auto& read: this->aligned_reads) {
        out << '\t' << read.raw_name << "_READ_" << read.start_idx + 1 << "_" << read.start_idx + 1 + read.aligned_string.size() - 1;
        if (USE_COLUMN_MERGING) {
            out << '_' << read.degree;
        }
    }

    for (int i = 0; i < (int) reference_seq.size(); ++i) {
        int first_index = -1;
        for (int j = 0; j < (int) GENOME_STRING.size(); ++j) {
            if (GENOME_STRING[j] != reference_seq[i] && collapsed_frequency_table[i][j]) {
                first_index = j;
                break;
            }
        }

        /* no substitutions */
        if (first_index == -1) {
            continue;
        }

        out << '\n' << CHROM << '\t' << i + 1 << "\t";

        std::string muts; 
        for (int j = first_index; j < (int) GENOME_STRING.size(); ++j) {
            if (GENOME_STRING[j] != reference_seq[i] && collapsed_frequency_table[i][j]) {  
                if (j > first_index) {
                    out << ',';
                }
                out << reference_seq[i] << i + 1 << GENOME_STRING[j];
                muts.push_back(GENOME_STRING[j]);
            }
        }
        /* alternate*/
        out << "\t" << reference_seq[i] << "\t";
         for (char c: muts) {
            if (c != muts[0]) out << ',';
            out << c;
        }
        out << "\t.\t.\t.\t.";
        for (int j = 0; j < (int) aligned_reads.size(); ++j) {
            int one_index = 0;
            int s = aligned_reads[j].start_idx;
            int e = aligned_reads[j].start_idx + aligned_reads[j].aligned_string.size() - 1;

            if (i < s || i > e) {
                one_index = 0;
            }
            else {
                char c = aligned_reads[j].aligned_string[i - s];
                size_t ind;
                /* if it had a mutation */
                if ((ind = muts.find(c)) != std::string::npos)
                {
                    one_index = 1 + (int)ind;
                }
            }

            out << '\t' << one_index;
        }
    }
}

void SAM::dump_freyja(std::ostream& dout, std::ostream& vout) {
    vout << "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    for (int i = 0; i < (int) reference_seq.size(); ++i) {
        int total_count = 0;

        bool found = false;
        for (const auto& indel: indel_frequency_table[i]) {
            found = true;
            total_count += indel.second;
        }

        for (int j = 0; j < (int) GENOME_STRING.size(); ++j) {
            if (GENOME_STRING[j] != reference_seq[i] && frequency_table[i][j]) {
                found = true;
            }
            total_count += frequency_table[i][j];
        }

        dout << "NC_045512v2\t" << i + 1 << '\t' << reference_seq[i] << '\t' << total_count << '\n';
        if (found) {
            /* mutation present */
            for (const auto& indel: indel_frequency_table[i]) {
                vout << "\nNC_045512v2\t" << 
                        i + 1 << '\t' << 
                        reference_seq.substr(i, indel.first.first) << i + 1 << indel.first.second << 
                        '\t' << reference_seq.substr(i, indel.first.first) << '\t' << indel.first.second <<
                        "\t.\t.\tAF=" << 
                        std::to_string((double) indel.second / (double) total_count);
            }

            for (int j = 0; j < (int)GENOME_STRING.size(); ++j)
            {
                if (GENOME_STRING[j] != reference_seq[i] && frequency_table[i][j])
                {
                    vout << "\nNC_045512v2\t" << 
                        i + 1 << '\t' << reference_seq[i] << i + 1 << GENOME_STRING[j] << 
                        '\t' << reference_seq[i] << '\t' << GENOME_STRING[j] << "\t.\t.\tAF=" << std::to_string((double)frequency_table[i][j] / (double)total_count);
                }
            }
        }
    }
}
