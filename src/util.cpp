#include <iostream>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include "panman/panmanUtils.hpp"

#include <tbb/parallel_sort.h>

#include "util.hpp"
#include "panman_bridge.hpp"

//Get the names of samples prsent in samples.vcf
std::vector<std::string> read_sample_vcf(const std::string& vcf_filename_samples) {
    std::vector<std::string> vcf_samples;

    // Boost library used to stream the contents of the input VCF file
    boost::filesystem::ifstream fileHandler(vcf_filename_samples);
    std::string s;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        std::stringstream ss(s);
        for (std::string tmp; ss >> tmp;) {
            words.push_back(tmp);
        }

        if (words.size() > 1) {
            //Checking for header
            if (words[1] == "POS") {
                //Leave certain fields based on our VCF format
                for (int j=9; j < (int)words.size(); j++)
                    vcf_samples.emplace_back(words[j]);
            }
        }
    }

    return vcf_samples;
}

boost::program_options::variables_map parseWBEcommand(boost::program_options::parsed_options parsed) {
    namespace po = boost::program_options;

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("Given Switch options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->default_value(""),
     "Input mutation-annotated tree file")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("output-files-prefix,v", po::value<std::string>()->default_value("my_vcf"),
    "Prefix to be used for dumping all intermediate files.")
    ("align-sam,s", po::value<std::string>()->default_value(""),
     "Input sam file representing reference sequence")
    ("distribution,d", po::value<std::string>()->default_value(""),
     "Give the distribution of samples, comma delimited.")
    ("haplotype-samples,w", po::value<int>()->default_value(10),
     "Give the number of haplotype samples")
    ("lineage,l", po::value<std::string>()->default_value(""),
     "Give lineage of samples, comma delimited.")
    ("prior-lineages,p", po::value<std::string>()->default_value(""),
     "Give lineage to be included, comma delimited.") 
    ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
    ("help,h", "Print help messages");
    // Collect all the unrecognized options from the first pass. This will include the
    // (positional) command name, so we need to erase that.
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    // Run the parser, with try/catch for help
    try {
        po::store(po::command_line_parser(opts)
                  .options(conv_desc)
                  .run(), vm);
        po::notify(vm);
    } catch(std::exception &e) {
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

std::vector<mutation> 
get_single_mutations(const std::string& ref, const panmanUtils::Node* node, const coord_converter &coord, bool ignore_root) {
    // if is root, then we treat it slightly different
    // because in panmat it's with respect to a blank sequence
    // but we do it with respect to a root sequence in the most case
    if (!node->parent && ignore_root) {
        return {};
    }

    std::vector<mutation> ret;
    std::unordered_set<size_t> covered;

    for (const panmanUtils::NucMut& mut: node->nucMutation) {
        size_t const count = (size_t) (mut.mutInfo >> 4);
        bool const is_deletion = (mut.mutInfo & 0xF) == panmanUtils::NucMutationType::ND || 
            (mut.mutInfo & 0xF) == panmanUtils::NucMutationType::NSNPD;

        for (size_t i = 0; i < count; ++i) {
            int offset_nuc_position = mut.nucGapPosition == -1 ? mut.nucPosition + i : mut.nucPosition;
            int offset_gap_position = mut.nucGapPosition == -1 ? -1 : mut.nucGapPosition + i; 
            int const start_pos = coord.query(mut.primaryBlockId, offset_nuc_position, offset_gap_position);

            mutation m;
            m.pos = start_pos;
            m.ref = nuc_from_char(ref[m.pos - 1]);
            int raw = (mut.nucs >> (4 * (5 - i))) & 0xF;
            m.mut = is_deletion ? NUC_GAP : nuc_from_pannuc(raw);
            ret.push_back(m);

            covered.insert(m.pos);
        }
    }

    for (const panmanUtils::BlockMut& block_mut : node->blockMutation) {
        assert(!block_mut.inversion);
        
        const auto [start, end] = coord.block_range(block_mut.primaryBlockId);
        if (block_mut.blockMutInfo) {
            // insertion
            for (size_t i = start; i <= end; ++i) {
                if (covered.find(i) != covered.end()) {
                    continue;
                }

                mutation mut;
                mut.pos = i;
                mut.ref = mut.mut = nuc_from_char(ref[i - 1]);
                ret.push_back(mut);

                covered.insert(i);
            }
        }
        else {
            // deletion
            for (size_t i = start; i <= end; ++i) {
                if (covered.find(i) != covered.end()) {
                    continue;
                }

                mutation mut;
                mut.pos = i;
                mut.ref = nuc_from_char(ref[i - 1]);
                mut.mut = NUC_GAP;
                ret.push_back(mut);

                covered.insert(i);
            }
        }
    }

    // if it's the root, whatever is not mentioned originally should be treated as gap
    // (if consensus is not already gap)
    if (node->parent == nullptr) {
        for (size_t i = 1; i <= ref.size(); ++i) {
            if (covered.find(i) == covered.end() && ref[i - 1] != '_') {
                mutation mut;
                mut.pos = i;
                mut.ref = nuc_from_char(ref[i - 1]);
                mut.mut = NUC_GAP;
                ret.push_back(mut);
            }
        }
    }

    return ret;
}

// we duplicate this logic a bunch...
int mutation_distance(std::vector<mutation> node1_mutations, std::vector<mutation> node2_mutations) {
    auto compareMutations = [](const mutation &a, const mutation &b)
    {
        if (a.pos != b.pos)
            return a.pos < b.pos;
        else
            return a.mut < b.mut;
    };

    tbb::parallel_sort(node1_mutations.begin(), node1_mutations.end(), compareMutations);
    tbb::parallel_sort(node2_mutations.begin(), node2_mutations.end(), compareMutations);

    int muts = 0;

    int i = 0;
    int last_i = node1_mutations.size();
    int j = 0;

    while (i < last_i || j < (int) node2_mutations.size()) {
        if (i == last_i) {                
            if (node2_mutations[j].mut != NUC_N) {
                ++muts;
            }
            ++j;
        }
        else if (j == (int) node2_mutations.size()) {
            if (node1_mutations[i].mut != NUC_N) {
                ++muts;
            }
            ++i;
        }
        else if (node1_mutations[i].pos < node2_mutations[j].pos) {
            if (node1_mutations[i].mut != NUC_N) {
                ++muts;
            }
            ++i;
        }
        else if (node1_mutations[i].pos > node2_mutations[j].pos) {
            if (node2_mutations[j].mut != NUC_N) {
                ++muts;
            }
            ++j;
        }
        else if (node1_mutations[i].pos == node2_mutations[j].pos && node1_mutations[i].mut != node2_mutations[j].mut && node1_mutations[i].mut != NUC_N && node2_mutations[j].mut != NUC_N) {
            ++muts;
            ++i; ++j;
        }
        else {
            ++i; ++j;
        }
    }

    return muts;
}
