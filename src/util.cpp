#include <iostream>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include "panman/panmanUtils.hpp"

#include <tbb/parallel_sort.h>

#include "util.hpp"
#include "panman_bridge.hpp"
#include "config.hpp"

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
        std::cerr << conv_desc << e.what() << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help"))
            exit(0);
        else
            exit(1);
    }
    return vm;
}

void
get_single_mutations(std::vector<mutation>& ret, const std::string& ref, const panmanUtils::Node* node, const coord_converter &coord) {
    std::unordered_set<size_t> covered;

    for (const panmanUtils::NucMut& mut: node->nucMutation) {
        size_t const count = (size_t) (mut.mutInfo >> 4);
        bool const is_deletion = (mut.mutInfo & 0xF) == panmanUtils::NucMutationType::ND || 
            (mut.mutInfo & 0xF) == panmanUtils::NucMutationType::NSNPD;
        for (size_t i = 0; i < count; ++i) {
            int offset_nuc_position = mut.nucGapPosition == -1 ? mut.nucPosition + i : mut.nucPosition;
            int offset_gap_position = mut.nucGapPosition == -1 ? -1 : mut.nucGapPosition + i; 
            int const pos = coord.query(mut.primaryBlockId, offset_nuc_position, offset_gap_position);

            mutation m;
            m.pos = pos;
            m.ref = nuc_from_char(ref[m.pos - 1]);
            int raw = (mut.nucs >> (4 * (5 - i))) & 0xF;
            m.mut = is_deletion ? NUC_GAP : nuc_from_pannuc(raw);
            if (!IGNORE_N_MUTS || m.mut != NUC_N) {
                ret.push_back(m);
                covered.insert(m.pos);
            }
        }
    }

    // blockMutation should only be present in root
    assert((node->parent == nullptr) == (node->blockMutation.size() != 0));

    return;
}

// precondition, both mutation lists are sorted
float mutation_distance(std::vector<mutation> const& node1_mutations, std::vector<mutation> const& node2_mutations) {
    int dels = 0, subs = 0, i = 0, j = 0;
    std::vector<int> deletions;

    while (i < (int) node1_mutations.size() || j < (int) node2_mutations.size()) {
        if (i == (int) node1_mutations.size()) {                
            if (node2_mutations[j].mut != NUC_N) {
                if (node2_mutations[j].mut == NUC_GAP)
                    deletions.emplace_back(node2_mutations[j].pos);
                    //++dels;
                else
                    ++subs;
            }
            ++j;
        }
        else if (j == (int) node2_mutations.size()) {
            if (node1_mutations[i].mut != NUC_N) {
                if (node1_mutations[i].mut == NUC_GAP)
                    deletions.emplace_back(node1_mutations[i].pos);
                    //++dels;
                else
                    ++subs;
            }
            ++i;
        }
        else if (node1_mutations[i].pos < node2_mutations[j].pos) {
            if (node1_mutations[i].mut != NUC_N) {
                if (node1_mutations[i].mut == NUC_GAP)
                    deletions.emplace_back(node1_mutations[i].pos);
                    //++dels;
                else
                    ++subs;
            }
            ++i;
        }
        else if (node1_mutations[i].pos > node2_mutations[j].pos) {
            if (node2_mutations[j].mut != NUC_N) {
                if (node2_mutations[j].mut == NUC_GAP)
                    deletions.emplace_back(node2_mutations[j].pos);
                    //++dels;
                else
                    ++subs;
            }
            ++j;
        }
        else if (node1_mutations[i].pos == node2_mutations[j].pos && node1_mutations[i].mut != node2_mutations[j].mut && node1_mutations[i].mut != NUC_N && node2_mutations[j].mut != NUC_N) {
            if (node1_mutations[i].mut == NUC_GAP || node2_mutations[j].mut == NUC_GAP)
                deletions.emplace_back(node1_mutations[i].pos);
                //++dels;
            else
                ++subs;

            ++i; ++j;
        }
        else {
            ++i; ++j;
        }
    }
    
    // Account Deletions
    if (deletions.size())
    {
        int start = deletions.front();
        for (size_t i = 0; i < deletions.size(); i++)
        {
            if (!i) 
            {
                start = deletions.front();
            }
            else if ((deletions[i] - deletions[i - 1]) > 1) 
            {
                dels += std::ceil((double)(deletions[i - 1] - start + 1) / 3);
                start = deletions[i];
            }
        }
        dels += std::ceil((double)(deletions.back() - start + 1) / 3);
    }

    return subs + (dels * DEL_SUBS_RATIO);
}