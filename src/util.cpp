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

namespace po = boost::program_options;
po::options_description conv_desc("Arguments");
void initializeConvDesc() {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message =
        "Number of threads to use when possible [DEFAULT uses all available cores, " +
        std::to_string(num_cores) + " detected on this machine]";

    conv_desc.add_options()
        ("input-mat,i", po::value<std::string>()->default_value(""), "Input mutation-annotated tree.")
        ("dataset,d", po::value<std::string>()->default_value(""), "Data folder containing reads.")
        ("max-reads,m", po::value<uint32_t>()->default_value(1e9), "Maximum number of reads.")
        ("file-prefix,p", po::value<std::string>()->default_value(""), "Prefix for intermediate files.")
        ("min-af,a", po::value<std::string>()->default_value("0.005"), "Allele Frequency threshold for masking errorneous alleles.")
        ("min-phred,q", po::value<u_int32_t>()->default_value(20), "Phred Score threshold for masking low quality alleles.")
        ("threads,T", po::value<uint32_t>()->default_value(num_cores), num_threads_message.c_str())
        ("help,h", "Print help messages");
}

po::variables_map parseWEPPcommand(po::parsed_options parsed) {
    initializeConvDesc();

    po::variables_map vm;
    std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
    opts.erase(opts.begin());

    try {
        po::store(po::command_line_parser(opts).options(conv_desc).run(), vm);
        po::notify(vm);
        if (vm.count("help")) {
            std::cout << conv_desc << std::endl;
            exit(0);
        }
    } catch (std::exception &e) {
        std::cerr << conv_desc << std::endl;
        if (vm.count("help")) {
            std::cout << conv_desc << std::endl;
            exit(0);
        }
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

    while (i < (int) node1_mutations.size() || j < (int) node2_mutations.size()) {
        if (i == (int) node1_mutations.size()) {                
            if (node2_mutations[j].mut != NUC_N) {
                if (node2_mutations[j].mut == NUC_GAP)
                    ++dels;
                else
                    ++subs;
            }
            ++j;
        }
        else if (j == (int) node2_mutations.size()) {
            if (node1_mutations[i].mut != NUC_N) {
                if (node1_mutations[i].mut == NUC_GAP)
                    ++dels;
                else
                    ++subs;
            }
            ++i;
        }
        else if (node1_mutations[i].pos < node2_mutations[j].pos) {
            if (node1_mutations[i].mut != NUC_N) {
                if (node1_mutations[i].mut == NUC_GAP)
                    ++dels;
                else
                    ++subs;
            }
            ++i;
        }
        else if (node1_mutations[i].pos > node2_mutations[j].pos) {
            if (node2_mutations[j].mut != NUC_N) {
                if (node2_mutations[j].mut == NUC_GAP)
                    ++dels;
                else
                    ++subs;
            }
            ++j;
        }
        else if (node1_mutations[i].pos == node2_mutations[j].pos && node1_mutations[i].mut != node2_mutations[j].mut && node1_mutations[i].mut != NUC_N && node2_mutations[j].mut != NUC_N) {
            if (node1_mutations[i].mut == NUC_GAP || node2_mutations[j].mut == NUC_GAP)
                ++dels;
            else
                ++subs;

            ++i; ++j;
        }
        else {
            ++i; ++j;
        }
    }
    
    return subs + (dels * DEL_SUBS_RATIO);
}

std::vector<mutation> mutation_distance_vector(std::vector<mutation> const& node1_mutations, std::vector<mutation> const& node2_mutations) {
    int i = 0, j = 0;
    std::vector<mutation> mutations;

    while (i < (int) node1_mutations.size() || j < (int) node2_mutations.size()) {
        if (i == (int) node1_mutations.size()) {                
            if (node2_mutations[j].mut != NUC_N) {
                mutations.emplace_back(node2_mutations[j]);
            }
            ++j;
        }
        else if (j == (int) node2_mutations.size()) {
            if (node1_mutations[i].mut != NUC_N) {
                mutations.emplace_back(node1_mutations[i]);
            }
            ++i;
        }
        else if (node1_mutations[i].pos < node2_mutations[j].pos) {
            if (node1_mutations[i].mut != NUC_N) {
                mutations.emplace_back(node1_mutations[i]);
            }
            ++i;
        }
        else if (node1_mutations[i].pos > node2_mutations[j].pos) {
            if (node2_mutations[j].mut != NUC_N) {
                mutations.emplace_back(node2_mutations[j]);
            }
            ++j;
        }
        else if (node1_mutations[i].pos == node2_mutations[j].pos && node1_mutations[i].mut != node2_mutations[j].mut && node1_mutations[i].mut != NUC_N && node2_mutations[j].mut != NUC_N) {
            if (node1_mutations[i].mut == NUC_GAP)
                mutations.emplace_back(node1_mutations[i]);
            else if (node2_mutations[j].mut == NUC_GAP)
                mutations.emplace_back(node2_mutations[j]);
            else
                mutations.emplace_back(node1_mutations[i]);

            ++i; ++j;
        }
        else {
            ++i; ++j;
        }
    }
    
    return mutations;
}