#include <iostream>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include "panman/panman.hpp"

#include <tbb/parallel_sort.h>

#include "util.hpp"

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

panmanUtils::Tree create_condensed_tree(panmanUtils::Node* ref_root, const std::vector<raw_read> &read_map, std::unordered_map<panmanUtils::Node*, std::vector<panmanUtils::Node*>> &node_mappings) {
    //MAP of site_read
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map[i];
        for (int j = rp.start; j <= rp.end; j++)
            site_read_map.insert(j);
    }

    //REMOVE sites not covered by reads
    std::queue<std::pair<panmanUtils::Node*, panmanUtils::Node*>> remaining_nodes;

    panmanUtils::Tree T;
    auto new_node = T.create_node("DUMMY-CONDENSED", -1.0, ref_root->clade_annotations.size());
    node_mappings[new_node] = std::vector<panmanUtils::Node*>();
    remaining_nodes.push(std::pair<panmanUtils::Node*, panmanUtils::Node*>(ref_root, new_node));
    
    //Add the new_node to the node_mappings
    while(remaining_nodes.size() > 0) {
        auto r_curr_node = remaining_nodes.front().first;
        auto n_parent_node = remaining_nodes.front().second;
        remaining_nodes.pop();
        std::vector<MAT::Mutation> covered_mutations;
        for (const auto& mut: r_curr_node->mutations) {
            if (site_read_map.find(mut.position) != site_read_map.end())
                covered_mutations.emplace_back(mut);
        }

        //Add a new node to tree if mutation found in range
        if (!covered_mutations.empty()) {
            //Create a new_node
            auto new_node = T.create_node(r_curr_node->identifier, n_parent_node, -1.0);
            //Add mutations to new_node
            new_node->mutations.reserve(covered_mutations.size());
            new_node->mutations.insert(new_node->mutations.end(), covered_mutations.begin(), covered_mutations.end());
            covered_mutations.clear();
            //Add new_node to the node_mappings
            node_mappings[new_node] = {r_curr_node};
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<panmanUtils::Node*, panmanUtils::Node*>(child, new_node));
        }
        //Without any mutation, can only be Placed as a child if r_curr_node is NOT leaf node
        else {
            //Add current_node to the n_parent_node's list in node_mappings
            node_mappings[n_parent_node].emplace_back(r_curr_node);
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<panmanUtils::Node*, panmanUtils::Node*>(child, n_parent_node));
        }
    }

    return T;
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
    ("ref-fasta,f", po::value<std::string>()->default_value(""),
     "Input fasta file representing reference sequence")
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
    auto n1_iterator = node1_mutations.begin();
    auto n2_iterator = node2_mutations.begin();
    int distance = 0;
    while (n1_iterator != node1_mutations.end() && n2_iterator != node2_mutations.end()) {
        if (n1_iterator->pos == n2_iterator->pos) {
            if (n1_iterator->mut != n2_iterator->mut) {
                distance++;
            }
            n1_iterator++;
            n2_iterator++;
        } else if (n1_iterator->pos < n2_iterator->pos) {
            distance++;
            n1_iterator++;
        } else {
            distance++;
            n2_iterator++;
        }
    }

    distance += std::distance(n1_iterator, node1_mutations.end());
    distance += std::distance(n2_iterator, node2_mutations.end());

    return distance;
}

std::vector<MAT::Mutation> get_mutations(const MAT::Tree& T, const std::string sample) {
    std::vector<MAT::Mutation> sample_mutations;
    for (auto anc: T.rsearch(sample, true)) { //Checking all ancestors of a node
        for (auto mut: anc->mutations) {
            auto sm_itr = sample_mutations.begin();
            while (sm_itr != sample_mutations.end()) {
                if (sm_itr->position == mut.position)
                    break;
                sm_itr++;
            }
            if (sm_itr == sample_mutations.end())
                sample_mutations.emplace_back(mut);
        }
    }
    //Remove Back-Mutations
    auto sm_itr = sample_mutations.begin();
    while (sm_itr != sample_mutations.end()) {
        if (sm_itr->ref_nuc == sm_itr->mut_nuc)
            sm_itr = sample_mutations.erase(sm_itr);
        else {
            sm_itr->par_nuc = sm_itr->ref_nuc; 
            sm_itr++;
        }
    }
    return sample_mutations;
}
