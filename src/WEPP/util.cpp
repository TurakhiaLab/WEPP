#include <iostream>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include "src/usher_graph.hpp"

#include "util.hpp"

//Get the names of samples prsent in samples.vcf
std::vector<std::pair<std::string, std::vector<MAT::Mutation>>> read_sample_vcf(const std::string& vcf_filename_samples) {
    std::vector<std::pair<std::string, std::vector<MAT::Mutation>>> vcf_samples;
    // Boost library used to stream the contents of the input VCF file
    boost::filesystem::ifstream fileHandler(vcf_filename_samples);
    std::string s;
    bool header_found = false;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            //Checking for header
            if (words[1] == "POS") {
                header_found = true;
                //Leave certain fields based on our VCF format
                for (int j = 9; j < (int)words.size(); j++)
                    vcf_samples.emplace_back(std::make_pair(words[j], std::vector<MAT::Mutation>{}));
            }
            else if (header_found) {
                std::vector<std::string> alleles;
                //Checking for different alleles at a site
                MAT::string_split(words[4], ',', alleles);
                for (int j = 9; j < (int)words.size(); j++) {
                    int idx = j - 9;
                    MAT::Mutation m;
                    m.chrom = words[0];
                    m.position = std::stoi(words[1]);
                    //Checking the mutating allele value within the allele sizes
                    if (std::stoi(words[j]) > int(alleles.size())) {
                        fprintf(stderr, "\n\nPosition: %d, idx = %d,\n", m.position, idx);
                        fprintf(stderr, "Allele_id: %d, Alleles_size: %ld\n\n",std::stoi(words[j]), alleles.size());
                    }
                    m.ref_nuc = MAT::get_nuc_id(words[3][0]);
                    assert((m.ref_nuc & (m.ref_nuc-1)) == 0); //check if it is power of 2
                    m.par_nuc = m.ref_nuc;
                    // Alleles such as '.' should be treated as missing
                    // data. if the word is numeric, it is an index to one
                    // of the alleles
                    if (isdigit(words[j][0])) {
                        int allele_id = std::stoi(words[j]);
                        if (allele_id > 0) {
                            std::string allele = alleles[allele_id-1];
                            if (allele[0] == 'N') {
                                m.is_missing = true;
                                m.mut_nuc = MAT::get_nuc_id('N');
                            } 
                            else {
                                auto nuc = MAT::get_nuc_id(allele[0]);
                                if (nuc == MAT::get_nuc_id('N')) {
                                    m.is_missing = true;
                                } 
                                else {
                                    m.is_missing = false;
                                }
                                m.mut_nuc = nuc;
                            }
                            vcf_samples[idx].second.emplace_back(m);
                        }
                    } 
                    else {
                        m.is_missing = true;
                        m.mut_nuc = MAT::get_nuc_id('N');
                        vcf_samples[idx].second.emplace_back(m);
                    }
                }
            }
        }
    }

    return vcf_samples;
}

MAT::Tree create_condensed_tree(MAT::Node* ref_root, const std::unordered_set<int>& site_read_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings) {
    std::queue<std::pair<MAT::Node*, MAT::Node*>> remaining_nodes;
    MAT::Tree T;

    // Create ref_root as root of new tree
    auto new_node = T.create_node(ref_root->identifier, -1.0, ref_root->clade_annotations.size());
    //Add mutations to new_node
    for (const auto& mut: ref_root->mutations) {
        if (site_read_map.find(mut.position) != site_read_map.end())
            new_node->mutations.emplace_back(mut);
    }
    //Add new_node to the node_mappings
    node_mappings[new_node] = {ref_root};
    //Add children to remaining_nodes    
    for (auto child: ref_root->children)
        remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, new_node));

    
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
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, new_node));
        }
        //Without any mutation, can only be Placed as a child if r_curr_node is NOT leaf node
        else {
            //Add current_node to the n_parent_node's list in node_mappings
            node_mappings[n_parent_node].emplace_back(r_curr_node);
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, n_parent_node));
        }
    }

    return T;
}

boost::program_options::variables_map parseWEPPcommand(boost::program_options::parsed_options parsed) {
    namespace po = boost::program_options;

    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("WEPP Arguments");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->default_value(""),
     "Input mutation-annotated tree file.")
    ("dataset,d", po::value<std::string>()->default_value(""),
      "Data folder containing reads.")
    ("max-reads,m", po::value<uint32_t>()->default_value(1e9),
     "Maximum number of reads to use.")
    ("file-prefix,p", po::value<std::string>()->default_value(""),
    "Prefix to be used for dumping all intermediate files.")
    ("ref-fasta,f", po::value<std::string>()->default_value(""),
     "Input fasta file representing reference sequence.")
    ("min-af,a", po::value<std::string>()->default_value("0.005"),
     "Allele Frequency threshold for masking errorneous alleles: Illumina: 0.005, Ion Torrent: 0.015, ONT: 0.02.")
    ("min-phred,q", po::value<u_int32_t>()->default_value(20),
     "Phred Score threshold for masking low quality alleles.")
    ("clade-idx,c", po::value<u_int32_t>()->default_value(1),
     "Index used for inferring lineage proportions from haplotypes.")
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
        // If help is requested, show it and exit
        if (vm.count("help")) {
            std::cout << conv_desc << std::endl;
            exit(0);
        }
    } catch(std::exception &e) {
        std::cerr << conv_desc << std::endl;
        // Return with error code 1 unless the user specifies help
        if (vm.count("help")) {
            std::cout << conv_desc << std::endl;
            exit(0);
        }
        else
            exit(1);
    }
    return vm;
}

int mutation_distance(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations) {
    auto compareMutations = [](const MAT::Mutation &a, const MAT::Mutation &b)
    {
        if (a.position != b.position)
            return a.position < b.position;
        else
            return a.mut_nuc < b.mut_nuc;
    };

    tbb::parallel_sort(node1_mutations.begin(), node1_mutations.end(), compareMutations);
    tbb::parallel_sort(node2_mutations.begin(), node2_mutations.end(), compareMutations);
    auto n1_iterator = node1_mutations.begin();
    auto n2_iterator = node2_mutations.begin();
    int distance = 0;
    while (n1_iterator != node1_mutations.end() && n2_iterator != node2_mutations.end()) {
        if (n1_iterator->position == n2_iterator->position) {
            if (n1_iterator->mut_nuc != n2_iterator->mut_nuc) {
                distance++;
            }
            n1_iterator++;
            n2_iterator++;
        } else if (n1_iterator->position < n2_iterator->position) {
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

std::vector<MAT::Mutation> mutation_distance_vector(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations, const std::string& reference) {
    auto compareMutations = [](const MAT::Mutation &a, const MAT::Mutation &b)
    {
        if (a.position != b.position)
            return a.position < b.position;
        else
            return a.mut_nuc < b.mut_nuc;
    };

    tbb::parallel_sort(node1_mutations.begin(), node1_mutations.end(), compareMutations);
    tbb::parallel_sort(node2_mutations.begin(), node2_mutations.end(), compareMutations);
    auto n1_iterator = node1_mutations.begin();
    auto n2_iterator = node2_mutations.begin();
    std::vector<MAT::Mutation> mut_diff;
    while (n1_iterator != node1_mutations.end() && n2_iterator != node2_mutations.end()) {
        if (n1_iterator->position == n2_iterator->position) {
            if (n1_iterator->mut_nuc != n2_iterator->mut_nuc) {
                mut_diff.emplace_back(*n1_iterator);
            }
            n1_iterator++;
            n2_iterator++;
        } else if (n1_iterator->position < n2_iterator->position) {
            mut_diff.emplace_back(*n1_iterator);
            n1_iterator++;
        } else {
            auto mut = *n2_iterator;
            mut.mut_nuc = MAT::get_nuc_id(reference[mut.position - 1]);
            mut_diff.emplace_back(mut);
            n2_iterator++;
        }
    }

    while (n1_iterator != node1_mutations.end()) {
        mut_diff.emplace_back(*n1_iterator);
        n1_iterator++;
    }
    
    while (n2_iterator != node2_mutations.end()) {
        auto mut = *n2_iterator;
        mut.mut_nuc = MAT::get_nuc_id(reference[mut.position - 1]);
        mut_diff.emplace_back(mut);
        n2_iterator++;
    }

    return mut_diff;
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

size_t get_num_leaves(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, MAT::Node* condensed_node) {
    size_t leaves = 0;
    if (!condensed_node_mappings.find(condensed_node)->second.empty()) {
        auto front_node = condensed_node_mappings.find(condensed_node)->second.front();
        std::queue<MAT::Node*> remaining_nodes;
        remaining_nodes.push(front_node);
        while (remaining_nodes.size() > 0) {
            MAT::Node* curr_node = remaining_nodes.front();
            remaining_nodes.pop();
            if (curr_node->children.empty())
                leaves++;
            else {
                for (auto c: curr_node->children)
                    remaining_nodes.push(c);
            }
        }
    }
    return leaves;
}