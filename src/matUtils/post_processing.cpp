#include "post_processing.hpp"

po::variables_map parse_post_processing_command(po::parsed_options parsed) {
    po::variables_map vm;
    po::options_description conv_desc("post_processing options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("read-vcf,v", po::value<std::string>()->default_value(""),
     "Output VCF file representing selected subtree. Default is full tree")
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

void post_processing(po::parsed_options parsed) {
    po::variables_map vm = parse_post_processing_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    std::string sample_vcf_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_samples.vcf";
    std::string hap_vcf_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_haplotypes.vcf";
    std::string hap_csv_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_hap_abundance.csv";
    std::string hap_clade_csv_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_hap_clade.csv";

    // Load input MAT and uncondense tree
    MAT::Tree T;
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    //Depth first expansion to get all nodes in the tree and 
    // comparison with given lineage to get all nodes of the required lineage 
    std::vector<MAT::Node*> dfs = T.depth_first_expansion(T.root); 
    
    //Checking how close are input samples with peaks
    std::unordered_map<int, struct read_info*> read_map;
    std::vector<std::string> vcf_samples;
    std::unordered_map<std::string, double> hap_abun_map;
    std::unordered_map<std::string, std::string> hap_clade_map;
    read_sample_vcf(vcf_samples, sample_vcf_filename);
    //read_vcf(T, dfs, read_map, hap_vcf_filename);
    read_csv(hap_abun_map, hap_csv_filename);
    read_csv(hap_clade_map, hap_clade_csv_filename);
    compute_abundance(T, read_map, hap_abun_map, hap_clade_map);
    compute_distance(T, dfs, read_map, vcf_samples);
}

void compute_distance(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, const std::unordered_map<int, struct read_info*> &read_map, const std::vector<std::string> &vcf_samples) {
    fprintf(stderr, "Haplotypes: %d\n\n", (int)read_map.size());
    printf("\nMUTATION DISTANCE NEW:\n");
    //Get Mutations of samples
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = get_mutations(T, sample);
        std::string best_node = "";
        auto itr = read_map.begin();
        while (itr != read_map.end()) {
            int curr_dist = mutation_distance(sample_mutations, itr->second->mutations);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = itr->second->read;
            }
            itr++;
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node.c_str(), min_dist);
    }

    //Samples' avg mut distance from other haplotypes of lineage
    printf("\n\nINTRA-LINEAGE AVG MUTATION DISTANCE:\n");
    for (auto sample: vcf_samples) {
        std::vector<MAT::Node*> hap_list;
        //Find all haplotypes of the lineage
        auto curr_clade = get_clade(T, T.get_node(sample));
        for (auto n: dfs) {
            if (n->clade_annotations[1] == curr_clade) {
                std::queue<MAT::Node*> remaining_nodes;
                remaining_nodes.push(n);
                while (remaining_nodes.size() > 0) {
                    MAT::Node* curr_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    if ((curr_node->clade_annotations[1] == "") || (curr_node->clade_annotations[1] == curr_clade)) {
                        if (curr_node->children.size() == 0)
                            hap_list.emplace_back(curr_node);
                        else {
                            for (auto c: curr_node->children) {
                                remaining_nodes.push(c);
                            }
                        }
                    }
                }
                break;
            }
        }
        //Compute Intra-lineage mutation distance
        double avg_distance = 0.0;
        auto sample_mutations = get_mutations(T, sample);
        for (auto hap: hap_list) {
            std::string hap_name = hap->identifier;
            if (hap_name != sample) {
                auto hap_mutations = get_mutations(T, hap_name);
                avg_distance += mutation_distance(sample_mutations, hap_mutations);
            }
        }
        avg_distance /= (int)hap_list.size() - 1;
        hap_list.clear();
        printf("Node: %s, avg_mutation_distance: %f\n", sample.c_str(), avg_distance);
    }
}

//Function to calculation distance between two nodes
int mutation_distance(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations) {
    std::sort(node1_mutations.begin(), node1_mutations.end(), compare_mutations);
    std::sort(node2_mutations.begin(), node2_mutations.end(), compare_mutations);
    auto n1_itr = node1_mutations.begin();
    while (n1_itr != node1_mutations.end()) {
        bool found_both = false;
        auto n2_itr = node2_mutations.begin();
        while (n2_itr != node2_mutations.end()) {
            if ((n2_itr->position == n1_itr->position) && (n2_itr->mut_nuc == n1_itr->mut_nuc)) {
                node2_mutations.erase(n2_itr);
                found_both = true;
                break;
            }
            else if (n2_itr->position == n1_itr->position) {
                node2_mutations.erase(n2_itr);
                break;
            }
            n2_itr++;
        }
        if (found_both)
            n1_itr = node1_mutations.erase(n1_itr);
        else
            n1_itr++;
    }
    return (int)(node1_mutations.size() + node2_mutations.size());
}

//Function to read csv containing abundance of haplotypes
void read_csv(std::unordered_map<std::string, double>& hap_abun_map, const std::string hap_csv_filename) {
    std::ifstream file(hap_csv_filename);
    if (!file.is_open()) {
        std::cout << "Failed to open the csv file" << std::endl;
        return;
    }

    std::string line;
    bool isFirstLine = true;  // To skip the header

    while (std::getline(file, line)) {
        if (isFirstLine) {
            isFirstLine = false;
            continue;  // Skip the header line
        }

        std::istringstream iss(line);
        std::string key, value;
        std::getline(iss, key, ',');
        std::getline(iss, value, ',');
        hap_abun_map[key] = std::stod(value);
    }
    file.close();
}

//Function to read csv containing lineage of haplotypes
void read_csv(std::unordered_map<std::string, std::string>& hap_clade_map, const std::string hap_clade_csv_filename) {
    std::ifstream file(hap_clade_csv_filename);
    if (!file.is_open()) {
        std::cout << "Failed to open the csv file" << std::endl;
        return;
    }

    std::string line;
    bool isFirstLine = true;  // To skip the header

    while (std::getline(file, line)) {
        if (isFirstLine) {
            isFirstLine = false;
            continue;  // Skip the header line
        }

        std::istringstream iss(line);
        std::string key, value;
        std::getline(iss, key, ',');
        std::getline(iss, value, ',');
        hap_clade_map[key] = value;
    }
    file.close();
}

//Computes Lineage Abundance
void compute_abundance(const MAT::Tree& T, const std::unordered_map<int, struct read_info*> &read_map, const std::unordered_map<std::string, double>& hap_abun_map, const std::unordered_map<std::string, std::string>& hap_clade_map) {
    std::unordered_map<std::string, double> clade_abun_map; 
    //Iterate through haplotype-abundance map
    auto h_abun_itr = hap_abun_map.begin();
    while (h_abun_itr != hap_abun_map.end()) {
        //Find the haplotype in haplotype-clade map 
        auto h_clade_itr = hap_clade_map.find(h_abun_itr->first);
        //Insert the clade in clade-abundance map if not present
        auto c_abun_itr = clade_abun_map.find(h_clade_itr->second);
        if (c_abun_itr == clade_abun_map.end())
            clade_abun_map.insert({h_clade_itr->second, h_abun_itr->second});
        //Else add the value to the map element
        else 
            clade_abun_map[h_clade_itr->second] += h_abun_itr->second;
        h_abun_itr++;
    }

    printf("\nLINEAGE ABUNDANCE (New Peaks)\n");
    auto c_abun_itr = clade_abun_map.begin();
    while (c_abun_itr != clade_abun_map.end()) {
        printf("%s: %f\n", c_abun_itr->first.c_str(), c_abun_itr->second);
        c_abun_itr++;
    }
}

//Getting all mutations of a haplotype
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
        else
            sm_itr++;
    }
    return sample_mutations;
}
