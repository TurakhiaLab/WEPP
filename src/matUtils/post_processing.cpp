#include "post_processing.hpp"

po::variables_map parse_post_processing_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("post_processing options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("read-vcf,v", po::value<std::string>()->default_value(""),
     "Output VCF file representing selected subtree. Default is full tree")
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

void post_processing(po::parsed_options parsed) {
    po::variables_map vm = parse_post_processing_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    std::string sample_vcf_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_samples.vcf";
    std::string hap_vcf_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_haplotypes.vcf";
    std::string read_vcf_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_reads.vcf";
    std::string hap_csv_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_hap_abundance.csv";
    std::string hap_clade_csv_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_hap_clade.csv";
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    //Initialize TBB threads
    tbb::task_scheduler_init init(num_threads);
    
    // Load input MAT and uncondense tree
    MAT::Tree T;
    T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();

    //////////////////////////////////////
    //std::unordered_map<int, struct read_info*> read_map, hap_map;
    
    //std::vector<std::string> vcf_samples;
    //read_sample_vcf(vcf_samples, sample_vcf_filename);
    //for (auto sample: vcf_samples) {
    //    auto sample_mutations = get_mutations(T, sample);
    //    auto clade = get_clade(T, T.get_node(sample));
    //    printf("TREE %s %s -> %d\n", clade.c_str(), sample.c_str(), (int)sample_mutations.size());
    //}
    
    //read_vcf(hap_map, sample_vcf_filename);
    //for (auto hap: hap_map)
    //    printf("VCF %s -> %d\n", hap.second->read.c_str(), (int)hap.second->mutations.size());

    ////std::unordered_map<int, int> read_mut_count;
    ////read_vcf(read_map, read_vcf_filename);
    ////for (auto read: read_map) {
    ////    int mut_count = read.second->mutations.size();
    ////    auto rm_itr = read_mut_count.find(mut_count);
    ////    if (rm_itr == read_mut_count.end())
    ////        read_mut_count.insert({mut_count, 1});
    ////    else
    ////        read_mut_count[mut_count] += 1;
    ////}

    ////for (auto rm_count: read_mut_count) {
    ////    printf("%d -> %d\n", rm_count.first, rm_count.second);
    ////}
    ////////////////////////////////////

    //Checking how close are input samples with peaks
    std::unordered_map<size_t, struct read_info*> read_map, hap_map;
    tbb::concurrent_hash_map<std::string, std::vector<size_t>> hap_read_map;
    std::vector<std::string> vcf_samples;
    std::unordered_map<std::string, double> hap_abun_map;
    read_sample_vcf(vcf_samples, sample_vcf_filename);
    read_vcf(hap_map, hap_vcf_filename);
    read_csv(hap_abun_map, hap_csv_filename);
    
    ////read_vcf(read_map, read_vcf_filename);
    ////place_reads(read_map, hap_map, hap_abun_map, hap_read_map);
    
    //std::unordered_map<std::string, std::string> hap_clade_map;
    //read_csv(hap_clade_map, hap_clade_csv_filename);
    //compute_abundance(hap_abun_map, hap_clade_map);
    compute_distance(T, hap_map, vcf_samples);


///////////////////////NEW TEST
    //std::unordered_map<int, struct read_info*> read_map;
    //read_vcf(read_map, read_vcf_filename);
    //std::vector<int> read_idx;
    //auto rd_itr = read_map.begin();
    //while (rd_itr != read_map.end()) {
    //    if (rd_itr->second->mutations.size() == 0)
    //        read_idx.emplace_back(rd_itr->first);
    //    rd_itr++;
    //}
    //printf("\n Repeated Reads: %d\n", (int)read_idx.size());

}


void place_reads(const std::unordered_map<int, struct read_info*>& read_map, const std::unordered_map<int, struct read_info*>& hap_map, const std::unordered_map<std::string, double>& hap_abun_map, tbb::concurrent_hash_map<std::string, std::vector<size_t>>& hap_read_map) {
    //Calculating node score for remaining reads
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_map.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                //Find the most parsimonious positions
                auto read_id = read_map.find(i)->second;
                std::vector<std::string> epp_list;
                int min_dist = std::numeric_limits<int>::max();
                auto hap_itr = hap_map.begin();
                while (hap_itr != hap_map.end()) {
                    std::vector<MAT::Mutation> mutations_in_range;
                    for (auto m: hap_itr->second->mutations)
                        if ((m.position >= read_id->start) && (m.position <= read_id->end))
                            mutations_in_range.emplace_back(m);
                    int curr_dist = mutation_distance(read_id->mutations, mutations_in_range);
                    mutations_in_range.clear();
                    if (curr_dist <= min_dist) {
                        if (curr_dist < min_dist) {
                            epp_list.clear();
                            min_dist = curr_dist;
                        }
                        epp_list.emplace_back(hap_itr->second->read);
                    }
                    hap_itr++;
                }
                //Get probabilities of all the EPP entries
                std::vector<double> prob_list;
                double total_prob = 0.0;
                for (auto epp: epp_list) {
                    double prob = hap_abun_map.find(epp)->second;
                    prob_list.emplace_back(prob);
                    total_prob += prob;
                }

                // Create a random number generator
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> dis(0.0, 1.0);
                // Generate a random number between 0 and 1
                double randomValue = dis(gen);
                // Assign the 'read' variable based on probabilities
                double cumulativeProbability = 0.0;
                int selectedIdx = -1;
                for (size_t j = 0; j < prob_list.size(); ++j) {
                    cumulativeProbability += (prob_list[j] / total_prob);
                    if (randomValue < cumulativeProbability) {
                        selectedIdx = j;
                        break;
                    }
                }

                //Assign Read -> Haplotype
                std::string haplotype = epp_list[selectedIdx]; 
                std::vector<size_t> new_vec = {i};
                tbb::concurrent_hash_map<std::string, std::vector<size_t>>::accessor ac;
                auto created = hap_read_map.insert(ac, std::make_pair(haplotype, new_vec));
                if (!created)
                    ac->second.emplace_back(i);
                ac.release();
            }
        },
    ap);
    
    //Check mutation region of each haplotype
    auto hrm_itr = hap_read_map.begin();
    while (hrm_itr != hap_read_map.end()) {
        std::unordered_map<int, std::vector<std::unordered_map<int8_t, int>>> read_mut_map;
        for (auto r_idx: hrm_itr->second) {
            for (auto mut: read_map.find(r_idx)->second->mutations) {
                //Mutation seen at this position
                if (read_mut_map.find(mut.position) != read_mut_map.end()) {
                    auto& rm_vec = read_mut_map[mut.position];
                    bool found = false;
                    for (auto& rm_map: rm_vec) {
                        // Increment the mutation count if found
                        if (rm_map.find(mut.mut_nuc) != rm_map.end()) {
                            rm_map[mut.mut_nuc]++;
                            found = true;
                            break;
                        }
                    }
                    // Add new mutation to the list if not found
                    if (!found) {
                        std::unordered_map<int8_t, int> mut_count = {{mut.mut_nuc, 1}};
                        rm_vec.emplace_back(mut_count);
                    }
                    
                }
                //Mutation not present in map
                else {
                    std::unordered_map<int8_t, int> mut_count = {{mut.mut_nuc, 1}};
                    std::vector<std::unordered_map<int8_t, int>> mut_count_list = {mut_count};
                    read_mut_map[mut.position] = mut_count_list;
                }
            }
        }

        //Get Haplotype mutations
        std::vector<MAT::Mutation> hap_mut;
        for (const auto& h_m: hap_map) {
            if (h_m.second->read == hrm_itr->first) {
                hap_mut = h_m.second->mutations;
                break;
            }
        }
        printf("\n\nHAP: %s, Total_reads: %d, hap_muts: %d\n", hrm_itr->first.c_str(), (int)hrm_itr->second.size(), (int)hap_mut.size());
        for (auto mut: hap_mut)
            printf(" %d%c", mut.position, MAT::get_nuc(mut.mut_nuc));

        //Print mutations
        for (const auto& r_mut: read_mut_map) {
            //Count all reads that have this position
            int total_reads = 0;
            for (const auto r_idx: hrm_itr->second){
                if ((r_mut.first >= read_map.find(r_idx)->second->start) && (r_mut.first <= read_map.find(r_idx)->second->end))
                    total_reads++;
            }
            //Get the mutations present tat this position
            for (const auto& allele_freq_list: r_mut.second) {
                for (const auto& allele_freq: allele_freq_list) {
                    bool present = false;
                    for (auto m: hap_mut)
                        if ((m.position == r_mut.first) && (m.mut_nuc == allele_freq.first)) {
                            present = true;
                            break;
                        } 
                    if (present)
                        printf("\n%d%c -> %f (%d) PRESENT", r_mut.first, MAT::get_nuc(allele_freq.first), ((double)allele_freq.second / (double)total_reads), allele_freq.second);
                    else
                        printf("\n%d%c -> %f (%d) ABSENT", r_mut.first, MAT::get_nuc(allele_freq.first), ((double)allele_freq.second / (double)total_reads), allele_freq.second);
                }
            }
        }
        hrm_itr++;
        read_mut_map.clear();
    }
}


void compute_distance(const MAT::Tree &T, const std::unordered_map<size_t, struct read_info*> &hap_map, const std::vector<std::string> &vcf_samples) {
    fprintf(stderr, "Haplotypes: %d\n\n", (int)hap_map.size());
    printf("\nMUTATION DISTANCE NEW:\n");
    //Get Mutations of samples
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = get_mutations(T, sample);
        std::string best_node = "";
        auto itr = hap_map.begin();
        while (itr != hap_map.end()) {
            int curr_dist = mutation_distance(sample_mutations, itr->second->mutations);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = itr->second->read;
            }
            itr++;
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node.c_str(), min_dist);
        
    }

    //Depth first expansion to get all nodes in the tree and 
    std::vector<MAT::Node*> dfs = T.depth_first_expansion(); 
    
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
void compute_abundance(const std::unordered_map<std::string, double>& hap_abun_map, const std::unordered_map<std::string, std::string>& hap_clade_map) {
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