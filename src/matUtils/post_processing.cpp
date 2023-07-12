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
    std::string hap_csv_filename = dir_prefix + vm["read-vcf"].as<std::string>() + "_hap_abundance.csv";
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    tbb::task_scheduler_init init(num_threads);

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
    read_sample_vcf(vcf_samples, sample_vcf_filename);
    read_vcf(T, dfs, read_map, hap_vcf_filename);
    read_csv(hap_abun_map, hap_csv_filename);
    compute_abundance(T, read_map, hap_abun_map);
    compute_distance(T, read_map, vcf_samples);
}

void compute_distance(const MAT::Tree &T, const std::unordered_map<int, struct read_info*> &read_map, const std::vector<std::string> &vcf_samples) {
    fprintf(stderr, "Haplotypes: %d\n\n", (int)read_map.size());
    printf("\nMUTATION DISTANCE NEW:\n");
    //Get Mutations of samples
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
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

//Computes Lineage Abundance
void compute_abundance(const MAT::Tree& T, const std::unordered_map<int, struct read_info*> &read_map, const std::unordered_map<std::string, double>& hap_abun_map) {
    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    tbb::concurrent_hash_map<std::string, double> lineage_score;
    std::vector<MAT::Node*> peak_nodes_dummy, dfs = T.depth_first_expansion(T.root); ;
    std::vector<int> mismatch_vector_dummy;
    
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_map.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                auto read_id = read_map.find(i)->second;
                place_reads(T, dfs, read_id, NULL, node_score, peak_nodes_dummy, mismatch_vector_dummy, hap_abun_map, lineage_score, 0);
            }
        },
    ap);

    double total_score = 0.0;
    auto ls_itr = lineage_score.begin();
    while (ls_itr != lineage_score.end()) {
        total_score += ls_itr->second;
        ls_itr++;
    }
    printf("\nLINEAGE ABUNDANCE (New Peaks)\n");
    ls_itr = lineage_score.begin();
    while (ls_itr != lineage_score.end()) {
        printf("Lineage: %s, Abundance: %f\n", ls_itr->first.c_str(), ls_itr->second / total_score);
        ls_itr++;
    }
}