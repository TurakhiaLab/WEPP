#include "wbe.hpp"

po::variables_map parseWBEcommand(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("place_read options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("output-files-prefix,v", po::value<std::string>()->default_value(""),
     "Output files prefix.")
    ("ref-fasta,f", po::value<std::string>()->default_value(""),
     "Input Fasta file representing reference sequence")
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

void filterLineages (po::parsed_options parsed) {
    //main argument for the complex extract command
    po::variables_map vm = parseWBEcommand(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
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
    std::string mismatch_matrix_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_mismatch_matrix.csv";
    std::string barcode_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_barcode.csv";
    std::string read_mutation_depth_vcf = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_peaks.vcf";
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    fprintf(stderr, "\nNum Cores: %d\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);
    
    timer.Start();
    std::ifstream fasta_f(ref_fasta);
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }
    fprintf(stderr, "\nLoading input MAT files %s\n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));
    fprintf(stderr,"Leaves in tree: %ld\n\n", T.get_num_leaves());

    //Get the input reads data
    std::unordered_map<size_t, struct read_info*> read_map;
    readVCF(read_map, vcf_filename_reads, ref_seq.size(), false);
    
    //Get the peak_nodes
    timer.Start();
    int tree_range = 600, tree_increment = 400, node_lim = 10, m_dist_thresh = 5;
    std::vector<MAT::Node*> peak_nodes;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    tbb::concurrent_vector<std::pair<MAT::Node*, double>> node_score_vector;
    std::vector<size_t> remaining_read_ids(read_map.size()), remove_read_ids;
    
    //Initializing remaining_read_ids
    for (size_t i = 0; i < read_map.size(); i++)
        remaining_read_ids[i] = i;
    //MAP reads to nodes
    placeReadHelper(T.root, read_map, remaining_read_ids, peak_nodes, node_score_map, remove_read_ids, ref_seq.size(), tree_increment, tree_range);
    fprintf(stderr, "Read mapping completed in %ld min \n\n", (timer.Stop() / (60 * 1000)));
    
    timer.Start();
    //SORT nodes
    sortNodeScore(T, node_score_map, node_score_vector);
    
    //Only taking non-neighborhood nodes from top_n_node_scores
    std::vector<MAT::Node*> prohibited_nodes, curr_prohibited_nodes;
    std::map<std::string, int> lineage_node_map;
    for (int idx = 0; idx < (int)node_score_vector.size(); idx++) {
        auto ref_n_s = node_score_vector[idx];
        if (std::find(prohibited_nodes.begin(), prohibited_nodes.end(), ref_n_s.first) != prohibited_nodes.end())
            continue;
        //Proceed further if ref_n_s not in nieghborhood
        auto curr_clade = getLineage(T, ref_n_s.first);
        auto cn_itr = lineage_node_map.find(curr_clade);
        if (cn_itr == lineage_node_map.end()){
            lineage_node_map.insert({curr_clade, 1});
            peak_nodes.emplace_back(ref_n_s.first);
        }
        else if (cn_itr->second < node_lim) {
            cn_itr->second++;
            peak_nodes.emplace_back(ref_n_s.first);
        }
        
        //Do neighbor check
        updateProhibitedNodes(T, {ref_n_s.first}, curr_prohibited_nodes, m_dist_thresh);
        tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
        for (const auto& curr_node: curr_prohibited_nodes) {
            if (node_score_map.find(ac, curr_node)) {
                node_score_map.erase(ac);
                prohibited_nodes.emplace_back(curr_node);
            }       
        }
        curr_prohibited_nodes.clear();
    }
    printf("\nPeak Nodes: %lu, Lineages: %lu\n\n", peak_nodes.size(), lineage_node_map.size());
    node_score_vector.clear();
    node_score_map.clear();
    prohibited_nodes.clear();
    lineage_node_map.clear();
    fprintf(stderr, "Lineage selection completed in %ld min \n\n", (timer.Stop() / (60 * 1000)));
    
    generateFilteringData(T, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf);
}