#include "wbe.hpp"

void filterLineages (po::parsed_options parsed) {
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
    std::string input_mat_filename = dir_prefix + vm["input-mat"].as<std::string>();
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string mismatch_matrix_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_mismatch_matrix.csv";
    std::string barcode_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_barcode.csv";
    std::string condensed_nodes_csv = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_condensed_nodes.csv";
    std::string read_mutation_depth_vcf = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_read_data.vcf";
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

    // Load input MAT and uncondense tree
    fprintf(stderr, "\nLoading input MAT files %s\n", input_mat_filename.c_str());
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));
    fprintf(stderr,"Leaves in tree: %ld\n\n", T.get_num_leaves());

    //Get the input reads data
    std::unordered_map<size_t, struct read_info*> read_map;
    readVCF(read_map, vcf_filename_reads, ref_seq.size(), true);
    
    //Get the peak_nodes
    timer.Start();
    int tree_range = 600, tree_increment = 400, node_lim = 10, m_dist_thresh = 2;
    MAT::Tree T_condensed;
    std::vector<MAT::Node*> peak_nodes;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    tbb::concurrent_vector<std::pair<MAT::Node*, double>> node_score_vector;
    std::vector<size_t> remaining_read_ids(read_map.size()), remove_read_ids;
    
    //Initializing remaining_read_ids
    for (size_t i = 0; i < read_map.size(); i++)
        remaining_read_ids[i] = i;

    //Create smaller tree with only sites covered by reads 
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);
    fprintf(stderr, "Condensed tree created in %ld sec \n\n", (timer.Stop() / 1000));
    
    //MAP reads to nodes
    timer.Start();
    placeReadHelper(T_condensed.root, condensed_node_mappings, read_map, remaining_read_ids, peak_nodes, node_score_map, remove_read_ids, ref_seq.size(), tree_increment, tree_range);
    fprintf(stderr, "Read mapping completed in %ld min \n\n", (timer.Stop() / (60 * 1000)));
    
    timer.Start();
    //SORT nodes
    sortNodeScore(condensed_node_mappings, node_score_map, node_score_vector);
    node_score_map.clear();
    
    //Only taking non-neighborhood nodes from top_n_node_scores
    std::unordered_set<MAT::Node*> prohibited_nodes;
    std::unordered_map<std::string, int> lineage_nodes_count;
    
    for (int idx = 0; idx < (int)node_score_vector.size(); idx++) {
        auto n_s = node_score_vector[idx];
        if (prohibited_nodes.find(n_s.first) != prohibited_nodes.end())
            continue;
        
        //Proceed further if n_s not in nieghborhood
        if (!condensed_node_mappings[n_s.first].empty()) {
            bool selected = false;
            //Search for lineages in uncondensed nodes of original tree
            std::unordered_set<std::string> curr_lineages; 
            for (const auto& curr_node: condensed_node_mappings[n_s.first])
                curr_lineages.insert(getLineage(T, curr_node));
            for (const auto& curr_lin: curr_lineages) {
                auto ln_itr = lineage_nodes_count.find(curr_lin);
                if (ln_itr == lineage_nodes_count.end()){
                    lineage_nodes_count.insert({curr_lin, 1});
                    selected = true;
                }
                else if (ln_itr->second < node_lim) {
                    ln_itr->second++;
                    selected = true;
                }
            }
            curr_lineages.clear();
            //Add current node to peak if lineages were selected and add its neighbors to prohibited list
            if (selected) {
                std::vector<MAT::Node*> curr_prohibited_nodes;
                peak_nodes.emplace_back(n_s.first);
                //Do neighbor check
                updateProhibitedNodes(T_condensed, {n_s.first}, curr_prohibited_nodes, m_dist_thresh);
                for (const auto& curr_node: curr_prohibited_nodes)
                    prohibited_nodes.insert(curr_node);
                curr_prohibited_nodes.clear();
            }
        }
    }
    printf("\nCondensed Peak Nodes: %lu, Lineages: %lu\n\n", peak_nodes.size(), lineage_nodes_count.size());

    node_score_vector.clear();
    prohibited_nodes.clear();
    lineage_nodes_count.clear();
    fprintf(stderr, "Lineage selection completed in %ld sec \n\n", (timer.Stop() / 1000));

    generateFilteringData(T, T_condensed, condensed_node_mappings, ref_seq, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
}