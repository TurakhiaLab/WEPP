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
    std::string hap_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotype_abundance.csv";
    std::string hap_vcf_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotypes.vcf";
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
    
    //Get the curr_peak_nodes
    timer.Start();
    int tree_range = 600, tree_increment = 400, node_lim = 5, prohibited_dist_thresh = 3, neighbor_dist_thresh = 5;
    MAT::Tree T_condensed;
    std::vector<MAT::Node*> prev_peak_nodes, curr_peak_nodes, peaks_and_neighbors;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    tbb::concurrent_vector<std::pair<MAT::Node*, double>> node_score_vector;
    std::vector<size_t> remaining_read_ids(read_map.size()), remove_read_ids;
    std::unordered_map<std::string, double> hap_abun_map;
    std::unordered_map<std::string, std::vector<std::string>> condensed_nodeNames_map;
    
    //Initializing remaining_read_ids
    for (size_t i = 0; i < read_map.size(); i++)
        remaining_read_ids[i] = i;

    //Create smaller tree with only sites covered by reads 
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);
    fprintf(stderr, "Condensed tree created in %ld sec \n\n", (timer.Stop() / 1000));
    
    //MAP reads to nodes
    timer.Start();
    
    placeReadHelper(T_condensed.root, condensed_node_mappings, read_map, remaining_read_ids, curr_peak_nodes, node_score_map, remove_read_ids, ref_seq.size(), tree_increment, tree_range, true);
    fprintf(stderr, "Read mapping completed in %ld min \n\n", (timer.Stop() / (60 * 1000)));
    
    timer.Start();
    //SORT nodes
    sortNodeScore(condensed_node_mappings, node_score_map, node_score_vector);
    //node_score_map.clear();
    
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
                if (ln_itr == lineage_nodes_count.end()) {
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
                curr_peak_nodes.emplace_back(n_s.first);
                //Do neighbor check
                getProhibitedNodes(T_condensed, T, condensed_node_mappings, n_s.first, curr_prohibited_nodes, prohibited_dist_thresh);
                for (const auto& p_node: curr_prohibited_nodes)
                    prohibited_nodes.insert(p_node);
                curr_prohibited_nodes.clear();
            }
        }
    }

    node_score_vector.clear();
    lineage_nodes_count.clear();
    fprintf(stderr, "Lineage selection completed in %ld min \n\n", (timer.Stop() / (60 * 1000)));


    //ITERATIVE Filtering and adding neighbors
    //Sort curr_peak_nodes
    tbb::parallel_sort(curr_peak_nodes.begin(), curr_peak_nodes.end(), 
        [](const auto& a, const auto& b) {
            return  a->identifier > b->identifier;
    });
    peaks_and_neighbors = curr_peak_nodes;
    curr_peak_nodes.clear();

    timer.Start();
    int loop_count = 0;
    while(loop_count++ < 20) {
        if (loop_count > 1) {
            //REMOVE previous files
            std::string command = "rm " + std::string(barcode_file) + " " + std::string(read_mutation_depth_vcf) + " " + std::string(condensed_nodes_csv) + " " + std::string(hap_csv_filename) + " " + std::string(hap_vcf_filename);
            int result = std::system(command.c_str());
            if (result)
                fprintf(stderr, "\nCannot remove files\n");
        }
        
        generateFilteringData(T, T_condensed, condensed_node_mappings, ref_seq, peaks_and_neighbors, read_map, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
        peaks_and_neighbors.clear();

        //Run peaks_filtering
        std::string command = "python src/WBE/peaks_filtering.py " + vm["output-files-prefix"].as<std::string>() + " " + std::string(dir_prefix);
        int result = std::system(command.c_str());
        if (result)
            fprintf(stderr, "\nCannot run peak_filetring.py\n");
    
        //Read haplotype and condensed node names
        readCSV(hap_abun_map, hap_csv_filename);
        readCSV(condensed_nodeNames_map, condensed_nodes_csv);

        //Convert hap_map to curr_peak_nodes list and include CONDENSED nodes
        std::string check_string = "CONDENSED";
        for (const auto& hap_abun: hap_abun_map) {
            std::string node_name;
            if (hap_abun.first.find(check_string) != std::string::npos)
                node_name = condensed_nodeNames_map[hap_abun.first].front();
            else
                node_name = hap_abun.first;
            
            size_t last_underscore = node_name.find_last_of('_');
            std::string real_node_name = node_name.substr(0, last_underscore);
            while (T_condensed.get_node(real_node_name) == NULL) {
                last_underscore = real_node_name.find_last_of('_');
                real_node_name = real_node_name.substr(0, last_underscore);
            }
            curr_peak_nodes.emplace_back(T_condensed.get_node(real_node_name));
        }
        hap_abun_map.clear();
        condensed_nodeNames_map.clear();

        //Sort curr_peak_nodes
        tbb::parallel_sort(curr_peak_nodes.begin(), curr_peak_nodes.end(), 
            [](const auto& a, const auto& b) {
                return  a->identifier > b->identifier;
            });

        //Check CONVERGENCE by comparing curr_peak_nodes with prev_peak_nodes
        bool convergence = false;
        if (curr_peak_nodes.size() == prev_peak_nodes.size()) {
            convergence = true;
            for (int i = 0; i < (int)curr_peak_nodes.size(); i++) {
                if (curr_peak_nodes[i]->identifier != prev_peak_nodes[i]->identifier) {
                    convergence = false;
                    break;
                }
            }
        }
        prev_peak_nodes.clear();
        if (convergence)
            break;

        //Add neighboring nodes to peak_nodes
        prev_peak_nodes = curr_peak_nodes;
        peaks_and_neighbors = curr_peak_nodes;
        curr_peak_nodes.clear();
        addNeighborLeaves(T_condensed, T, condensed_node_mappings, node_score_map, peaks_and_neighbors, prohibited_dist_thresh, neighbor_dist_thresh);
        fprintf(stderr, "Iter - %d completed \n\n", loop_count);
    }
    prohibited_nodes.clear();
    node_score_map.clear();
    condensed_node_mappings.clear();
    fprintf(stderr, "Iterative peaks selection completed in %ld min \n\n", (timer.Stop() / (60*1000)));
}
