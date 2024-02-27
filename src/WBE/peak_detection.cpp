#include "wbe.hpp"

void detectPeaks (po::parsed_options parsed) {
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
    std::string prior_lineages = vm["prior-lineages"].as<std::string>();
    std::string vcf_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.vcf";
    std::string proto_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_sam.pb";
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string hap_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotype_abundance.csv";
    std::string barcode_file = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_barcode.csv";
    std::string condensed_nodes_csv = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_condensed_nodes.csv";
    std::string read_mutation_depth_vcf = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_read_data.vcf";
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    fprintf(stderr, "\nNum Cores: %d\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);
    
    timer.Start();
    //Get prior lineages
    std::vector<std::string> selected_lineages; 
    std::stringstream lin_str(prior_lineages);
    std::string str;
    while (std::getline(lin_str,str,',')) {
        selected_lineages.emplace_back(str);
    }
    //Loading reference genome
    std::ifstream fasta_f(ref_fasta);
    if (!fasta_f.is_open()) {
        std::cerr << "Error: Unable to open file " << ref_fasta << std::endl;
        exit(1); 
    }
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }
    fprintf(stderr, "\nLoading input MAT files %s \n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));

    //Read samples.vcf to check how close are peaks to samples 
    std::vector<std::string> vcf_samples;
    readSampleVCF(vcf_samples, vcf_filename_samples);
    //Get the input reads data

    std::unordered_map<size_t, struct read_info*> read_map;
    std::unordered_map<std::string, std::vector<std::string>> reverse_merge;
    load_reads_from_proto(proto_reads, read_map, reverse_merge);
    
    //Get haplotype abundances and condensed node names
    timer.Start();
    std::unordered_map<std::string, double> hap_abun_map;
    std::unordered_map<std::string, std::vector<std::string>> condensed_nodeNames_map;
    std::vector<MAT::Node*> curr_peak_nodes;
    readCSV(hap_abun_map, hap_csv_filename);
    readCSV(condensed_nodeNames_map, condensed_nodes_csv);

    //CREATE condensed tree for neighbor lineage search
    int neighbor_dist_thresh = 4;
    MAT::Tree T_condensed;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);
    
    //EXTRACT lineages from haplotype names
    std::string check_string = "CONDENSED";
    for (const auto& hap_abun: hap_abun_map) {
        if (hap_abun.first.find(check_string) != std::string::npos) {
            auto node_names_list = condensed_nodeNames_map[hap_abun.first];
            for (int i = 0; i < (int)node_names_list.size(); i++) {
                auto node_name = node_names_list[i];
                size_t last_underscore = node_name.find_last_of('_');
                std::string real_node_name = node_name.substr(0, last_underscore);
                while (T.get_node(real_node_name) == NULL) {
                    last_underscore = real_node_name.find_last_of('_');
                    real_node_name = real_node_name.substr(0, last_underscore);
                }
                if (!i)
                    curr_peak_nodes.emplace_back(T_condensed.get_node(real_node_name));
                std::string curr_lineage = node_name.substr(last_underscore + 1);
                if (std::find(selected_lineages.begin(), selected_lineages.end(), curr_lineage) == selected_lineages.end())
                    selected_lineages.emplace_back(curr_lineage);
            }

        }
        else {
            size_t last_underscore = hap_abun.first.find_last_of('_');
            std::string real_node_name = hap_abun.first.substr(0, last_underscore);
            while (T.get_node(real_node_name) == NULL) {
                last_underscore = real_node_name.find_last_of('_');
                real_node_name = real_node_name.substr(0, last_underscore);
            }
            curr_peak_nodes.emplace_back(T_condensed.get_node(real_node_name));
            std::string curr_lineage = hap_abun.first.substr(last_underscore + 1);
            if (std::find(selected_lineages.begin(), selected_lineages.end(), curr_lineage) == selected_lineages.end())
                selected_lineages.emplace_back(curr_lineage);
        }
    }

    //ADD Neighboring LINEAGES
    addNeighborLineages(T_condensed, T, condensed_node_mappings, curr_peak_nodes, selected_lineages, neighbor_dist_thresh);
    curr_peak_nodes.clear();
    condensed_node_mappings.clear();
    MAT::clear_tree(T_condensed);
    fprintf(stderr, "Added Lineages in %ld sec\n\n", (timer.Stop() / 1000));
    
    //CREATE new tree containg only selected_lineages
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    MAT::Tree T_lin;
    createLineageTree(T.root, selected_lineages, T_lin);
    selected_lineages.clear();

    analyzeReads(T, T_lin, ref_seq, read_map, node_score_map, vcf_samples, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
    
    //RUN peaks_filtering
    std::string command = "python src/WBE/peaks_filtering.py " + vm["output-files-prefix"].as<std::string>() + " " + std::string(dir_prefix) + " 100";
    int result = std::system(command.c_str());
    if (result)
        fprintf(stderr, "\nCannot run peak_filtering.py\n");
}

//Main peak search algorithm
void analyzeReads(const MAT::Tree &T_ref, const MAT::Tree &T, const std::string &ref_seq, const std::unordered_map<size_t, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, const std::vector<std::string> &vcf_samples, const std::string &barcode_file, const std::string &read_mutation_depth_vcf, const std::string &condensed_nodes_csv) {
    timer.Start();
    int top_n = 25, prohibited_dist_thresh = 1, neighbor_dist_thresh = 4, neighbor_peaks_thresh = 100, tree_increment, tree_range = 600, tree_overlap = 200, range_factor = 1;
    std::vector<MAT::Node*> peak_nodes, curr_peak_nodes, prohibited_nodes, neighbor_nodes, curr_neighbor_nodes;
    std::vector<size_t> remaining_reads(read_map.size()), remove_reads;
    MAT::Tree T_condensed;
    std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> condensed_node_mappings;
    
    //GREEDY ALGORITHM for peak nodes
    //    1. Place remaining reads on tree
    //    2. Remove nodes not to be considered (present in Prohibited nodes)
    //    3. Sort and Find top scoring nodes
    //    4. Don't consider top scoring nodes in the neighborhood of other top nodes
    //    5. Remove reads mapped to these nodes
    //    6. Add neighbors of selected top scoring nodes to Prohibited nodes 
    
    
    //Create smaller tree with only sites covered by reads 
    createCondensedTree(T.root, read_map, condensed_node_mappings, T_condensed);
    
    //Initializing remaining_reads
    for (size_t i = 0; i < read_map.size(); i++)
        remaining_reads[i] = i;

    //ITERATE till no reads left
    while ((int)remaining_reads.size() > 0) {
        printf("\n");
        //Adjusting tree_range based on remaining_reads
        if ((2 * range_factor * remaining_reads.size()) <= read_map.size()) {
            range_factor *= 2;
            tree_range = std::min(tree_range * 2, (int)ref_seq.size());
        }
        tree_increment = tree_range - tree_overlap;
        //MAP reads to nodes
        placeReadHelper(T_condensed.root, condensed_node_mappings, read_map, remaining_reads, curr_peak_nodes, node_score_map, remove_reads, ref_seq.size(), tree_increment, tree_range, false);
        
        //REMOVE prohibited_nodes from node_score_map
        static tbb::affinity_partitioner ap;
        using my_mutex_t = tbb::queuing_mutex;
        my_mutex_t my_mutex;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, prohibited_nodes.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto p_node = prohibited_nodes[i];
                    tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
                    if (node_score_map.find(ac, p_node))
                        node_score_map.erase(ac);
                    ac.release();
                }
            },
        ap);
        //STOP search it no new node found
        if (node_score_map.empty())
            break;
        
        //SORT nodes
        tbb::concurrent_vector<std::pair<MAT::Node*, double>> node_score_vector;
        sortNodeScore(condensed_node_mappings, node_score_map, node_score_vector);

        //GET 4*top_n top_score values from node_score_vector into top_n_node_scores
        std::vector<std::pair<MAT::Node*, double>> top_n_node_scores;
        top_n_node_scores.reserve(std::min((4*top_n), (int)node_score_vector.size()));
        auto top_score = node_score_vector.begin()->second;
        //Take 4*top_n nodes for now as some of them would be removed during neighbor peaks removal
        for (int i = 0; i < std::min((4*top_n), (int)node_score_vector.size()); ++i) {
            auto n_s = node_score_vector[i];
            //Only consider 4*top_n nodes that have score equal to top node
            if (abs(top_score - n_s.second) < 1e-9)  
                top_n_node_scores.emplace_back(n_s);
            else 
                break;
        }
        node_score_vector.clear();

        //REMOVE peak from top_n_node_scores that are in each other's neighborhood
        //Finding top_n_node_scores in neighborhood
        int nodes_found = 0;
        std::vector<MAT::Node*> top_n_node_scores_remove_nodes;
        for (int idx = 0; idx < (int)top_n_node_scores.size(); idx++) {
            auto ref_n_s = top_n_node_scores[idx];
            if (std::find(top_n_node_scores_remove_nodes.begin(), top_n_node_scores_remove_nodes.end(), ref_n_s.first) != top_n_node_scores_remove_nodes.end())
                continue;
            //Proceed further if ref_n_s not in nieghborhood
            curr_peak_nodes.emplace_back(ref_n_s.first);
            auto curr_clade = getLineage(T, condensed_node_mappings[ref_n_s.first].front());
            printf("PEAK: %s, Score: %f, Clade:%s, reads: %d\n",ref_n_s.first->identifier.c_str(), ref_n_s.second, curr_clade.c_str(), (int)remaining_reads.size());
            //fprintf(stderr,"PEAK: %s, Score: %f, Clade:%s, reads: %d\n",ref_n_s.first->identifier.c_str(), ref_n_s.second, curr_clade.c_str(), (int)remaining_reads.size());
            if ((++nodes_found) == top_n)
                break;
            //Do neighbor check
            int num_check_nodes = (int)top_n_node_scores.size() - idx - 1;
            tbb::parallel_for(tbb::blocked_range<int>(0, num_check_nodes),
                [&](tbb::blocked_range<int> k) {
                    for (int i = k.begin(); i < k.end(); ++i) {
                        auto curr_n_s = top_n_node_scores[i + idx + 1];
                        //Only consider curr_n_s if mutation_distance > prohibited_dist_thresh
                        int m_dist = mutationDistance(T_condensed, T_condensed, ref_n_s.first, curr_n_s.first);
                        if (m_dist <= prohibited_dist_thresh) {
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            top_n_node_scores_remove_nodes.emplace_back(curr_n_s.first);
                        }
                    }
                },
            ap);
        }
        top_n_node_scores_remove_nodes.clear();
        top_n_node_scores.clear();
        
        //ADD neighbors of curr_peak_nodes to neighbor_nodes

        curr_neighbor_nodes = updateNeighborNodes(T_condensed, condensed_node_mappings, curr_peak_nodes, peak_nodes, node_score_map, neighbor_nodes, neighbor_dist_thresh, neighbor_peaks_thresh);
        node_score_map.clear();
        neighbor_nodes.reserve(neighbor_nodes.size() + curr_neighbor_nodes.size());
        neighbor_nodes.insert(neighbor_nodes.end(), curr_neighbor_nodes.begin(), curr_neighbor_nodes.end());

        //ADD curr_neighbor_nodes to prohibited_nodes
        for (const auto &neighbor: curr_neighbor_nodes) {
            if (std::find(prohibited_nodes.begin(), prohibited_nodes.end(), neighbor) == prohibited_nodes.end())
                prohibited_nodes.emplace_back(neighbor);
        }
        curr_neighbor_nodes.clear();
        
        //ADD nodes to prohibited nodes list
        updateProhibitedNodes(T_condensed, curr_peak_nodes, prohibited_nodes, prohibited_dist_thresh);
        
        //REMOVE reads mapped to curr_peak_nodes
        placeReadHelper(T_condensed.root, condensed_node_mappings, read_map, remaining_reads, curr_peak_nodes, node_score_map, remove_reads, ref_seq.size(), tree_increment, tree_range, false);

        //ERASE remove_reads from remaining_reads
        tbb::parallel_sort(remove_reads.begin(), remove_reads.end());
        auto rmr_itr = remaining_reads.begin();
        auto rr_itr = remove_reads.begin();
        while ((rmr_itr != remaining_reads.end()) && (rr_itr != remove_reads.end())) {
            //Remove from remaining_reads if equal
            if (*rmr_itr == *rr_itr) {
                rr_itr++;
                rmr_itr = remaining_reads.erase(rmr_itr);
            }
            //Move to next remove_read if current remove_read not found in remaining_reads
            else if (*rmr_itr > *rr_itr) {
                rr_itr++;
                fprintf(stderr,"remove_read not present in remaining_reads!!!");
            }
            //Move to next remaining_reads if it is smaller than current remove_read
            else
                rmr_itr++;
        }
        remove_reads.clear();

        //ADD curr_peak_nodes to peak_nodes
        peak_nodes.reserve(peak_nodes.size() + curr_peak_nodes.size());
        peak_nodes.insert(peak_nodes.end(), curr_peak_nodes.begin(), curr_peak_nodes.end());

        //Clear vectors needed for next iteration
        curr_peak_nodes.clear();
    }

    fprintf(stderr,"\nRemaining Reads: %d\n", (int)remaining_reads.size());
    remaining_reads.clear();
    prohibited_nodes.clear();
    printf("\nInital PEAK nodes: %d\n", (int)peak_nodes.size());
    fprintf(stderr,"\nPeak search took %ld min\n\n", (timer.Stop() / (60 * 1000)));
    
    //Verify Recovery of Input Samples
    printf("MUTATION DISTANCE ORIG:\n"); 
    timer.Start();
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map.find(i)->second;
        for (int j = rp->start; j <= rp->end; j++)
            site_read_map.insert(j);
    }

    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = getMutations(T_ref, sample);
        //Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end()) {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        MAT::Node* best_node = NULL;
        for (const auto &pn: peak_nodes) {
            for (const auto &node: condensed_node_mappings.find(pn)->second) {
                //Getting node_mutations from the Tree
                auto node_mutations = getMutations(T_ref, node->identifier);
                //Remove mutations from node_mutations that are not present in site_read_map
                auto mut_itr = node_mutations.begin();
                while (mut_itr != node_mutations.end()) {
                    if (site_read_map.find(mut_itr->position) == site_read_map.end())
                        mut_itr = node_mutations.erase(mut_itr);
                    else
                        mut_itr++;
                }

                int curr_dist = mutationDistance(sample_mutations, node_mutations);
                if (curr_dist < min_dist) {
                    min_dist = curr_dist;
                    best_node = pn;
                }
            }
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld sec\n\n", (timer.Stop() / 1000));
    
    //ADD neighbor_nodes to peak_nodes
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();
    
    generateFilteringData(T, T_condensed, condensed_node_mappings, ref_seq, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf, condensed_nodes_csv);
    condensed_node_mappings.clear();
    MAT::clear_tree(T_condensed);
}