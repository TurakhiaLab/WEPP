#include "wbe.hpp"

void detectPeaks (po::parsed_options parsed) {
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
    std::string vcf_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.vcf";
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string hap_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotype_abundance.csv";
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
    readVCF(read_map, vcf_filename_reads, ref_seq.size(), false);
    
    //CREATE new tree containg only selected_lineage_list
    //Get haplotype abundances
    std::unordered_map<std::string, double> hap_abun_map;
    readCSV(hap_abun_map, hap_csv_filename);
    //Extract lineages from haplotype names
    std::vector<std::string> selected_lineage_list; 
    for (const auto& hap_abun: hap_abun_map) {
        size_t last_underscore = hap_abun.first.find_last_of('_');
        std::string curr_lineage = hap_abun.first.substr(last_underscore + 1);
        if (std::find(selected_lineage_list.begin(), selected_lineage_list.end(), curr_lineage) == selected_lineage_list.end())
            selected_lineage_list.emplace_back(curr_lineage);
    }
    
    //CREATE smaller lineage tree
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    MAT::Tree T_new;
    createLineageTree(T.root, selected_lineage_list, T_new);
    analyzeReads(T_new, T, read_map, node_score_map, vcf_samples, barcode_file, read_mutation_depth_vcf, ref_seq.size());
}

//Main peak search algorithm
void analyzeReads(const MAT::Tree &T, const MAT::Tree &T_ref, const std::unordered_map<size_t, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, const std::vector<std::string> &vcf_samples, const std::string &barcode_file, const std::string &read_mutation_depth_vcf, const size_t &seq_len) {
    timer.Start();
    int top_n = 25, m_dist_thresh = 2, neighbor_dist_thresh = 7, neighbor_peaks_thresh = 100, tree_range = 600, tree_increment = 400;
    std::vector<MAT::Node*> peak_nodes, curr_peak_nodes, prohibited_nodes, neighbor_nodes, curr_neighbor_nodes;
    std::vector<size_t> remaining_reads(read_map.size()), remove_reads;
    
    //GREEDY ALGORITHM for peak nodes
    //    1. Place remaining reads on tree
    //    2. Remove nodes not to be considered (present in Prohibited nodes)
    //    3. Sort and Find top scoring nodes
    //    4. Don't consider top scoring nodes in the neighborhood of other top nodes
    //    5. Remove reads mapped to these nodes
    //    6. Add neighbors of selected top scoring nodes to Prohibited nodes 
    
    //Initializing remaining_reads
    for (size_t i = 0; i < read_map.size(); i++)
        remaining_reads[i] = i;
    
    //PREVENTS "DUMMY-OG" from becoming a Peak 
    prohibited_nodes.emplace_back(T.root);

    //ITERATE till no reads left
    while ((int)remaining_reads.size() > 0) {
        printf("\n");
        fprintf(stderr, "\n");
        
        //MAP reads to nodes
        placeReadHelper(T.root, read_map, remaining_reads, curr_peak_nodes, node_score_map, remove_reads, seq_len, tree_increment, tree_range);
        
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
        sortNodeScore(T, node_score_map, node_score_vector);

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
            auto curr_clade = getLineage(T, ref_n_s.first);
            printf("PEAK: %s, Score: %f, Clade:%s, reads: %d\n",ref_n_s.first->identifier.c_str(), ref_n_s.second, curr_clade.c_str(), (int)remaining_reads.size());
            fprintf(stderr,"PEAK: %s, Score: %f, Clade:%s, reads: %d\n",ref_n_s.first->identifier.c_str(), ref_n_s.second, curr_clade.c_str(), (int)remaining_reads.size());
            if ((++nodes_found) == top_n)
                break;
            //Do neighbor check
            int num_check_nodes = (int)top_n_node_scores.size() - idx - 1;
            tbb::parallel_for(tbb::blocked_range<int>(0, num_check_nodes),
                [&](tbb::blocked_range<int> k) {
                    for (int i = k.begin(); i < k.end(); ++i) {
                        auto curr_n_s = top_n_node_scores[i + idx + 1];
                        //Only consider curr_n_s if mutation_distance > m_dist_thresh
                        int m_dist = mutationDistance(T, T, ref_n_s.first, curr_n_s.first);
                        if (m_dist <= m_dist_thresh) {
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
        curr_neighbor_nodes = updateNeighborNodes(T, curr_peak_nodes, peak_nodes, node_score_map, neighbor_nodes, neighbor_dist_thresh, neighbor_peaks_thresh);
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
        updateProhibitedNodes(T, curr_peak_nodes, prohibited_nodes, m_dist_thresh);
        
        //REMOVE reads mapped to curr_peak_nodes
        placeReadHelper(T.root, read_map, remaining_reads, curr_peak_nodes, node_score_map, remove_reads, seq_len, tree_increment, tree_range);

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
    
    printf("MUTATION DISTANCE ORIG:\n"); 
    //Verify Recovery of Input Samples
    timer.Start();
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_node = T_ref.get_node(sample);
        MAT::Node* best_node = NULL;
        for (auto pn: peak_nodes) {
            int curr_dist = mutationDistance(T, T_ref, pn, sample_node);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = pn;
            }
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld msec\n\n", timer.Stop());
    
    //ADD neighbor_nodes to peak_nodes
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();
    
    generateFilteringData(T, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf);
}

//Get MAT having ONLY desired lineages
void createLineageTree(MAT::Node* ref_root, const std::vector<std::string> &lineage_list, MAT::Tree &T) {
    //Have a queue of current_node (from ref_Tree)and parent_node (from new_Tree) pair
    std::queue<struct node_pair_clade> remaining_nodes;
    auto new_node = T.create_node("DUMMY-OG", -1.0, ref_root->clade_annotations.size());
    struct node_pair_clade cur_node_pair_clade = {ref_root, new_node, ""};
    remaining_nodes.push(cur_node_pair_clade);
        
    //Include the curr_node to the tree if it belongs to the lineage_list
    while(remaining_nodes.size() > 0) {
        cur_node_pair_clade = remaining_nodes.front();
        auto r_curr_node = cur_node_pair_clade.ref_tree_node;
        auto n_parent_node = cur_node_pair_clade.new_tree_parent;
        auto r_par_lineage = cur_node_pair_clade.ref_tree_parent_lineage;       
        remaining_nodes.pop();
        
        std::string curr_clade = r_curr_node->clade_annotations[1];
        if(curr_clade == "")
            curr_clade = r_par_lineage;
        //Add a new node to tree if curr_clade belongs to lineage_list 
        if (std::find(lineage_list.begin(), lineage_list.end(), curr_clade) != lineage_list.end()) {
            //Create a new_node
            auto new_node = T.create_node(r_curr_node->identifier, n_parent_node, -1.0);
            //Add lineage name
            new_node->clade_annotations[1] = r_curr_node->clade_annotations[1];
            //Add mutations to new_node
            new_node->mutations.reserve(new_node->mutations.size() + r_curr_node->mutations.size());
            new_node->mutations.insert(new_node->mutations.end(), r_curr_node->mutations.begin(), r_curr_node->mutations.end());

            //Add mutations between parent node and current node if this is lineage-root
            if (curr_clade != r_par_lineage) {
                auto n_parent_name = n_parent_node->identifier;
                auto anc = r_curr_node->parent;
                while (anc != NULL) {
                    //Check if anc_clade is same as n_parent_name
                    if (anc->identifier == n_parent_name)
                        break;
                    //Iterate through all mutations of anc
                    for (const auto& anc_mut: anc->mutations) {
                        auto curr_mut_itr = new_node->mutations.begin();
                        while (curr_mut_itr != new_node->mutations.end()) {
                            if (curr_mut_itr->position == anc_mut.position) {
                                curr_mut_itr->par_nuc = anc_mut.par_nuc;
                                break;
                            }
                            curr_mut_itr++;
                        }
                        if (curr_mut_itr == new_node->mutations.end())
                            new_node->mutations.emplace_back(anc_mut);
                    }
                    anc = anc->parent;
                }
                //Remove mutations having same mut_nuc and par_nuc
                auto curr_mut_itr = new_node->mutations.begin();
                while (curr_mut_itr != new_node->mutations.end()) {
                    if (curr_mut_itr->mut_nuc == curr_mut_itr->par_nuc)
                        curr_mut_itr = new_node->mutations.erase(curr_mut_itr);
                    else
                        curr_mut_itr++;
                }
            }
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children) {
                struct node_pair_clade child_node_pair_clade = {child, new_node, curr_clade};
                remaining_nodes.push(child_node_pair_clade);
            }
        }
        //Continue search if not beloning to lineage_list
        else {
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children) {
                struct node_pair_clade child_node_pair_clade = {child, n_parent_node, curr_clade};
                remaining_nodes.push(child_node_pair_clade);
            }
        }
    }
}

//Add neighboring nodes to prohibited nodes
void updateProhibitedNodes(const MAT::Tree &T, const std::vector<MAT::Node*> &curr_peak_nodes, std::vector<MAT::Node*> &prohibited_nodes, const int& m_dist_thresh) {
    ////ADDING current_peak_nodes and its neighbors to prohibited_nodes
    ////  1. Find the farthest ancestor of every peak node within m_dist_thresh
    ////  2. Recursively only analyze its children within m_dist_thresh
    ////  3. Only include unseen neighborhood nodes to prohibited_nodes
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, curr_peak_nodes.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                auto peak = curr_peak_nodes[i];
                //Find farthest ancestor with mut distance from peak <= m_dist_thresh
                MAT::Node* anc_node = NULL;
                for (const auto& n: T.rsearch(peak->identifier, true)) {
                    int m_dist = mutationDistance(T, T, n, peak);
                    if (m_dist <= m_dist_thresh)
                        anc_node = n;
                    else
                        break;
                }
                //Adding neighborhood peaks to prohibited_nodes
                std::queue<MAT::Node*> remaining_nodes;
                remaining_nodes.push(anc_node);
                while (remaining_nodes.size() > 0) {
                    MAT::Node* present_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    //Only add if mutation distance between peak and present_node <= m_dist_thresh
                    int m_dist = mutationDistance(T, T, present_node, peak);
                    if (m_dist <= m_dist_thresh) {
                        //Only add if present_node NOT already present in prohibited_nodes
                        //Looking for exact node matches here as we want to remove same mutation nodes
                        {   
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            if (std::find(prohibited_nodes.begin(), prohibited_nodes.end(), present_node) == prohibited_nodes.end())
                                prohibited_nodes.emplace_back(present_node);
                        }
                        //Add present_node's children in list for checking
                        for (const auto& c: present_node->children)
                            remaining_nodes.push(c);
                    }
                }
            }
        },
    ap);
}

//Add neighboring nodes to current peak
std::vector<MAT::Node*> updateNeighborNodes(const MAT::Tree &T, const std::vector<MAT::Node*> &curr_peak_nodes, const std::vector<MAT::Node*> &peak_nodes, const tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, std::vector<MAT::Node*> &neighbor_nodes, const int& neighbor_dist_thresh, const int& neighbor_peaks_thresh) {
    std::vector<MAT::Node*> final_neighbors; 
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex_outer;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, curr_peak_nodes.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                std::vector<std::pair<MAT::Node*, double>> potential_neighbor_node_scores;
                std::vector<MAT::Node*> potential_neighbor_nodes; 
                std::queue<MAT::Node*> remaining_nodes;
                auto peak = curr_peak_nodes[i];
                //Find farthest ancestor with mut distance from peak <= neighbor_dist_thresh
                MAT::Node* anc_node = NULL;
                for (const auto& n: T.rsearch(peak->identifier, true)) {
                    int m_dist = mutationDistance(T, T, n, peak);
                    if (m_dist <= neighbor_dist_thresh)
                        anc_node = n;
                    else
                        break;
                }

                //Adding neighborhood peaks to potential_neighbor_nodes
                remaining_nodes.push(anc_node);
                while (remaining_nodes.size() > 0) {
                    MAT::Node* present_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    //Only add if present_node is unique AND mutation distance between peak and present_node <= neighbor_dist_thresh
                    int m_dist = mutationDistance(T, T, present_node, peak);
                    if (m_dist <= neighbor_dist_thresh) {
                        if (m_dist)
                            potential_neighbor_nodes.emplace_back(present_node);
                        //ADD present_node's children in list for checking
                        for (const auto& c: present_node->children)
                            remaining_nodes.push(c);
                    }
                }
                
                //Selecting unseen nodes from potential_neighbor_nodes to potential_neighbors_node_score
                using my_mutex_t = tbb::queuing_mutex;
                my_mutex_t my_mutex;
                static tbb::affinity_partitioner ap_1;
                tbb::parallel_for(tbb::blocked_range<size_t>(0, potential_neighbor_nodes.size()),
                    [&](tbb::blocked_range<size_t> l) {
                        for (size_t j = l.begin(); j < l.end(); ++j) {
                            using rwmutex_t = tbb::queuing_rw_mutex;
                            rwmutex_t my_mutex_rw;
                            static tbb::affinity_partitioner ap_2;
                            bool found = false;
                            auto present_node = potential_neighbor_nodes[j];

                            //Check if similar present_node NOT already included in curr_peak_nodes
                            for (const auto &hap: curr_peak_nodes) {
                                int mut_dist = mutationDistance(T, T, hap, present_node);
                                if (!mut_dist) {
                                    found = true;
                                    break;
                                }
                            }

                            //Check if similar present_node NOT already included in peak_nodes
                            if (!found) {
                                tbb::parallel_for(tbb::blocked_range<size_t>(0, peak_nodes.size()),
                                    [&](tbb::blocked_range<size_t> m) {
                                        for (size_t h = m.begin(); h < m.end(); ++h) {
                                            {
                                                rwmutex_t::scoped_lock my_lock{my_mutex_rw, false};
                                                if (found)
                                                    continue;

                                            }
                                            auto hap = peak_nodes[h];
                                            int mut_dist = mutationDistance(T, T, hap, present_node);
                                            if (!mut_dist) {
                                                rwmutex_t::scoped_lock my_lock{my_mutex_rw, true};
                                                found = true;
                                            }
                                        }
                                    },
                                ap_2);
                            }

                            //Check if similar present_node NOT already included in neighbor_nodes
                            if (!found) {
                                tbb::parallel_for(tbb::blocked_range<size_t>(0, neighbor_nodes.size()),
                                    [&](tbb::blocked_range<size_t> m) {
                                        for (size_t h = m.begin(); h < m.end(); ++h) {
                                            {
                                                rwmutex_t::scoped_lock my_lock{my_mutex_rw, false};
                                                if (found)
                                                    continue;

                                            }
                                            auto hap = neighbor_nodes[h];
                                            int mut_dist = mutationDistance(T, T, hap, present_node);
                                            if (!mut_dist) {
                                                rwmutex_t::scoped_lock my_lock{my_mutex_rw, true};
                                                found = true;
                                            }
                                        }
                                    },
                                ap_2);
                            }
                                
                            //Only ADD if similar node NOT seen before
                            if (!found) {
                                tbb::concurrent_hash_map<MAT::Node*, double>::const_accessor k_ac;
                                if (node_score_map.find(k_ac, present_node)) {
                                    my_mutex_t::scoped_lock my_lock{my_mutex};
                                    potential_neighbor_node_scores.emplace_back(std::make_pair(present_node, k_ac->second));
                                }
                            }

                        }
                    },
                ap_1);
                potential_neighbor_nodes.clear();
                
                //SORT potential_neighbor_node_scores
                tbb::parallel_sort(potential_neighbor_node_scores.begin(), potential_neighbor_node_scores.end(), 
                     [&T](const auto& a, const auto& b) {
                        return compareNodeScore(T, a, b);
                });
                
                //ADD Unique top nodes to final_neighbors
                int neighbors_added = 0, idx = 0;
                while (neighbors_added < neighbor_peaks_thresh) {
                    if (idx == (int)potential_neighbor_node_scores.size())
                        break;
                    auto neighbor_node = potential_neighbor_node_scores[idx++].first;
                    {
                        my_mutex_t::scoped_lock my_lock{my_mutex_outer};
                        if (std::find(final_neighbors.begin(), final_neighbors.end(), neighbor_node) == final_neighbors.end()) {
                            final_neighbors.emplace_back(neighbor_node);
                            neighbors_added++;
                        }
                    }
                }
            }
        },
    ap);

    return final_neighbors;
}

//Function to read csv containing abundance of haplotypes
void readCSV(std::unordered_map<std::string, double>& hap_abun_map, const std::string &hap_csv_filename) {
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

//Get the names of samples prsent in samples.vcf
void readSampleVCF(std::vector<std::string> &vcf_samples, const std::string vcf_filename_samples) {
    // Boost library used to stream the contents of the input VCF file
    boost::filesystem::ifstream fileHandler(vcf_filename_samples);
    std::string s;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            //Checking for header
            if (words[1] == "POS") {
                //Leave certain fields based on our VCF format
                for (int j=9; j < (int)words.size(); j++)
                    vcf_samples.emplace_back(words[j]);
            }
        }
    }
}