#include "peak_search.hpp"

po::variables_map parse_peak_search_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("place_read options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("ref-mat,j", po::value<std::string>()->required(),
     "Ref mutation-annotated tree file [REQUIRED]")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("write-file,v", po::value<std::string>()->default_value(""),
     "Output file start names. Default is full tree")
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

void peak_search(po::parsed_options parsed) {
    //main argument for the complex extract command
    //uses included code from multiple modules
    //specifically, the modules select, describe, convert, and filter support this command
    po::variables_map vm = parse_peak_search_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string ref_mat_filename = vm["ref-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    std::string vcf_filename_samples = dir_prefix + vm["write-file"].as<std::string>() + "_samples.vcf";
    std::string vcf_filename_reads = dir_prefix + vm["write-file"].as<std::string>() + "_reads.vcf";
    std::string barcode_file = dir_prefix + vm["write-file"].as<std::string>() + "_barcode.csv";
    std::string read_mutation_depth_vcf = dir_prefix + vm["write-file"].as<std::string>() + "_read_mutation_depth.vcf";
    uint32_t num_threads = vm["threads"].as<uint32_t>();

    fprintf(stderr, "\nNum Cores: %d\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);
    
    timer.Start();
    fprintf(stderr, "\nLoading input MAT files %s and %s.\n", input_mat_filename.c_str(), ref_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T, T_ref;
    if ((input_mat_filename.find(".pb\0") != std::string::npos) || (ref_mat_filename.find(".pb\0") != std::string::npos)) {
        T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T_ref = MAT::load_mutation_annotated_tree(ref_mat_filename);
        T.uncondense_leaves();
        T_ref.uncondense_leaves();
    } else if ((input_mat_filename.find(".json\0") != std::string::npos) || (ref_mat_filename.find(".json\0") != std::string::npos)) {
        T = load_mat_from_json(input_mat_filename);
        T_ref = load_mat_from_json(ref_mat_filename);
    } else {
        fprintf(stderr, "ERROR: Input file ending not recognized. Must be .json or .pb\n");
        exit(1);
    }
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));

    fprintf(stderr,"Leaves in tree: %ld\n\n", T.get_num_leaves());

    std::unordered_map<size_t, struct read_info*> read_map;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    std::vector<std::string> vcf_samples;
    //read samples.vcf to check how close are peaks to samples 
    read_sample_vcf(vcf_samples, vcf_filename_samples);
    //Get the reads.vcf data
    read_vcf(read_map, vcf_filename_reads);
    //Core algorithm 
    analyze_reads(T, T_ref, read_map, node_score, vcf_samples, barcode_file, read_mutation_depth_vcf);
}

//Get the names of samples prsent in samples.vcf
void read_sample_vcf(std::vector<std::string> &vcf_samples, const std::string vcf_filename_samples) {
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

//Stoe the Reads from the VCF in read_map
void read_vcf(std::unordered_map<size_t, struct read_info*> &read_map, const std::string vcf_filename_reads) {
    // Boost library used to stream the contents of the input VCF file
    timer.Start();
    std::vector<size_t> missing_idx;
    std::vector<struct read_info*> read_ids;
    std::string s;
    boost::filesystem::ifstream fileHandler(vcf_filename_reads);
    bool header_found = false;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            //Checking for header
            if (words[1] == "POS") {
                header_found = true;
                //Leave certain fields based on our VCF format
                for (size_t j=9; j < words.size(); j++) {
                    struct read_info * rp = new struct read_info;
                    rp->read = words[j];
                    //Get start and end positions of the Read
                    std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                    std::smatch match;
                    if (std::regex_search(words[j], match, rgx)) {
                        rp->start = std::stoi(match[1]);
                        rp->end = std::stoi(match[2]);
                    }
                    //Consider the entire genome
                    else {
                        rp->start = 1;
                        rp->end = 29903;
                    }
                    read_ids.emplace_back(rp);
                    missing_idx.emplace_back(j);
                }
            }
            else if (header_found) {
                std::vector<std::string> alleles;
                alleles.clear();
                //Checking for different alleles at a site
                MAT::string_split(words[4], ',', alleles);
                size_t k = 0;
                while (k < missing_idx.size()) {
                    size_t j = missing_idx[k];
                    auto iter = read_ids.begin();
                    std::advance(iter, k);
                    if (iter != read_ids.end()) {
                        read_map.insert({k, (*iter)});
                        MAT::Mutation m;
                        m.chrom = words[0];
                        m.position = std::stoi(words[1]);
                        //Checking the mutating allele value within the allele sizes
                        if (std::stoi(words[j]) > int(alleles.size())) {
                            fprintf(stderr, "\n\nPosition: %d, k = %d,\n", m.position, k);
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
                                } else {
                                    auto nuc = MAT::get_nuc_id(allele[0]);
                                    if (nuc == MAT::get_nuc_id('N')) {
                                        m.is_missing = true;
                                    } else {
                                        m.is_missing = false;
                                    }
                                    m.mut_nuc = nuc;
                                }
                                (*iter)->mutations.emplace_back(m);
                            }
                        } else {
                            m.is_missing = true;
                            m.mut_nuc = MAT::get_nuc_id('N');
                            (*iter)->mutations.emplace_back(m);
                        }
                    } 
                    k++;
                }
            }
        }
    }
    fprintf(stderr,"%s parsed in %ld sec\n\n", vcf_filename_reads.c_str(), (timer.Stop() / 1000));   
}

//Parsimonious placement search for reads
int place_reads(const std::vector<MAT::Node*> &dfs, struct read_info* rp, const MAT::Node* check_node, tbb::concurrent_hash_map<MAT::Node*, double> &node_score) {
    std::stack<struct parsimony> parsimony_stack;
    struct min_parsimony min_par;
    //Checking all nodes of the tree
    for (int i = 0; i < (int)dfs.size(); i++) {
        std::vector<MAT::Mutation> uniq_curr_node_mut, common_node_mut, curr_node_par_mut;
        struct parsimony curr_par;
        //Nothing in stack for first node
        if (i) {
            //Get the parsimony vector from parent
            auto parent_parsimony = parsimony_stack.top();
            while ((dfs[i]->parent != parent_parsimony.curr_node) && (!parsimony_stack.empty()))   {
                parsimony_stack.pop();
                if (!parsimony_stack.empty())
                    parent_parsimony = parsimony_stack.top();
                else 
                    fprintf(stderr, "\nERROR in Parsimony Stack!!!!\n");
            } 
            if (parsimony_stack.empty())
                fprintf(stderr, "\nERROR in Parsimony Stack!!!!\n");
            curr_node_par_mut = parent_parsimony.p_node_par;
        }
        //Checking all mutations of current node
        for (auto node_mut: dfs[i]->mutations) {
            //Only look at mutations within the read range
            if ((node_mut.position >= rp->start) && (node_mut.position <= rp->end)) {
                bool found = false;
                //Check in Mutation position found in parsimony of parent node
                for (auto par_node_mut: curr_node_par_mut) {
                    if (par_node_mut.position == node_mut.position) {
                        //mut_nuc matches => remove mutation from parsimony
                        if ((par_node_mut.mut_nuc == node_mut.mut_nuc) || (par_node_mut.mut_nuc == 0b1111)) {
                           auto itr = curr_node_par_mut.begin();
                           while (!((itr->position == par_node_mut.position) && (itr->mut_nuc == par_node_mut.mut_nuc) && (itr->ref_nuc == par_node_mut.ref_nuc))) {
                                itr++; 
                           }
                           common_node_mut.emplace_back(par_node_mut);
                           curr_node_par_mut.erase(itr);
                           //Did not work because == not defined for MAT::Mutation
                           //curr_node_par_mut.erase(std::remove(curr_node_par_mut.begin(), curr_node_par_mut.end(), par_node_mut), curr_node_par_mut.end());
                        }
                        //Update the par_nuc in parsimony if found
                        else {
                            par_node_mut.par_nuc = node_mut.mut_nuc;    
                        }
                        found = true;
                        break;
                    }
                }
                if (found)
                    continue;

                //Not found in parent parsimony
                for (auto read_mut: rp->mutations) {
                    if (read_mut.position == node_mut.position) {
                        //For root 
                        if (!i) {
                            //Mutation found in read, add to common_node_mut
                            if (read_mut.mut_nuc == node_mut.mut_nuc)
                                common_node_mut.emplace_back(read_mut);
                            //read_mut is 'N', add to common_node_mut
                            else if (read_mut.mut_nuc == 0b1111)
                                common_node_mut.emplace_back(read_mut);
                            //If mutation is different, handle the mutation as child of current node
                            else {
                                struct MAT::Mutation new_mut;
                                new_mut.position = read_mut.position;
                                new_mut.ref_nuc = read_mut.ref_nuc;
                                new_mut.par_nuc = node_mut.mut_nuc;
                                new_mut.mut_nuc = read_mut.mut_nuc;
                                curr_node_par_mut.emplace_back(new_mut);
                                //Placing it in common_mut so don't add this mut again
                                common_node_mut.emplace_back(new_mut);
                            }
                        }
                        //Otherwise, if mut_nuc is same it would have been detcted above in parent parsimony check
                        //At this step assume the mut_nuc is different and add it as a uniq mutation
                        //Reverse par_nuc and mut_nuc as it will again get flipped in uniq_curr_node_mut
                        else {
                            struct MAT::Mutation new_mut;
                            new_mut.position = read_mut.position;
                            new_mut.ref_nuc = read_mut.ref_nuc;
                            new_mut.par_nuc = read_mut.mut_nuc;
                            new_mut.mut_nuc = node_mut.mut_nuc;
                            uniq_curr_node_mut.emplace_back(new_mut);
                        }
                        found = true;
                        break;
                    }
                }
                if (found)
                    continue;

                //Niether in parent parsimony nor in read_mut
                uniq_curr_node_mut.emplace_back(node_mut);
            }
        } 

        //Adding only unseen read_mut to node parsimony for root node
        if (!i) {
            for (auto read_mut: rp->mutations) {
                bool present = false;
                //Check if mut present in common_node_mut
                auto c_itr = common_node_mut.begin();
                while (c_itr != common_node_mut.end()) {
                    if (c_itr->position == read_mut.position) {
                        if (c_itr->mut_nuc != read_mut.mut_nuc)
                            fprintf(stderr,"common_node mut does not match read_mut!!!\n");
                        present = true;
                        break;
                    }
                    c_itr++;
                }
                if (present)
                    continue;

                //Else add it to curr_node_mut if mut_nuc != 'N'
                if (read_mut.mut_nuc != 0b1111)
                    curr_node_par_mut.emplace_back(read_mut);
            }
        }


        bool placed_child = false;
        auto curr_node = dfs[i];
        
        //Place as a sibling if common_node_mut is not empty
        if ((common_node_mut.size()) && (!curr_node->is_root())) {
            //Checking min_parsimony
            int new_min_par = -1; 
            // If best_par_score is empty and curr_par_score >= limit -> CHANGE
            if ((min_par.par_list.empty()) && ((int)curr_node_par_mut.size() >= par_score_lim))
                new_min_par = 1;
            // If curr_par_score < best_par_score and curr_par_score >= limit -> CHANGE
            else if ((curr_node_par_mut.size() < min_par.par_list[0].size()) && ((int)curr_node_par_mut.size() >= par_score_lim))
                new_min_par = 1;
            // If cur_par_score == best_par_score -> APPEND
            else if ((curr_node_par_mut.size() == min_par.par_list[0].size()))
                new_min_par = 0;
        
            if (new_min_par == 1) {
                min_par.idx_list.clear();
                min_par.par_list.clear();
                min_par.idx_list.emplace_back(i);
                min_par.par_list.emplace_back(curr_node_par_mut);
            }
            else if (new_min_par == 0) {
                min_par.idx_list.emplace_back(i);
                min_par.par_list.emplace_back(curr_node_par_mut);
            }   
        }
        //Place as a child if current node is NOT leaf node or a ROOT node
        else if (!curr_node->is_leaf()) {
            placed_child = true;
            //Updating parsimony to be stored as a child
            for (auto uniq_mut: uniq_curr_node_mut) {
                //Just reverse mut_nuc and par_nuc and add it to parsimony
                int8_t temp = uniq_mut.par_nuc;
                uniq_mut.par_nuc = uniq_mut.mut_nuc;
                uniq_mut.mut_nuc = temp;
                curr_node_par_mut.emplace_back(uniq_mut);
            }
            //Checking min_parsimony
            int new_min_par = -1; 
            // If best_par_score is empty and curr_par_score >= limit -> CHANGE
            if ((min_par.par_list.empty()) && ((int)curr_node_par_mut.size() >= par_score_lim))
                new_min_par = 1;
            // If curr_par_score < best_par_score and curr_par_score >= limit -> CHANGE
            else if ((curr_node_par_mut.size() < min_par.par_list[0].size()) && ((int)curr_node_par_mut.size() >= par_score_lim))
                new_min_par = 1;
            // If cur_par_score == best_par_score -> APPEND
            else if ((curr_node_par_mut.size() == min_par.par_list[0].size()))
                new_min_par = 0;
        
            if (new_min_par == 1) {
                min_par.idx_list.clear();
                min_par.par_list.clear();
                min_par.idx_list.emplace_back(i);
                min_par.par_list.emplace_back(curr_node_par_mut);
            }
            else if (new_min_par == 0) {
                min_par.idx_list.emplace_back(i);
                min_par.par_list.emplace_back(curr_node_par_mut);
            }
        }

        //Only update if not placed as child
        if (!placed_child) {
            //Updating parsimony to be stored as a child
            for (auto uniq_mut: uniq_curr_node_mut) {
                //Just reverse mut_nuc and par_nuc and add it to parsimony
                int8_t temp = uniq_mut.par_nuc;
                uniq_mut.par_nuc = uniq_mut.mut_nuc;
                uniq_mut.mut_nuc = temp;
                curr_node_par_mut.emplace_back(uniq_mut);
            }
        }

        //Updating curr_node_mut for having read as child
        curr_par.p_node_par = curr_node_par_mut;
        curr_par.curr_node = dfs[i];
        parsimony_stack.push(curr_par);
        
        uniq_curr_node_mut.clear();
        common_node_mut.clear();
        curr_node_par_mut.clear();
    }

    //Clear Stack
    while (!parsimony_stack.empty())
        parsimony_stack.pop();
    
    //Compute node_score
    if (check_node == NULL) {
        min_par.par_list.clear();
        //Keeping tab on weighted read score being mapped to each parsimonious node
        //Give score to parsimonious node: Score -> 1/N^2
        double score = (1.0 / pow(min_par.idx_list.size(), 2));
        for (auto n_idx: min_par.idx_list) {
            auto node = dfs[n_idx];
            tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
            auto created = node_score.insert(ac, std::make_pair(node, score));
            if (!created)
                ac->second += score;
            ac.release();
        }
        min_par.idx_list.clear();
        return 0;
    }
    //Check if given node is present in parsimonious list of read
    else {
        min_par.par_list.clear();
        for (auto n_idx: min_par.idx_list) {
            auto node = dfs[n_idx];
            if (node == check_node) {
                min_par.idx_list.clear();
                return 1;
            }
        }
        min_par.idx_list.clear();
        return 0;
    }
}

//Main peak search algorithm
void analyze_reads(const MAT::Tree &T, const MAT::Tree &T_ref, const std::unordered_map<size_t, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score, const std::vector<std::string> &vcf_samples, const std::string &barcode_file, const std::string &read_mutation_depth_vcf) {
    timer.Start();
    int top_n = 25, m_dist_thresh = 2, neighbor_dist_thresh = 7, neighbor_peaks_thresh = 100;
    std::vector<MAT::Node*> dfs, peak_nodes, curr_peak_nodes, prohibited_nodes;
    tbb::concurrent_vector<size_t> remaining_reads(read_map.size());
    
    //Depth first expansion to get all nodes in the tree and 
    dfs = T.depth_first_expansion(); 
    
    //GREEDY ALGORITHM for peak nodes
    /*
        1. Place remaining reads on tree
        2. Remove nodes not to be considered (present in Prohibited nodes)
        3. Find top scoring nodes
        4. Only consider highest scoring node(s) not seen before
        5. Remove reads mapped to these nodes
        6. Don't consider top scoring nodes in the neighborhood of other top nodes
        7. Add neighbors of selected top scoring nodes to Prohibited nodes 
    */
    //Initializing remaining_reads
    for (size_t i = 0; i < read_map.size(); i++)
        remaining_reads[i] = i;
    //Iterate till no reads left
    while ((int)remaining_reads.size() > 0) {
        //Calculating node score for remaining reads
        printf("\n");
        static tbb::affinity_partitioner ap;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    size_t rm_idx = remaining_reads[i];
                    auto read_id = read_map.find(rm_idx)->second;
                    place_reads(dfs, read_id, NULL, node_score);
                }
            },
        ap);
        //REMOVE prohibited_nodes from node_score
        tbb::parallel_for( tbb::blocked_range<size_t>(0, prohibited_nodes.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto p_node = prohibited_nodes[i];
                    tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
                    if (node_score.find(ac, p_node))
                        node_score.erase(ac);
                }
            },
        ap);
        //Get top_n node_scores
        std::vector<std::pair<MAT::Node*, double>> top_n_node_score(top_n);
        std::partial_sort_copy(node_score.begin(),
                            node_score.end(),
                            top_n_node_score.begin(),
                            top_n_node_score.end(),
                            [](std::pair<const MAT::Node*, double> const& l,
                               std::pair<const MAT::Node*, double> const& r)
                            {
                                return l.second > r.second;
                            });
        node_score.clear();

        //Only consider nodes not seen before
        std::vector<bool> peak_vec(top_n_node_score.size(), true);
        //Find the reads not mapped to the best nodes
        auto top_n_itr = top_n_node_score.begin();
        auto top_score = top_n_itr->second;
        while (top_n_itr != top_n_node_score.end()) {
            auto curr_node = top_n_itr->first;
            //Only add unique nodes in curr_peak_nodes with score == top_score
            if (abs(top_score - top_n_itr->second) < 1e-9) {
                //Check if current node is not in neighborhood of peaks seen till now
                if (peak_vec[top_n_itr - top_n_node_score.begin()]) {
                    curr_peak_nodes.emplace_back(curr_node);
                    auto curr_clade = get_clade(T, curr_node);
                    printf("PEAK: %s, %s\n",curr_node->identifier.c_str(), curr_clade.c_str());
                    //Remove reads mapped to curr_node
                    std::vector<size_t> remove_reads;
                    using my_mutex_t = tbb::queuing_mutex;
                    my_mutex_t my_mutex;
                    static tbb::affinity_partitioner ap;
                    //Parallel_for loop for each remaining read
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
                        [&](tbb::blocked_range<size_t> k) {
                            for (size_t i = k.begin(); i < k.end(); ++i) {
                                size_t rm_idx = remaining_reads[i];
                                auto read_id = read_map.find(rm_idx)->second;
                                int read_present = place_reads(dfs, read_id, curr_node, node_score);
                                node_score.clear();
                                //If Read contains current node -> REMOVE
                                if (read_present) {
                                    my_mutex_t::scoped_lock my_lock{my_mutex};
                                    remove_reads.emplace_back(rm_idx);
                                }
                            }
                        },
                    ap);
                    //Remove reads from remaining_reads present in remove_reads
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, remove_reads.size()),
                        [&](tbb::blocked_range<size_t> k) {
                            for (size_t i = k.begin(); i < k.end(); ++i) {
                                const size_t rm_idx = remove_reads[i];
                                // Use std::remove_if and lambda function to remove the element
                                remaining_reads.erase(std::remove_if(remaining_reads.begin(), remaining_reads.end(),
                                                           [&](const size_t& element) {
                                                               return element == rm_idx;
                                                           }),
                                                        remaining_reads.end());
                            }
                        },
                    ap);
                    remove_reads.clear();
                }
                //If present in neighbourhood of current peaks, move to next node
                else {
                    top_n_itr++;
                    continue;
                }
            }
            //node score is less than top score
            else
                break;
            
            //Remove nearby peaks from further analysis
            auto top_n_peak_cmp_itr = top_n_itr;
            std::advance(top_n_peak_cmp_itr, 1);
            while (top_n_peak_cmp_itr != top_n_node_score.end()) {
                auto cmp_node = top_n_peak_cmp_itr->first;
                //Stop checking if score < top_score
                if (abs(top_score - top_n_peak_cmp_itr->second) > 1e-9)
                    break;
                //Check next peak if this peak is already in neighborhood of some other peak
                else if (!peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()]) {
                    top_n_peak_cmp_itr++;
                    continue;
                }
                //Don't consider peaks within mutation distance limit 
                int m_dist = mutation_distance(T, T, curr_node, cmp_node);
                if (m_dist <= m_dist_thresh)
                    peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = false;
                top_n_peak_cmp_itr++;
            }

            top_n_itr++;
        }
        top_n_node_score.clear();
        peak_vec.clear();

        //ADDING current_peak_nodes and its neighbors to prohibited_nodes
        /*
            1. Find the farthest ancestor of every peak node within m_dist_thresh
            2. Recursively only analyze its children within m_dist_thresh
            3. Only include unique neighborhood nodes to prohibited_nodes
        */
        using my_mutex_t = tbb::queuing_mutex;
        my_mutex_t my_mutex;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, curr_peak_nodes.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto peak = curr_peak_nodes[i];
                    //Find farthest ancestor with mut distance from peak <= m_dist_thresh
                    MAT::Node* anc_node;
                    for (const auto& n: T.rsearch(peak->identifier, true)) {
                        int m_dist = mutation_distance(T, T, n, peak);
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
                        int m_dist = mutation_distance(T, T, present_node, peak);
                        if (m_dist <= m_dist_thresh) {
                            //Only add if present_node NOT already present in prohibited_nodes
                            //Looking for exact node matches here as we want to remove same mutation nodes
                            my_mutex_t::scoped_lock my_lock;
                            my_lock.acquire(my_mutex);
                            if (!std::binary_search(prohibited_nodes.begin(), prohibited_nodes.end(), present_node))
                                prohibited_nodes.emplace_back(present_node)
                            my_lock.release();
                            //Add present_node's children in list for checking
                            for (const auto& c: present_node->children)
                                remaining_nodes.push(c);
                        }
                    }
                }
            },
        ap);

        //Add curr_peak_nodes to peak_nodes
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
        MAT::Node* best_node;
        for (auto pn: peak_nodes) {
            int curr_dist = mutation_distance(T, T_ref, pn, sample_node);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = pn;
            }
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld msec\n\n", timer.Stop());

    add_neighbor_peaks(T, peak_nodes, neighbor_dist_thresh, neighbor_peaks_thresh);
    generate_regression_abundance_data(T, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf);
}


//Add neighboring nodes to peaks
void add_neighbor_peaks(const MAT::Tree &T, std::vector<MAT::Node*> &peak_nodes, const int neighbor_dist_thresh, const int neighbor_peaks_thresh) {
    std::vector<MAT::Node*> neighbor_nodes;
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, peak_nodes.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                MAT::Node* peak = peak_nodes[i];
                //Keep looking for neighbors until neighbors added == neighbor_dist_thresh OR mut_distance > neighbor_dist_thresh
                int neighbor_added = 0, branch_length_limit = 1;
                std::vector<MAT::Node*> prev_visited_nodes, curr_peak_neighbors;
                while (neighbor_added < neighbor_peaks_thresh) {
                    // 1. Search children nodes within branch_length_limit
                    auto curr_node = peak;
                    int curr_branch_length = 0;
                    child_nodes_addition(my_mutex, peak_nodes, neighbor_nodes, neighbor_added, prev_visited_nodes, peak, curr_node, curr_branch_length, branch_length_limit, neighbor_peaks_thresh, neighbor_dist_thresh);
                    // 2. Search parent nodes within branch_length_limit
                    while ((curr_branch_length < branch_length_limit) && (neighbor_added < neighbor_peaks_thresh)) {
                        curr_node = curr_node->parent;
                        curr_branch_length++;
                        child_nodes_addition(my_mutex, peak_nodes, neighbor_nodes, neighbor_added, prev_visited_nodes, peak, curr_node, curr_branch_length, branch_length_limit, neighbor_peaks_thresh, neighbor_dist_thresh);
                    }
                    branch_length_limit++;
                }




                ///////////////////////////////////////////////////////////////////////////////////////////CHANGE
                //Find farthest ancestor with mut distance from peak <= neighbor_dist_thresh
                MAT::Node* anc_node;
                for (auto n: T.rsearch(peak->identifier, true)) {
                    //Get mutation distance between ancestor and curr_node
                    int m_dist = mutation_distance(T, T, n, peak);
                    if (m_dist <= neighbor_dist_thresh)
                        anc_node = n;
                    else
                        break;
                }
                //Add neighborhood peaks to peak_nodes
                std::queue<MAT::Node*> remaining_nodes;
                remaining_nodes.push(anc_node);
                int neighbor_added = 0;
                while (remaining_nodes.size() > 0) {
                    MAT::Node* present_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    //Only add present_node if within neighbor_dist_thresh
                    int m_dist = mutation_distance(T, T, present_node, peak);
                    if ((m_dist <= neighbor_dist_thresh) && (neighbor_added < neighbor_peaks_thresh)) {
                        bool found = false;
                        //Check if a peak with same mutations is NOT already present
                        //Not going for exact node match as we want single representative of mutation vector
                        for (const auto& hap: peak_nodes) {
                            int mut_dist = mutation_distance(T, T, hap, present_node);
                            if (!mut_dist) {
                                found = true;
                                break;
                            }
                        }
                        //Check if neighbor peak with same mutations is already accounted
                        if (!found) {
                            my_mutex_t::scoped_lock my_lock;
                            my_lock.acquire(my_mutex);
                            for (const auto& hap: neighbor_nodes) {
                                int mut_dist = mutation_distance(T, T, hap, present_node);
                                if (!mut_dist) {
                                    found = true;
                                    break;
                                }
                            }
                            my_lock.release();
                        }
                        //Only add if similar node hasn't been accounted
                        if (!found) {
                            my_mutex_t::scoped_lock my_lock;
                            my_lock.acquire(my_mutex);
                            neighbor_nodes.emplace_back(present_node);
                            my_lock.release();
                            neighbor_added++;
                        }
                        //Add current node's children in list for checking
                        for (auto c: present_node->children)
                            remaining_nodes.push(c);
                    }
                    //Stop checking for neighbors if neighbor_added == neighbor_peak_thresh
                    else if (neighbor_added == neighbor_peaks_thresh) {
                        while (!remaining_nodes.empty())
                            remaining_nodes.pop();
                    }
                }

                ///////////////////////////////////////////////////////////////////////////////////////////CHANGE
            }
        },
    ap);
        //Reserve more memory to ensure there are no segmentation faults
        peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
        peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
        neighbor_nodes.clear();
        printf("\nPost neighbor adding PEAK nodes: %d\n", (int)peak_nodes.size());
}

//Function to search for child nodes within given branch length
void child_nodes_addition(my_mutex_t& my_mutex, const std::vector<MAT::Node*>& peak_nodes, std::vector<MAT::Node*>& neighbor_nodes, int& neighbor_added, std::vector<MAT::Node*>& prev_visited_nodes, const MAT::Node* peak, const MAT::Node* init_node, const int& init_branch_length, const int& branch_length_thresh, const int& neighbor_peaks_thresh, const int& neighbor_dist_thresh) {
    int curr_branch_length = init_branch_length;
    auto curr_node = init_node;
    std::vector<MAT::Node*> selected_neighbors;
    std::queue<struct node_branch> remaining_nodes_branch;
    remaining_nodes_branch.push(std::make_pair(curr_node, curr_branch_length));
    
    //Conditions to check: neighbor_peaks_thresh, neighbor_dist_thresh, prev_visited_nodes, neighbor_nodes, branch_length_thresh 
    while(remaining_nodes_branch.size() > 0) {
        struct node_branch curr_node_branch = remaining_nodes_branch.front();
        remaining_nodes_branch.pop();
        //Only analyze if neighbor_peaks_thresh NOT reached
        if (neighbor_added < neighbor_peaks_thresh) {
            int m_dist = mutation_distance(T, T, curr_node_branch.node, peak);
            //Only analyze nodes and their children within neighbor_dist_thresh
            if (m_dist <= neighbor_dist_thresh) {
                //Do further checks on this node if NOT analyzed before
                if (std::find(prev_visited_nodes.begin(), prev_visited_nodes.end(), curr_node_branch.node) == prev_visited_nodes.end()) {
                    prev_visited_nodes.emplace_back(curr_node_branch.node);
                    bool found = false;
                    //Check if a peak with same mutations is NOT already present
                    //Not going for exact node match as we want single representative of mutation vector
                    for (const auto& hap: peak_nodes) {
                        int mut_dist = mutation_distance(T, T, hap, curr_node_branch.node);
                        if (!mut_dist) {
                            found = true;
                            break;
                        }
                    }
                    //Check if neighbor peak with same mutations is already accounted
                    if (!found) {
                        my_mutex_t::scoped_lock my_lock{my_mutex};
                        for (const auto& hap: neighbor_nodes) {
                            int mut_dist = mutation_distance(T, T, hap, curr_node_branch.node);
                            if (!mut_dist) {
                                found = true;
                                break;
                            }
                        }
                    }
                    //Only add if similar node hasn't been accounted
                    if (!found) {
                        my_mutex_t::scoped_lock my_lock;
                        my_lock.acquire(my_mutex);
                        neighbor_nodes.emplace_back(present_node);
                        my_lock.release();
                        neighbor_added++;
                    }
                }
                
                //Only add current node's children in list if within branch_length_thresh
                curr_branch_length = curr_node_branch.branch_length + 1;
                if (curr_branch_length <= branch_length_thresh) {
                    for (auto child: present_node->children)
                        remaining_nodes_branch.push(std::make_pair(child, curr_branch_length));
                }
            }
        }
        //Stop checking for neighbors if neighbor_added == neighbor_peak_thresh
        else {
            while (!remaining_nodes_branch.empty())
                remaining_nodes_branch.pop();
        }
    }
}

//Function to calculation distance between two nodes
int mutation_distance(const MAT::Tree &T1, const MAT::Tree &T2, const MAT::Node* N1, const MAT::Node* N2) {
    std::vector<MAT::Mutation> node1_mutations, node2_mutations;
    //Checking all ancestors of a node
    for (auto anc: T1.rsearch(N1->identifier, true)) { 
        for (auto mut: anc->mutations) {
            //If mutation at same position is seen closer to leaf then consider that mutation
            auto n1_itr = node1_mutations.begin();
            while (n1_itr != node1_mutations.end()) {
                if (n1_itr->position == mut.position)
                    break;
                n1_itr++;
            }
            if (n1_itr == node1_mutations.end())
                node1_mutations.emplace_back(mut);
        }
    }
    //Checking all ancestors of a node
    for (auto anc: T2.rsearch(N2->identifier, true)) { 
        for (auto mut: anc->mutations) {
            //If mutation at same position is seen closer to leaf then consider that mutation
            auto n2_itr = node2_mutations.begin();
            while (n2_itr != node2_mutations.end()) {
                if (n2_itr->position == mut.position)
                    break;
                n2_itr++;
            }
            if (n2_itr == node2_mutations.end())
                node2_mutations.emplace_back(mut);
        }
    }
    //Remove Back Mutations -> mutations with same ref_nuc and mut_nuc
    auto n1_itr = node1_mutations.begin();
    while (n1_itr != node1_mutations.end()) {
        if (n1_itr->ref_nuc == n1_itr->mut_nuc)
            n1_itr = node1_mutations.erase(n1_itr);
        else
            n1_itr++;
    }
    auto n2_itr = node2_mutations.begin();
    while (n2_itr != node2_mutations.end()) {
        if (n2_itr->ref_nuc == n2_itr->mut_nuc)
            n2_itr = node2_mutations.erase(n2_itr);
        else
            n2_itr++;
    }
    //Sort the mutations to find difference between node1_mutations and node2_mutations faster
    std::sort(node1_mutations.begin(), node1_mutations.end(), compare_mutations);
    std::sort(node2_mutations.begin(), node2_mutations.end(), compare_mutations);
    //Finding the unique mutations between them
    n1_itr = node1_mutations.begin();
    while (n1_itr != node1_mutations.end()) {
        bool found_both = false;
        n2_itr = node2_mutations.begin();
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

//Get the clade name
std::string get_clade(const MAT::Tree &T, MAT::Node* n) {
    //Checking all ancestors of a node to get clade
    for (auto anc: T.rsearch(n->identifier, true)) {  
        //Don't consider the clades that start with numbers {Different clade naming}
        auto clade = anc->clade_annotations[1];
        if(clade != "")
            return clade;
    }
    return "";
}

//Generate regression based estimate algorithm data
void generate_regression_abundance_data(const MAT::Tree &T, const std::vector<MAT::Node*> &peak_nodes, const std::unordered_map<size_t, struct read_info*> &read_map, const std::string &barcode_file, const std::string &read_mutation_depth_vcf) {
    timer.Start();
    //Get Mutations of Peaks for deconvolution method
    //Create a map of peak node mutations 
    tbb::concurrent_hash_map<MAT::Node*, std::vector<MAT::Mutation>> peak_mut_map;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, peak_nodes.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto pk = peak_nodes[i];
                    std::vector<MAT::Mutation> mut_list;
                    for (auto anc: T.rsearch(pk->identifier, true)) {
                        for (auto c_mut: anc->mutations) {
                            //Only consider the mutation closest to root at a site
                            auto m_itr = mut_list.begin();
                            while (m_itr != mut_list.end()) {
                                if (m_itr->position == c_mut.position)
                                    break;
                                m_itr++;
                            }
                            if (m_itr == mut_list.end())
                                mut_list.emplace_back(c_mut);
                        }
                    }
                    //Remove mutations where ref_nuc is same as mut_nuc -> Back Mutation
                    auto m_itr = mut_list.begin();
                    while (m_itr != mut_list.end()) {
                        if (m_itr->mut_nuc == m_itr->ref_nuc)
                            m_itr = mut_list.erase(m_itr);
                        else
                            m_itr++;
                    }
                    tbb::concurrent_hash_map<MAT::Node*, std::vector<MAT::Mutation>>::accessor ac;
                    peak_mut_map.insert(ac, std::make_pair(pk, mut_list));
                    ac.release();
                    mut_list.clear();
                }
            },
        ap);

    //REMOVE the mutations that are common to all peak nodes
    auto ref_pk_itr = peak_mut_map.begin();
    //Iterating though mutations of ref_peak;
    auto ref_mut_itr = ref_pk_itr->second.begin();
    while (ref_mut_itr != ref_pk_itr->second.end()) {
        int common_mut_count = 1;
        //Iterating through rest of peak nodes
        auto peak_itr = std::next(ref_pk_itr,1);
        while (peak_itr != peak_mut_map.end()) {
            //Iterating though mutations of rest of peaks;
            auto cmp_mut_itr = peak_itr->second.begin();
            while (cmp_mut_itr != peak_itr->second.end()) {
                if ((cmp_mut_itr->position == ref_mut_itr->position) && (cmp_mut_itr->mut_nuc == ref_mut_itr->mut_nuc)) {
                    common_mut_count++;
                    break;
                }
                cmp_mut_itr++;
            }
            peak_itr++;
        }
        //Remove common mutation across all peaks
        if (common_mut_count == (int)peak_mut_map.size()) {
            //Iterating through rest of peak nodes
            peak_itr = std::next(ref_pk_itr,1);
            while (peak_itr != peak_mut_map.end()) {
                //Iterating though mutations of rest of peaks;
                auto cmp_mut_itr = peak_itr->second.begin();
                while (cmp_mut_itr != peak_itr->second.end()) {
                    if ((cmp_mut_itr->position == ref_mut_itr->position) && (cmp_mut_itr->mut_nuc == ref_mut_itr->mut_nuc)) {
                        peak_itr->second.erase(cmp_mut_itr);
                        break;
                    }
                    cmp_mut_itr++;    
                }
                peak_itr++;
            }
            //Erase mutation in ref peak as well
            ref_mut_itr = ref_pk_itr->second.erase(ref_mut_itr);
        }
        else 
            ref_mut_itr++;
    }

    //Initialize the buffers for writing files
    std::ofstream outfile_barcode(barcode_file, std::ios::out | std::ios::binary);
    std::ofstream outfile_vcf(read_mutation_depth_vcf, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_barcode, outbuf_vcf;
    if (barcode_file.find(".gz\0") != std::string::npos) {
            outbuf_barcode.push(boost::iostreams::gzip_compressor());
    }
    if (read_mutation_depth_vcf.find(".gz\0") != std::string::npos) {
            outbuf_vcf.push(boost::iostreams::gzip_compressor());
    }
    outbuf_barcode.push(outfile_barcode);
    outbuf_vcf.push(outfile_vcf);
    std::ostream barcode(&outbuf_barcode);
    std::ostream vcf(&outbuf_vcf);

    //ADD unique mutations captured from read_vcf first
    std::vector<MAT::Mutation> peak_mut_list;
    auto rm_itr = read_map.begin();
    while (rm_itr != read_map.end()) {
        for (auto read_mut: rm_itr->second->mutations) {
            auto pm_itr = peak_mut_list.begin();
            while (pm_itr != peak_mut_list.end()) {
                if ((pm_itr->position == read_mut.position) && (pm_itr->mut_nuc == read_mut.mut_nuc))
                    break;
                pm_itr++;
            }
            if (pm_itr == peak_mut_list.end())
                peak_mut_list.emplace_back(read_mut);
        }
        rm_itr++;
    }
    //Storing unique mutations from the peak nodes
    auto peak_mut_itr = peak_mut_map.begin();
    while (peak_mut_itr != peak_mut_map.end()) {
        //Iterating through mutations of current peak
        auto mut_itr = peak_mut_itr->second.begin();
        while (mut_itr != peak_mut_itr->second.end()) {
            //Only add unique mutations to the list
            auto pm_itr = peak_mut_list.begin();
            while (pm_itr != peak_mut_list.end()) {
                if ((pm_itr->position == mut_itr->position) && (pm_itr->mut_nuc == mut_itr->mut_nuc))
                    break;
                pm_itr++;                
            }
            if (pm_itr == peak_mut_list.end())
                peak_mut_list.emplace_back(*mut_itr);
            mut_itr++;
        }
        peak_mut_itr++;
    }

    //VCF needs mutations in sorted order
    std::sort(peak_mut_list.begin(), peak_mut_list.end(), compare_mutations);
    //WRITING the header mutations in barcode file and complete VCF
    std::string barcode_print, vcf_print;
    vcf_print += "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tDEPTH\n";
    for (auto mut: peak_mut_list) {
        //Write the mutations header in barcode file
        barcode_print += ",";
        barcode_print += MAT::get_nuc(mut.par_nuc) + std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc);
        //Calculate AF
        int total_reads = 0;
        int match_reads = 0;
        auto rm_itr = read_map.begin();
        while (rm_itr != read_map.end()) {
            if ((mut.position >= rm_itr->second->start) && (mut.position <= rm_itr->second->end)) {
                total_reads++;
                for (auto read_mut: rm_itr->second->mutations) {
                    if ((read_mut.position == mut.position) && (read_mut.mut_nuc == mut.mut_nuc)) {
                        match_reads++;
                        break;
                    }
                }
            }
            rm_itr++;
        }
        float af = 0.0;
        if (total_reads > 0)
            af = (float)match_reads / (float)total_reads;
        vcf_print += "NC_045512v2\t" + std::to_string(mut.position) + "\t" + MAT::get_nuc(mut.par_nuc) + std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc) + "\t" + MAT::get_nuc(mut.par_nuc) + "\t" + MAT::get_nuc(mut.mut_nuc) + "\t.\t.\tAF=";
        vcf_print += std::to_string(af) + "\t" + std::to_string(total_reads) + "\n";
    } 
    vcf << vcf_print;
    vcf_print.clear();

    //WRITING peak mutations in barcode
    peak_mut_itr = peak_mut_map.begin();
    while (peak_mut_itr != peak_mut_map.end()) {
        //Write peak name in barcode file
        barcode_print += "\n";
        barcode_print += peak_mut_itr->first->identifier + "_" + get_clade(T, peak_mut_itr->first);
        //Iterating through mutations of this peak to store indices its mutations from peak_mut_list
        std::vector<size_t> mut_idx_list;
        auto mut_itr = peak_mut_itr->second.begin();
        while (mut_itr != peak_mut_itr->second.end()) {
            //Search for mutation in the peak_mut_list
            auto pm_itr = peak_mut_list.begin();
            while (pm_itr != peak_mut_list.end()) {
                if ((pm_itr->position == mut_itr->position) && (pm_itr->mut_nuc == mut_itr->mut_nuc))
                    break;
                pm_itr++;                
            }
            mut_idx_list.emplace_back(pm_itr - peak_mut_list.begin());
            mut_itr++;
        }
        //Sort the mutation indexes before writing to barcode
        std::sort(mut_idx_list.begin(), mut_idx_list.end());
        
        //Write peak mutations in barcode file
        size_t idx = 0;
        for (const auto& m_idx: mut_idx_list) {
            while (idx < m_idx) {
                barcode_print += ",0";
                idx++;
            }
            if (idx == m_idx) {
                barcode_print += ",1";
                idx++;
            }
        }
        while (idx < peak_mut_list.size()) {
            barcode_print += ",0";
            idx++;
        }
        mut_idx_list.clear();
        peak_mut_itr++;
    }
    barcode << barcode_print;
    barcode_print.clear();
        
    fprintf(stderr,"Barcode file writing took %ld sec\n\n", (timer.Stop() / 1000));
}

//Comparing mutations for sorting a vector 
bool compare_mutations(const MAT::Mutation &a, const MAT::Mutation &b) {
    return a.position < b.position;
}