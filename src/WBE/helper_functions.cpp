#include "wbe.hpp"
//Paring the commands
po::variables_map parseWBEcommand(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("Given Switch options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->default_value(""),
     "Input mutation-annotated tree file")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("output-files-prefix,v", po::value<std::string>()->default_value(""),
    "Prefix to be used for dumping all intermediate files.")
    ("ref-fasta,f", po::value<std::string>()->default_value(""),
     "Input fasta file representing reference sequence")
    ("align-sam,s", po::value<std::string>()->default_value(""),
     "Input sam file representing reference sequence")
    ("distribution,d", po::value<std::string>()->default_value(""),
     "Give the distribution of samples, comma delimited.")
    ("haplotype-samples,w", po::value<int>()->default_value(10),
     "Give the number of haplotype samples")
    ("lineage,l", po::value<std::string>()->default_value(""),
     "Give lineage of samples, comma delimited.")
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

//Store the Reads from the VCF in read_map
void readVCF(std::unordered_map<size_t, struct read_info*> &read_map, const std::string &vcf_filename_reads, const size_t &seq_len, const bool &filter) {
    // Boost library used to stream the contents of the input VCF file
    tbb::concurrent_hash_map<int8_t, std::vector<size_t>> mut_read_map;
    double ignore_thresh = 0.005;

    timer.Start();
    std::string s;
    boost::filesystem::ifstream fileHandler(vcf_filename_reads);
    bool header_found = false;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            size_t field_offset = 9;
            //Checking for header
            if (words[1] == "POS") {
                header_found = true;
                //Leave certain fields based on our VCF format
                for (size_t j = field_offset; j < words.size(); j++) {
                    struct read_info* rp = new struct read_info;
                    rp->read = words[j];
                    //Get start-end positions and depth of the Read
                    std::regex rgx(".*READ_(\\w+)_(\\w+)");
                    std::smatch match;

                    if (std::regex_search(words[j], match, rgx)) {
                        rp->start = std::stoi(match[1]);
                        rp->end = std::stoi(match[2]);
                    }
                    //Consider the entire genome
                    else {
                        rp->start = 1;
                        rp->end = seq_len;
                    }
                    read_map.insert({(j-field_offset), rp});
                }
            }
            else if (header_found) {
                std::vector<std::string> alleles;
                //Checking for different alleles at a site
                MAT::string_split(words[4], ',', alleles);
                static tbb::affinity_partitioner ap;
                tbb::parallel_for(tbb::blocked_range<size_t>(0, read_map.size()),
                    [&](tbb::blocked_range<size_t> i) {
                        for (size_t k = i.begin(); k < i.end(); ++k) {
                            bool mut_present = false;
                            size_t j = k + field_offset;
                            auto rp = read_map[k];
                            MAT::Mutation m;
                            m.chrom = words[0];
                            m.position = std::stoi(words[1]);
                            //Checking the mutating allele value within the allele sizes
                            if (std::stoi(words[j]) > int(alleles.size())) {
                                fprintf(stderr, "\n\nPosition: %d, k = %ld,\n", m.position, k);
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
                                    mut_present = true;
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
                                    rp->mutations.emplace_back(m);
                                }
                                else if ((!allele_id) && (m.position >= rp->start) && (m.position <= rp->end)) {
                                    m.mut_nuc = m.ref_nuc;
                                    mut_present = true;
                                }

                            } else {
                                mut_present = true;
                                m.is_missing = true;
                                m.mut_nuc = MAT::get_nuc_id('N');
                                rp->mutations.emplace_back(m);
                            }
                            
                            //STORE the number of reads for each nucleotide
                            if (mut_present && filter) {
                                tbb::concurrent_hash_map<int8_t, std::vector<size_t>>::accessor ac;
                                auto created = mut_read_map.insert(ac, {m.mut_nuc, {k}});
                                if (!created) 
                                    ac->second.emplace_back(k);
                                ac.release();
                            }
                        }
                    },
                ap);

                //UPDATE mutations with read_coverage < ignore_thresh 
                size_t total_mut_reads = 0, max_mut_reads = 0;
                int8_t max_nuc = 0b1111;
                //Get total_mut_reads and max_mut
                for (const auto& mut_reads: mut_read_map) {
                    auto curr_reads = mut_reads.second.size();
                    total_mut_reads += curr_reads;
                    if (curr_reads > max_mut_reads) {
                        max_mut_reads = curr_reads;
                        max_nuc = mut_reads.first;
                    }
                }
                //Update mutation on reads
                for (const auto& mut_reads: mut_read_map) {
                    auto curr_reads = mut_reads.second.size();
                    double coverage_frac = static_cast<double>(curr_reads) / total_mut_reads;
                    if ((abs(coverage_frac - ignore_thresh) < 1e-9) || (coverage_frac < ignore_thresh)) {
                        /*
                            1. If max_nuc is ref_nuc then remove mutation from reads
                            2. If curr_nuc is ref_nuc then introduce mutation in reads
                            3. Replace mut_nuc in reads with max_nuc
                        */
                        int8_t ref_nuc = MAT::get_nuc_id(words[3][0]);
                        for (const auto& rd_idx: mut_reads.second) {
                            auto rp = read_map[rd_idx];
                            //Remove last added mutation
                            if (ref_nuc == max_nuc)
                                rp->mutations.pop_back();
                            //Introduce mutation
                            else if (ref_nuc == mut_reads.first) {
                                MAT::Mutation m;
                                m.chrom = words[0];
                                m.position = std::stoi(words[1]);
                                m.ref_nuc = ref_nuc;
                                m.par_nuc = m.ref_nuc;
                                m.is_missing = false;
                                m.mut_nuc = max_nuc;
                                rp->mutations.emplace_back(m);
                            }
                            //Replace mut_nuc
                            else {
                                auto &m = rp->mutations.back();
                                m.is_missing = false;
                                m.mut_nuc = max_nuc; 
                            }
                        }
                    }
                }
                alleles.clear();
                mut_read_map.clear();
            }
        }
    }
    
    fprintf(stderr,"%s parsed in %ld sec\n\n", vcf_filename_reads.c_str(), (timer.Stop() / 1000));   
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

//Function to read csv containing condensed nodes
void readCSV(std::unordered_map<std::string, std::vector<std::string>>& condensed_nodes_map, const std::string &condensed_nodes_csv) {
    std::ifstream file(condensed_nodes_csv);
    if (!file.is_open()) {
        std::cout << "Failed to open the csv file" << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> node_name_list;
        std::string key;
        if (std::getline(iss, key, ',')) {
            std::string node_name;
            while (std::getline(iss, node_name, ','))
                node_name_list.emplace_back(node_name);
        }
        condensed_nodes_map[key] = node_name_list;
        node_name_list.clear();
    }
    file.close();
}

//PLACING reads using small range trees
void placeReadHelper(MAT::Node* ref_root, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, const std::unordered_map<size_t, struct read_info*>& read_map, std::vector<size_t> remaining_read_ids, const std::vector<MAT::Node*>& peak_nodes, tbb::concurrent_hash_map<MAT::Node*, double>& node_score_map, std::vector<size_t>& remove_reads, const int& seq_len, const int& tree_increment, const int& tree_range) {
    static tbb::affinity_partitioner ap;
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
    //MAPPING reads to nodes
    for (int start = 0; start < seq_len; start += tree_increment) {
        int end = start + tree_range;
        std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> node_mappings;
        MAT::Tree range_Tree;
        std::vector<size_t> curr_read_ids;
        
        //Extract all the remaining_read_ids within range
        for (size_t i = 0; i < remaining_read_ids.size(); i++) {
            auto rm_idx = remaining_read_ids[i];
            auto read_id = read_map.find(rm_idx)->second;	
            if ((read_id->start >= start) && (read_id->end <= end))
                curr_read_ids.emplace_back(rm_idx);
        }
        //Sort curr_read_ids
        tbb::parallel_sort(curr_read_ids.begin(), curr_read_ids.end());

        //REMOVE curr_read_ids from remaining_read_ids
        auto rri_itr = remaining_read_ids.begin();
        auto cri_itr = curr_read_ids.begin();
        while ((rri_itr != remaining_read_ids.end()) && (cri_itr != curr_read_ids.end())) {
            //Remove from remaining_read_ids if equal
            if (*rri_itr == *cri_itr) {
                cri_itr++;
                rri_itr = remaining_read_ids.erase(rri_itr);
            }
            //Move to next curr_read_ids if current curr_read_ids not found in remaining_read_ids
            else if (*rri_itr > *cri_itr) {
                cri_itr++;
                fprintf(stderr,"curr_read_ids not present in remaining_read_ids!!!");
            }
            //Move to next remaining_read_ids if it is smaller than current curr_read_ids
            else
                rri_itr++;
        }
        //Create range tree
        if (!curr_read_ids.empty())
            createRangeTree(ref_root, condensed_node_mappings, start, end, node_mappings, range_Tree);

        //Tree Search
        tbb::parallel_for(tbb::blocked_range<size_t>(0, curr_read_ids.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    size_t rm_idx = curr_read_ids[i];
                    auto read_id = read_map.find(rm_idx)->second;
                    //MAP reads to nodes
                    if (peak_nodes.empty())	
                        placeReads(range_Tree, read_id, peak_nodes, node_score_map, node_mappings, condensed_node_mappings);
                    //GET list of reads mapped to peak_nodes
                    else {
                        int read_present = placeReads(range_Tree, read_id, peak_nodes, node_score_map, node_mappings, condensed_node_mappings);
                        //If Read contains current node -> REMOVE
                        if (read_present) {
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            remove_reads.emplace_back(rm_idx);
                        }
                    }
                }
            },
        ap);

        //Clear data structures
        node_mappings.clear();
        MAT::clear_tree(range_Tree);
        curr_read_ids.clear();
        //Do not check remaining range trees if no reads are left
        if (remaining_read_ids.empty())
            break;
    }
    if (!remaining_read_ids.empty())
       fprintf(stderr,"Read Left: %lu\nSmaller MAT Code NOT working properly!!!\n\n", remaining_read_ids.size());
}

//Parsimonious placement search for reads
int placeReads(const MAT::Tree &T, const struct read_info* rp, const std::vector<MAT::Node*> &check_nodes, tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings) {
    auto dfs = T.depth_first_expansion();
    std::stack<struct parsimony> parsimony_stack;
    struct min_parsimony min_par;
    std::unordered_set<size_t> common_mut_nodes;
    //Checking all nodes of the tree
    for (size_t i = 0; i < dfs.size(); i++) {
        auto curr_node = dfs[i];
        std::vector<MAT::Mutation> uniq_curr_node_mut, common_node_mut, curr_node_par_mut;
        struct parsimony curr_par;
        //Nothing in stack for first node
        if (i) {
            //Get the parsimony vector from parent
            auto parent_parsimony = parsimony_stack.top();
            while ((curr_node->parent != parent_parsimony.curr_node) && (!parsimony_stack.empty()))   {
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
        for (auto node_mut: curr_node->mutations) {
            //Only look at mutations within the read range
            if ((node_mut.position >= rp->start) && (node_mut.position <= rp->end)) {
                bool found = false;
                //Check mutation position in parsimony of parent node
                auto curr_node_par_itr = curr_node_par_mut.begin();
                while (curr_node_par_itr != curr_node_par_mut.end()) {
                    if (curr_node_par_itr->position == node_mut.position) {
                        //mut_nuc matches => remove mutation from parsimony
                        if (curr_node_par_itr->mut_nuc == node_mut.mut_nuc) {
                           common_node_mut.emplace_back(*curr_node_par_itr);
                           curr_node_par_mut.erase(curr_node_par_itr);
                        }
                        //Update the par_nuc in parsimony if only position matches
                        else
                            curr_node_par_itr->par_nuc = node_mut.mut_nuc;
                        found = true;
                        break;
                    }
                    else 
                        curr_node_par_itr++;
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
                            //Reverse par_nuc and mut_nuc as it will again get flipped in uniq_curr_node_mut
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
                        //Otherwise, if mut_nuc is same it would have been detcted above in parent_parsimony check
                        //At this step the mut_nuc is getting RE-INTRODUCED -> add it as a uniq mutation
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
        //Ensure NOT accessing empty vector (range tree root node empty possibility) 
        if (!node_mappings.find(curr_node)->second.empty()) {
            bool is_root_node = false, is_leaf_node = false;
            auto condensed_front_node = node_mappings.find(curr_node)->second.front();
            if (!((condensed_node_mappings.find(condensed_front_node)->second.empty()) && (node_mappings.find(curr_node)->second.size() == 1))) {
                //REMOVE condensed_front_node if it does not map to original tree
                if ((condensed_node_mappings.find(condensed_front_node)->second.empty())) {
                    node_mappings.find(curr_node)->second.erase(node_mappings.find(curr_node)->second.begin());
                    condensed_front_node = node_mappings.find(curr_node)->second.front();
                }
                //ROOT node check
                else if (condensed_front_node->is_root())
                    is_root_node = true;
                //LEAF node check
                else if (condensed_node_mappings.find(condensed_front_node)->second.front()->is_leaf())
                    is_leaf_node = true;

                //Place as a sibling if common_node_mut is not empty and NOT root in original tree
                if ((common_node_mut.size()) && (!is_root_node)) {
                    //ADD curr_node to common_mut_nodes if NOT leaf in original tree and uniq_curr_node_mut > 0
                    if ((!is_leaf_node) && (uniq_curr_node_mut.size()))
                        common_mut_nodes.insert(i);
                    
                    //Checking min_parsimony
                    int new_min_par = -1; 
                    // If best_par_score is empty and curr_par_score >= limit -> CHANGE
                    if (min_par.par_list.empty())
                        new_min_par = 1;
                    // If curr_par_score < best_par_score and curr_par_score >= limit -> CHANGE
                    else if ((curr_node_par_mut.size() < min_par.par_list[0].size()))
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
                //Place as a child if current node is NOT leaf node in original_tree or a ROOT node
                else if (!is_leaf_node) {
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
                    if (min_par.par_list.empty())
                        new_min_par = 1;
                    // If curr_par_score < best_par_score and curr_par_score >= limit -> CHANGE
                    else if (curr_node_par_mut.size() < min_par.par_list[0].size()) 
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
        curr_par.curr_node = curr_node;
        parsimony_stack.push(curr_par);

        uniq_curr_node_mut.clear();
        common_node_mut.clear();
        curr_node_par_mut.clear();
    }

    //Clear variables
    while (!parsimony_stack.empty())
        parsimony_stack.pop();
    min_par.par_list.clear();

    //Compute node_score_map
    if (check_nodes.empty()) {
        //Getting num_nodes in parallel to calculate score
        size_t num_nodes = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, min_par.idx_list.size()), 0.0, 
            [&](const tbb::blocked_range<size_t> &k, size_t block_num_nodes) -> size_t {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto n_idx = min_par.idx_list[i];
                    //If common_mut_nodes has the given node then it is internal node placed as sibling with other node mutations as well
                    if (common_mut_nodes.find(n_idx) == common_mut_nodes.end())
                        block_num_nodes += node_mappings.find(dfs[n_idx])->second.size();
                    else
                        block_num_nodes++;   
                }
                return block_num_nodes;
            },
            [&] (size_t x, size_t y) -> size_t {
                return x + y;
            }
        );
        
        double score = 1.0 / log2(num_nodes+1);
        
        //Update node_score_map in parallel
        static tbb::affinity_partitioner ap;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, min_par.idx_list.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto n_idx = min_par.idx_list[i];
                    //If common_mut_nodes has the given node then it is internal node placed as sibling with other node mutations as well
                    if (common_mut_nodes.find(n_idx) == common_mut_nodes.end()) {
                        for (auto const& node: node_mappings.find(dfs[n_idx])->second) {
                            tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
                            auto created = node_score_map.insert(ac, std::make_pair(node, score));
                            if (!created)
                                ac->second += score;
                            ac.release();
                        }
                    }
                    else {
                        auto node = node_mappings.find(dfs[n_idx])->second.front();
                        tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
                        auto created = node_score_map.insert(ac, std::make_pair(node, score));
                        if (!created)
                            ac->second += score;
                        ac.release();
                    }
                }
            },
        ap);
        
        min_par.idx_list.clear();
        return 0;
    }
    //Check if given node is present in parsimonious list of read
    else {
        using rwmutex_t = tbb::queuing_rw_mutex;
        rwmutex_t my_mutex_rw;
        static tbb::affinity_partitioner ap;
        bool found = false;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, min_par.idx_list.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    {
                        rwmutex_t::scoped_lock my_lock{my_mutex_rw, false};
                        if (found)
                            continue;
                    }
                    auto n_idx = min_par.idx_list[i];
                    //If common_mut_nodes has the given node then it is internal node placed as sibling with other node mutations as well
                    if (common_mut_nodes.find(n_idx) == common_mut_nodes.end()) {
                        for (auto const& node: node_mappings.find(dfs[n_idx])->second) {
                            if (std::find(check_nodes.begin(), check_nodes.end(), node) != check_nodes.end()) {
                                rwmutex_t::scoped_lock my_lock{my_mutex_rw, true};
                                found = true;
                                break;
                            }
                        }
                    }
                    else {
                        auto node = node_mappings.find(dfs[n_idx])->second.front();
                        if (std::find(check_nodes.begin(), check_nodes.end(), node) != check_nodes.end()) {
                            rwmutex_t::scoped_lock my_lock{my_mutex_rw, true};
                            found = true;
                            break;
                        }
                    }
                }
            },
        ap);

        min_par.idx_list.clear();
        if (found)
            return 1;
        else
            return 0;
    }
}

//Get the clade name
std::string getLineage(const MAT::Tree &T, MAT::Node* n) {
    //Checking all ancestors of a node to get clade
    for (auto anc: T.rsearch(n->identifier, true)) {  
        //Don't consider the clades that start with numbers {Different clade naming}
        auto clade = anc->clade_annotations[1];
        if(clade != "")
            return clade;
    }
    return "";
}

//Function to calculation distance between two nodes
int mutationDistance(const MAT::Tree &T1, const MAT::Tree &T2, const MAT::Node* N1, const MAT::Node* N2) {
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
    tbb::parallel_sort(node1_mutations.begin(), node1_mutations.end(), compareMutations);
    tbb::parallel_sort(node2_mutations.begin(), node2_mutations.end(), compareMutations);
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

//Function to calculation distance between two nodes
int mutationDistance(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations) {
    tbb::parallel_sort(node1_mutations.begin(), node1_mutations.end(), compareMutations);
    tbb::parallel_sort(node2_mutations.begin(), node2_mutations.end(), compareMutations);
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

//Getting all mutations of a haplotype
std::vector<MAT::Mutation> getMutations(const MAT::Tree& T, const std::string sample) {
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

//SORT node_scores
void sortNodeScore(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, const tbb::concurrent_hash_map<MAT::Node*, double>& node_score_map, tbb::concurrent_vector<std::pair<MAT::Node*, double>>& node_score_vector) {
    // CONVERT the concurrent hash map into a vector in parallel
    // Create a range for parallel processing
    tbb::blocked_range<size_t> range(0, node_score_map.size());
     //Use TBB parallel_for to iterate over the concurrent hash map and convert it to a vector
    tbb::parallel_for(range, [&node_score_map, &node_score_vector](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            auto it = node_score_map.begin();
            // Advance the iterator to the correct position
            std::advance(it, i); 
            // Extract the node and score from the concurrent hash map
            auto node = it->first;
            auto score = it->second;
            // Create a pair and add it to the concurrent vector
            node_score_vector.push_back(std::make_pair(node, score));
        }
    });
    node_score_vector.shrink_to_fit();
    //Sort the node_score_vector filled above
    tbb::parallel_sort(node_score_vector.begin(), node_score_vector.end(), 
         [&condensed_node_mappings](const auto& a, const auto& b) {
            return compareNodeScore(condensed_node_mappings, a, b);
    });
}

//Comparing different node_scores  
bool compareNodeScore (const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, const std::pair<MAT::Node*, double>& a, const std::pair<MAT::Node*, double>& b) {
    //Compare based on double values
    if (abs(a.second - b.second) > 1e-9)
        return a.second > b.second;
    else {
        size_t a_leaves = getNumLeaves(condensed_node_mappings, a.first);
        size_t b_leaves = getNumLeaves(condensed_node_mappings, b.first);
        //Compare based on number of leaf nodes
        if (a_leaves != b_leaves)
          return a_leaves > b_leaves;
        //Compare alphabetically
        else
          return a.first->identifier > b.first->identifier;
   }
};

//Get number of leaves
size_t getNumLeaves(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, MAT::Node* condensed_node) {
    size_t leaves = 0;
    if (!condensed_node_mappings.find(condensed_node)->second.empty()) {
        auto front_node = condensed_node_mappings.find(condensed_node)->second.front();
        std::queue<MAT::Node*> remaining_nodes;
        remaining_nodes.push(front_node);
        while (remaining_nodes.size() > 0) {
            MAT::Node* curr_node = remaining_nodes.front();
            remaining_nodes.pop();
            if (curr_node->children.size() == 0)
                leaves++;
            else {
                for (auto c: curr_node->children)
                    remaining_nodes.push(c);
            }
        }
    }
    return leaves;
}

//Comparing mutations for sorting a mutation vector 
bool compareMutations(const MAT::Mutation &a, const MAT::Mutation &b) {
    if (a.position != b.position)
        return a.position < b.position;
    else
        return a.mut_nuc < b.mut_nuc;
}

//Get MAT within range
void createRangeTree(MAT::Node* ref_root, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, const int &start, const int &end, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings, MAT::Tree &T) {
    //Have a queue of current_node (from ref_Tree)and parent_node (from new_Tree) pair
    std::queue<std::pair<MAT::Node*, MAT::Node*>> remaining_nodes;
    auto new_node = T.create_node("DUMMY", -1.0, ref_root->clade_annotations.size());
    node_mappings[new_node] = std::vector<MAT::Node*>();
    remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(ref_root, new_node));
    
    //Add the new_node to the node_mappings
    while(remaining_nodes.size() > 0) {
        auto r_curr_node = remaining_nodes.front().first;
        auto n_parent_node = remaining_nodes.front().second;
        remaining_nodes.pop();

        //Get mutations within specified range
        std::vector<MAT::Mutation> range_mutations;
        for (const auto& mut: r_curr_node->mutations) {
            if ((mut.position >= start) && (mut.position <= end))
                range_mutations.emplace_back(mut);
        }
        
        //Add a new node to tree if mutation found in range
        if (!range_mutations.empty()) {
            //Create a new_node
            auto new_node = T.create_node(r_curr_node->identifier, n_parent_node, -1.0);
            //Add mutations to new_node
            new_node->mutations.reserve(range_mutations.size());
            new_node->mutations.insert(new_node->mutations.end(), range_mutations.begin(), range_mutations.end());
            range_mutations.clear();
            //Add new_node to the node_mappings
            node_mappings[new_node] = {r_curr_node};
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, new_node));
        }
        //ONLY check children if r_curr_node doesn't map to original tree
        else if (condensed_node_mappings.find(r_curr_node)->second.empty()) {
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, n_parent_node));
        }
        //Without any mutation, can only be Placed as a child if r_curr_node is NOT leaf node
        else if (!condensed_node_mappings.find(r_curr_node)->second.front()->is_leaf()) {
            //Add current_node to the n_parent_node's list in node_mappings
            node_mappings[n_parent_node].emplace_back(r_curr_node);
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, n_parent_node));
        }
    }
}

//Get MAT having ONLY desired lineages
void createLineageTree(MAT::Node* ref_root, const std::vector<std::string> &lineage_list, MAT::Tree &T) {
    //Have a queue of current_node (from ref_Tree)and parent_node (from new_Tree) pair
    std::queue<struct node_pair_clade> remaining_nodes;
    auto new_node = T.create_node("DUMMY-LINEAGE", -1.0, ref_root->clade_annotations.size());
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
            new_node->mutations.reserve(r_curr_node->mutations.size());
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

//Condense tree nodes based on read coverage
void createCondensedTree(MAT::Node* ref_root, const std::unordered_map<size_t, struct read_info*> &read_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings, MAT::Tree &T) {
    //MAP of site_read
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map.find(i)->second;
        for (int j = rp->start; j <= rp->end; j++)
            site_read_map.insert(j);
    }

    //REMOVE sites not covered by reads
    std::queue<std::pair<MAT::Node*, MAT::Node*>> remaining_nodes;
    auto new_node = T.create_node("DUMMY-CONDENSED", -1.0, ref_root->clade_annotations.size());
    node_mappings[new_node] = std::vector<MAT::Node*>();
    remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(ref_root, new_node));
    
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
            //Add lineage name
            new_node->clade_annotations[1] = r_curr_node->clade_annotations[1];
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
        else if (!r_curr_node->is_leaf()) {
            //Add current_node to the n_parent_node's list in node_mappings
            node_mappings[n_parent_node].emplace_back(r_curr_node);
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, n_parent_node));
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
std::vector<MAT::Node*> updateNeighborNodes(const MAT::Tree &T, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, const std::vector<MAT::Node*> &curr_peak_nodes, const std::vector<MAT::Node*> &peak_nodes, const tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, std::vector<MAT::Node*> &neighbor_nodes, const int& neighbor_dist_thresh, const int& neighbor_peaks_thresh) {
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
                     [&condensed_node_mappings](const auto& a, const auto& b) {
                        return compareNodeScore(condensed_node_mappings, a, b);
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

//Add neighboring nodes to peaks
void addNeighborNodes(const MAT::Tree &T, std::vector<MAT::Node*> &peak_nodes, const int &neighbor_dist_thresh, const int &neighbor_peaks_thresh) {
    timer.Start();
    std::vector<MAT::Node*> neighbor_nodes;
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, peak_nodes.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                MAT::Node* peak = peak_nodes[i];
                int neighbor_added = 0, branch_length_limit = 0, feasible_branch_length = 0;
                std::vector<MAT::Node*> prev_visited_nodes;
                
                //Keep looking for neighbors until neighbors added == neighbor_peaks_thresh AND neighbors within branch_length_limit 
                while ((neighbor_added < neighbor_peaks_thresh) && (feasible_branch_length == branch_length_limit)) {
                    // 1. Search children nodes within branch_length_limit
                    branch_length_limit++;
                    auto curr_node = peak;
                    int curr_branch_length = 0;
                    feasible_branch_length = childNodesAddition(T, my_mutex, peak_nodes, neighbor_nodes, neighbor_added, prev_visited_nodes, peak, curr_node, curr_branch_length, branch_length_limit, neighbor_peaks_thresh, neighbor_dist_thresh);
                    
                    // 2. Search parent nodes within branch_length_limit
                    while ((curr_branch_length < branch_length_limit) && (neighbor_added < neighbor_peaks_thresh)) {
                        if (curr_node->parent == NULL)
                            break;
                        curr_node = curr_node->parent;
                        curr_branch_length++;
                        int temp_feasible_branch_length = childNodesAddition(T, my_mutex, peak_nodes, neighbor_nodes, neighbor_added, prev_visited_nodes, peak, curr_node, curr_branch_length, branch_length_limit, neighbor_peaks_thresh, neighbor_dist_thresh);
                        if (temp_feasible_branch_length > feasible_branch_length)
                            feasible_branch_length = temp_feasible_branch_length;
                    }
                }
                prev_visited_nodes.clear();
            }
        },
    ap);
    
    //Reserve more memory to ensure there are no segmentation faults
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();
    printf("\nPost neighbor adding PEAK nodes: %d\n", (int)peak_nodes.size());
    fprintf(stderr, "\nNeighbor addition done in %ld sec\n", (timer.Stop()/1000));
}

//Function to search for child nodes within given branch length
int childNodesAddition(const MAT::Tree &T, my_mutex_t& my_mutex, const std::vector<MAT::Node*>& peak_nodes, std::vector<MAT::Node*>& neighbor_nodes, int& neighbor_added, std::vector<MAT::Node*>& prev_visited_nodes, MAT::Node* peak, MAT::Node* init_node, const int& init_branch_length, const int& branch_length_thresh, const int& neighbor_peaks_thresh, const int& neighbor_dist_thresh) {
    std::vector<MAT::Node*> selected_neighbors;
    std::queue<struct node_branch> remaining_nodes_branch;
    struct node_branch init_node_branch = {init_node, init_branch_length};
    remaining_nodes_branch.push(init_node_branch);
    int last_branch_length = -1;
    
    //Conditions to check: neighbor_peaks_thresh, neighbor_dist_thresh, {prev_visited_nodes, neighbor_nodes}, branch_length_thresh 
    while(remaining_nodes_branch.size() > 0) {
        struct node_branch curr_node_branch = remaining_nodes_branch.front();
        remaining_nodes_branch.pop();
        //Only analyze if neighbor_peaks_thresh NOT reached
        if (neighbor_added < neighbor_peaks_thresh) {
            int m_dist = mutationDistance(T, T, curr_node_branch.node, peak);
            //Only analyze nodes and their children within neighbor_dist_thresh
            if (m_dist <= neighbor_dist_thresh) {
                last_branch_length = curr_node_branch.branch_length;
                //Do further checks on this node if NOT analyzed before
                if (std::find(prev_visited_nodes.begin(), prev_visited_nodes.end(), curr_node_branch.node) == prev_visited_nodes.end()) {
                    prev_visited_nodes.emplace_back(curr_node_branch.node);
                    bool found = false;
                    //Check if a peak with same mutations is NOT already present
                    //Not going for exact node match as we want single representative of mutation vector
                    for (const auto& hap: peak_nodes) {
                        int mut_dist = mutationDistance(T, T, hap, curr_node_branch.node);
                        if (!mut_dist) {
                            found = true;
                            break;
                        }
                    }
                    //Check if neighbor peak with same mutations is already accounted
                    if (!found) {
                        my_mutex_t::scoped_lock my_lock{my_mutex};
                        for (const auto& hap: neighbor_nodes) {
                            int mut_dist = mutationDistance(T, T, hap, curr_node_branch.node);
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
                        neighbor_nodes.emplace_back(curr_node_branch.node);
                        my_lock.release();
                        neighbor_added++;
                    }
                }
                //Only add current node's children in list if within branch_length_thresh
                auto child_branch_length = curr_node_branch.branch_length + 1;
                if (child_branch_length <= branch_length_thresh) {
                    for (auto child: curr_node_branch.node->children) {
                        struct node_branch child_node_branch = {child, child_branch_length};
                        remaining_nodes_branch.push(child_node_branch);
                    }
                }
            }
        }

        //Stop checking for neighbors if neighbor_added == neighbor_peak_thresh
        else {
            while (!remaining_nodes_branch.empty())
                remaining_nodes_branch.pop();
        }
    }
    return last_branch_length;
}

//Generate regression based estimate algorithm data
void generateFilteringData(const MAT::Tree &T_orig, const MAT::Tree &T, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, const std::string &ref_seq, const std::vector<MAT::Node*> &peak_nodes, const std::unordered_map<size_t, struct read_info*> &read_map, const std::string &barcode_file, const std::string &read_mutation_depth_vcf, const std::string &condensed_nodes_csv) {
    timer.Start();
    //MAP of peak_node mutations 
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

                //Sort mutations
                tbb::parallel_sort(mut_list.begin(), mut_list.end(), compareMutations);
                
		        //Insert peak and its mutations in peak_mut_map
                tbb::concurrent_hash_map<MAT::Node*, std::vector<MAT::Mutation>>::accessor ac;
                peak_mut_map.insert(ac, std::make_pair(pk, mut_list));
                ac.release();
                mut_list.clear();
            }
        },
    ap);
    

    //MAP of site_reads
    std::unordered_map<int, std::vector<std::pair<char, size_t>>> site_read_map;
    for (size_t i = 0; i < read_map.size(); ++i) {
        auto rp = read_map.find(i)->second;
        auto nuc_pos = rp->start;
        for (const auto& mut: rp->mutations) {
            //Store the coverage of nucleotides leading up to the mutating nucleotide
            while (nuc_pos < mut.position) {
                char curr_nuc = ref_seq[nuc_pos-1];
                if (site_read_map.find(nuc_pos) == site_read_map.end()) 
                    site_read_map.insert(std::make_pair(nuc_pos, std::vector<std::pair<char, size_t>>(1, {curr_nuc, 1})));
                else {
                    bool found = false;
                    for (size_t vec_idx = 0; vec_idx < site_read_map[nuc_pos].size(); vec_idx++) {
                        if (site_read_map[nuc_pos][vec_idx].first == curr_nuc) {
                            site_read_map[nuc_pos][vec_idx].second++;
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        site_read_map[nuc_pos].emplace_back(std::make_pair(curr_nuc, 1));
                }
                nuc_pos++;
            }
            //Store the coverage of mutating nucleotide
            char curr_nuc = MAT::get_nuc(mut.mut_nuc);
            if (site_read_map.find(mut.position) == site_read_map.end()) 
                site_read_map.insert(std::make_pair(mut.position, std::vector<std::pair<char, size_t>>(1, {curr_nuc, 1})));
            else {
                bool found = false;
                for (size_t vec_idx = 0; vec_idx < site_read_map[mut.position].size(); vec_idx++) {
                    if (site_read_map[mut.position][vec_idx].first == curr_nuc) {
                        site_read_map[mut.position][vec_idx].second++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    site_read_map[mut.position].emplace_back(std::make_pair(curr_nuc, 1));
            }
            nuc_pos++;
        }
        //Store the coverage of nucleotide remaining after mutating ones
        while (nuc_pos <= rp->end) {
            char curr_nuc = ref_seq[nuc_pos-1]; 
            if (site_read_map.find(nuc_pos) == site_read_map.end()) 
                site_read_map.insert(std::make_pair(nuc_pos, std::vector<std::pair<char, size_t>>(1, {curr_nuc, 1})));
            else {
                bool found = false;
                for (size_t vec_idx = 0; vec_idx < site_read_map[nuc_pos].size(); vec_idx++) {
                    if (site_read_map[nuc_pos][vec_idx].first == curr_nuc) {
                        site_read_map[nuc_pos][vec_idx].second++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    site_read_map[nuc_pos].emplace_back(std::make_pair(curr_nuc, 1));
            }
            nuc_pos++;
        }
    }

    //ADD mut_nuc in site_read_map from peak_nodes
    for (const auto& pk: peak_nodes) {
        tbb::concurrent_hash_map<MAT::Node*, std::vector<MAT::Mutation>>::const_accessor k_ac;
        if (peak_mut_map.find(k_ac, pk)) {
            //Iterate through mutations of each peak
            auto mut_itr = k_ac->second.begin();
            while (mut_itr != k_ac->second.end()) {
                auto sr_itr = site_read_map.find(mut_itr->position);
                bool found = false;
                for (size_t vec_idx = 0; vec_idx < sr_itr->second.size(); vec_idx++) {
                    if (sr_itr->second[vec_idx].first == MAT::get_nuc(mut_itr->mut_nuc)) {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    sr_itr->second.emplace_back(std::make_pair(MAT::get_nuc(mut_itr->mut_nuc), 0));
                mut_itr++;
            }
        } 
    }

    //Initialize the buffers for writing files
    std::ofstream outfile_barcode(barcode_file, std::ios::out | std::ios::binary);
    std::ofstream outfile_vcf(read_mutation_depth_vcf, std::ios::out | std::ios::binary);
    std::ofstream outfile_condensed(condensed_nodes_csv, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_barcode, outbuf_vcf, outbuf_condensed_nodes;
    if (barcode_file.find(".gz\0") != std::string::npos) {
            outbuf_barcode.push(boost::iostreams::gzip_compressor());
    }
    if (read_mutation_depth_vcf.find(".gz\0") != std::string::npos) {
            outbuf_vcf.push(boost::iostreams::gzip_compressor());
    }
    if (condensed_nodes_csv.find(".gz\0") != std::string::npos) {
            outbuf_condensed_nodes.push(boost::iostreams::gzip_compressor());
    }
    outbuf_barcode.push(outfile_barcode);
    outbuf_vcf.push(outfile_vcf);
    outbuf_condensed_nodes.push(outfile_condensed);
    std::ostream barcode(&outbuf_barcode);
    std::ostream vcf(&outbuf_vcf);
    std::ostream condensed(&outbuf_condensed_nodes);
    
    //WRITING the VCF and header mutations in the barcode file
    std::vector<std::pair<int, char>> read_mutations_list; 
    std::string barcode_print, vcf_print, condensed_print;
    vcf << "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tDEPTH\n";
    for (int i = 1; i < (int)ref_seq.size(); i++) {
        if (site_read_map.find(i) != site_read_map.end()) {
            //Calculate AF
            int total_reads = 0;
            //Iterate through all observed alleles at given position
            for (const auto& mut_count: site_read_map[i])
                total_reads += mut_count.second;
            for (const auto& mut_count: site_read_map[i]) {
                //Only consider mutating nucleotides
                auto ref_nuc = ref_seq[i-1];
                if (mut_count.first != ref_nuc) {
                    //Write the mutations header in barcode file
                    barcode_print += ",";
                    barcode_print += ref_nuc + std::to_string(i) + mut_count.first;
                    read_mutations_list.emplace_back(std::make_pair(i, mut_count.first));
                    float af = (float)mut_count.second / (float)total_reads;
                    vcf_print += "NC_045512v2\t" + std::to_string(i) + "\t" + ref_nuc + std::to_string(i) + mut_count.first + "\t" + ref_nuc + "\t" + mut_count.first + "\t.\t.\tAF=";
                    vcf_print += std::to_string(af) + "\t" + std::to_string(total_reads) + "\n";
                    vcf << vcf_print;
                    vcf_print.clear();
                }
            }
        }
    } 
    barcode << barcode_print;
    barcode_print.clear();

    //WRITING peak mutations in barcode
    size_t uniq_condensedPeak_count = 0;
    std::vector<size_t> mut_idx_list;
    for (const auto &pk_mut: peak_mut_map) {
        //Write peak name in barcode file
        auto uncondensed_peaks = condensed_node_mappings.find(pk_mut.first)->second;

        if (!uncondensed_peaks.empty()) {
            //DO NOT include DUMMY lineage node
            std::string check_string = "DUMMY";
            if (uncondensed_peaks.front()->identifier.find(check_string) != std::string::npos)
                uncondensed_peaks.erase(uncondensed_peaks.begin());
        }

        if (uncondensed_peaks.empty())
            continue;

        barcode_print += "\n";
        if (uncondensed_peaks.size() > 1) {
            std::string peak_name = "CONDENSED-" + std::to_string(++uniq_condensedPeak_count);
            barcode_print += peak_name;
            condensed_print += peak_name;
            for (const auto& node: uncondensed_peaks)
                condensed_print += "," + node->identifier + "_" + getLineage(T_orig, node);
            condensed_print += "\n";
        }
        else 
            barcode_print += uncondensed_peaks.front()->identifier + "_" + getLineage(T_orig, uncondensed_peaks.front());    

        //Iterating through mutations of pk_mut to store indices of its mutations from read_mutations_list
        for (const auto& mut: pk_mut.second) {
            //Search for mutation in read_mutations_list
            auto rm_itr = std::find(read_mutations_list.begin(), read_mutations_list.end(), std::make_pair(mut.position, MAT::get_nuc(mut.mut_nuc)));
            if (rm_itr != read_mutations_list.end())
                mut_idx_list.emplace_back(rm_itr - read_mutations_list.begin());
            else {
                if (site_read_map.find(mut.position) != site_read_map.end())
                    fprintf(stderr, "%d%c: NOT found in read_mutations_list \n", mut.position, MAT::get_nuc(mut.mut_nuc));
                else
                    fprintf(stderr, "%d%c: NOT found in site_read_map\n", mut.position, MAT::get_nuc(mut.mut_nuc));
            }
        }

        //WRITE peak mutations in barcode file
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
        while (idx < read_mutations_list.size()) {
            barcode_print += ",0";
            idx++;
        }

        barcode << barcode_print;
        condensed << condensed_print;
        barcode_print.clear();
        condensed_print.clear();
        mut_idx_list.clear();
    }
    
    peak_mut_map.clear();
    read_mutations_list.clear();    
    fprintf(stderr,"Barcode and VCF file writing took %ld sec\n\n", (timer.Stop() / 1000));
}

//Map reads to haplotypes
void placeReads(const MAT::Tree &T, const std::string &ref_seq, const std::unordered_map<size_t, struct read_info*>& read_map, const std::unordered_map<size_t, struct read_info*>& hap_map) {
    //1. STORE mutations of reads with parsimony > 0 in mut_AFDepth_map
    tbb::concurrent_hash_map<std::pair<int, char>, size_t> mut_AFDepth_map;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_map.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                auto rp = read_map.find(i)->second;
                //Find the most parsimonious positions
                std::vector<size_t> epp_list;
                int min_dist = std::numeric_limits<int>::max();
                std::vector<MAT::Mutation> mutations_in_range;
                auto hap_itr = hap_map.begin();
                while (hap_itr != hap_map.end()) {
                    for (auto m: hap_itr->second->mutations) {
                        if ((m.position >= rp->start) && (m.position <= rp->end))
                            mutations_in_range.emplace_back(m);
                    }
                    int curr_dist = mutationDistance(rp->mutations, mutations_in_range);
                    mutations_in_range.clear();
                    if (curr_dist <= min_dist) {
                        if (curr_dist < min_dist) {
                            epp_list.clear();
                            min_dist = curr_dist;
                        }
                        epp_list.emplace_back(hap_itr->first);
                    }
                    hap_itr++;
                }

                //STORE read mutations when parsimony > 0
                if (min_dist) {
                    tbb::concurrent_hash_map<std::pair<int, char>, size_t>::accessor ac;
                    for (const auto& mut: rp->mutations) {
                        auto created = mut_AFDepth_map.insert(ac, std::make_pair(std::make_pair(mut.position, MAT::get_nuc(mut.mut_nuc)), 1));
                        if (!created) 
                            ac->second++;
                        ac.release();
                    }
                    epp_list.clear();
                }
            }
        },
    ap);

    //2. CREATE site_reads_map
    std::unordered_map<int, std::vector<std::pair<char, size_t>>> site_reads_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map.find(i)->second;
        auto nuc_pos = rp->start;
        for (const auto& mut: rp->mutations) {
            //Store the coverage of nucleotides leading up to the mutating nucleotide
            while (nuc_pos < mut.position) {
                char curr_nuc = ref_seq[nuc_pos-1];
                if (site_reads_map.find(nuc_pos) == site_reads_map.end()) 
                    site_reads_map.insert(std::make_pair(nuc_pos, std::vector<std::pair<char, size_t>>(1, {curr_nuc, 1})));
                else {
                    bool found = false;
                    for (size_t vec_idx = 0; vec_idx < site_reads_map[nuc_pos].size(); vec_idx++) {
                        if (site_reads_map[nuc_pos][vec_idx].first == curr_nuc) {
                            site_reads_map[nuc_pos][vec_idx].second++;
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        site_reads_map[nuc_pos].emplace_back(std::make_pair(curr_nuc, 1));
                }
                nuc_pos++;
            }
            //Store the coverage of mutating nucleotide
            char curr_nuc = MAT::get_nuc(mut.mut_nuc);
            if (site_reads_map.find(mut.position) == site_reads_map.end()) 
                site_reads_map.insert(std::make_pair(mut.position, std::vector<std::pair<char, size_t>>(1, {curr_nuc, 1})));
            else {
                bool found = false;
                for (size_t vec_idx = 0; vec_idx < site_reads_map[mut.position].size(); vec_idx++) {
                    if (site_reads_map[mut.position][vec_idx].first == curr_nuc) {
                        site_reads_map[mut.position][vec_idx].second++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    site_reads_map[mut.position].emplace_back(std::make_pair(curr_nuc, 1));
            }
            nuc_pos++;
        }
        //Store the coverage of nucleotide remaining after mutating ones
        while (nuc_pos <= rp->end) {
            char curr_nuc = ref_seq[nuc_pos-1]; 
            if (site_reads_map.find(nuc_pos) == site_reads_map.end()) 
                site_reads_map.insert(std::make_pair(nuc_pos, std::vector<std::pair<char, size_t>>(1, {curr_nuc, 1})));
            else {
                bool found = false;
                for (size_t vec_idx = 0; vec_idx < site_reads_map[nuc_pos].size(); vec_idx++) {
                    if (site_reads_map[nuc_pos][vec_idx].first == curr_nuc) {
                        site_reads_map[nuc_pos][vec_idx].second++;
                        found = true;
                        break;
                    }
                }
                if (!found)
                    site_reads_map[nuc_pos].emplace_back(std::make_pair(curr_nuc, 1));
            }
            nuc_pos++;
        }
    }
    
    //COMPUTE average_read_depth and AF and read_depth of mutations in mut_AFDepth_map
    double average_read_depth = 0.0;
    for (const auto& sr_list: site_reads_map) {
        for (const auto& mut_count: sr_list.second)
            average_read_depth += mut_count.second;
    }
    average_read_depth /= site_reads_map.size();

    /////////////////////////////////////////
    //for (double af_thresh = 0.05; af_thresh < 0.5; af_thresh += 0.05) {
    //    for (double rel_depth_thresh = 0.1; rel_depth_thresh <= 1; rel_depth_thresh += 0.1) {
    //        size_t depth_thresh = rel_depth_thresh * average_read_depth;
    //        size_t tp = 0, fp = 0;
    //        //Iterating through mut_AFDepth_map to get mutations satisfying the thresholds
    //        for (const auto& m_ad: mut_AFDepth_map) {
    //            size_t total_reads = 0, mut_reads = 0;
    //            int pos = m_ad.first.first;
    //            char mut_nuc = m_ad.first.second;
    //            for (const auto& mut_count: site_reads_map[pos]) {
    //                if (mut_count.first == mut_nuc)
    //                    mut_reads = mut_count.second;
    //                total_reads += mut_count.second;
    //            }
    //            double curr_af = (double)mut_reads / (double)total_reads;
    //        }
    //        printf("AF_thresh: %f, Depth_thresh_rel_avg_depth: %f, TP: %lu, FP: %lu\n", af_thresh, rel_depth_thresh, tp, fp);
    //    }
    //}

    ////Iterate through hap_read_map to check Allele Frequency at every site of haplotype
    //for (const auto& hr: hap_read_map) {
    //    auto hap_idx = hr.first;
    //    auto reads_idx_list = hr.second;
    //    //Get Allele Frequency of read with non-zero parsimony
    //    tbb::concurrent_hash_map<size_t, std::vector<size_t>>::const_accessor k_ac;
    //    if (hap_nonZeroreads_map.find(k_ac, hap_idx)) {
    //        if (hap_map.find(hap_idx)->second->read == "CONDENSED-28_B.39") {
    //            printf("\nHAP: %s\n", hap_map.find(hap_idx)->second->read.c_str());
    //            for (const auto& mut: hap_map.find(hap_idx)->second->mutations)
    //                printf("%d%c\n", mut.position, MAT::get_nuc(mut.mut_nuc));
    //            printf("\n Read Mutations\n");
    //        }
    //        auto nonZeroread_idx_list = k_ac->second;
    //        //Get max_read_depth
    //        size_t max_read_depth = 0;
    //        for (const auto& sr_list: site_read_map) {
    //            size_t curr_read_depth = 0;
    //            for (const auto& mut_count: sr_list.second)
    //                curr_read_depth += mut_count.second;
    //            if (curr_read_depth > max_read_depth)
    //                max_read_depth = curr_read_depth;
    //        }
    //        //Analyze the mutations of non-zero reads
    //        for (const auto& rd_idx: nonZeroread_idx_list) {
    //            auto rp = read_map.find(rd_idx)->second;
    //            for (const auto& mut: rp->mutations) {
    //                size_t total_reads = 0;
    //                for (const auto& mut_count: site_read_map[mut.position])
    //                    total_reads += mut_count.second;
    //                for (const auto& mut_count: site_read_map[mut.position]) {
    //                    if (mut_count.first == MAT::get_nuc(mut.mut_nuc)) {
    //                        if (hap_map.find(hap_idx)->second->read == "CONDENSED-28_B.39")
    //                        printf("%d%c -> %f, %f\n", mut.position, mut_count.first,  ((double)mut_count.second / (double)total_reads), ((double)mut_count.second) / (double)max_read_depth);
    //                        break;
    //                    }
    //                }
    //            }
    //        }    
    //    }
    //    site_read_map.clear();
    //}
}

//Comparing rd_idx for sorting a rd_idx_vector 
bool compareIdx(const std::pair<int, size_t> &a, const std::pair<int, size_t> &b) {
    return a.second < b.second;
}

//Store the Reads from the SAM in read_map and site_MutReads_map
void readSAM(const std::string &sam_file, const std::string &ref_seq, std::unordered_map<size_t, std::string> &read_map, std::unordered_map<int, std::vector<std::tuple<std::string, std::string, std::vector<size_t>>>> &site_MutReads_map) {
    //NOTE: Only considers I,D,N,M in CIGAR
    boost::filesystem::ifstream fileHandler(sam_file);
    std::string s;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            //Skip header
            if (words[0][0] == '@')
                continue;
            else {
                //Skip reads if UNMAPPED
                if (std::stoi(words[1]) == 4)
                    continue;

                //Capturing CIGAR string
                std::regex cigar_pattern(R"(\d+[A-Za-z])"); // Matches one or more digits followed by a character
                std::sregex_iterator r_iter(words[5].begin(), words[5].end(), cigar_pattern);
                std::sregex_iterator end;
                std::vector<std::pair<int, char>> cigar_chunks;
                int cigar_len = 0;
                while (r_iter != end) {
                    std::string sub_cigar = (*r_iter).str();
                    //Capturing INDELS
                    std::regex pos_char("([0-9]+)([A-Za-z])");
                    std::sregex_iterator regex_iter(sub_cigar.begin(), sub_cigar.end(), pos_char);
                    std::smatch pc_match = *regex_iter;
                    //Convert the captured substrings to int and char and store in cigar_chunks
                    cigar_chunks.emplace_back(std::make_pair(std::stoi(pc_match.str(1)), pc_match.str(2)[0]));
                    //INSERTIONS don't increase length of sequence
                    if (pc_match.str(2) != "I")
                        cigar_len += std::stoi(pc_match.str(1));
                    r_iter++;
                } 

                //Find start and end idx in read_seq
                int start_idx = std::stoi(words[3]);
                std::string read_seq = words[9];
                int end_idx = start_idx + cigar_len - 1;

                //Add read to read_map
                auto rd_idx = read_map.size();
                read_map.insert({rd_idx, (words[0] + "_READ_" + std::to_string(start_idx) + "_" + std::to_string(end_idx))});
                
                //Update site_MutReads_map
                int curr_sub_len, seq_idx = 0, ins_count = 0, del_count = 0;
                std::string alt_nuc, ref_nuc;
                for (size_t i = 0; i < cigar_chunks.size(); i++) {
                    auto [cig_len, cig_val] = cigar_chunks[i];
                    curr_sub_len = 0;
                    while (curr_sub_len < cig_len) {
                        int nuc_pos = start_idx + seq_idx - ins_count + del_count;
                        switch (cig_val) {
                            case 'I':
                                //Insertion at pos 1 requires ref_nuc at pos 1 after the alt_nuc
                                if (nuc_pos == 1) {
                                    ref_nuc = ref_seq[0];
                                    alt_nuc = read_seq.substr(seq_idx, cig_len) + ref_nuc;
                                }
                                //Otherwise the ref_nuc needs to be from prev pos 
                                else {
                                    alt_nuc = ref_nuc + read_seq.substr(seq_idx, cig_len);
                                    nuc_pos -= 1; 
                                }
                                curr_sub_len += cig_len;   
                                seq_idx += cig_len;
                                ins_count += cig_len;
                                break;

                            case 'D':
                                //Deletion at pos 1 requires alt_nuc at pos 1 after the ref_nuc
                                if (nuc_pos == 1) {
                                    alt_nuc = ref_seq[cig_len];
                                    ref_nuc = ref_seq.substr(0, cig_len) + alt_nuc;
                                }
                                else { 
                                    alt_nuc = ref_nuc;
                                    ref_nuc += ref_seq.substr(nuc_pos - 1, cig_len);
                                    nuc_pos -= 1;
                                }
                                curr_sub_len += cig_len;   
                                del_count += cig_len;
                                break;
                            
                            case 'N': 
                                ref_nuc = ref_seq[nuc_pos - 1]; 
                                alt_nuc = "N";
                                curr_sub_len++;   
                                seq_idx++;
                                break;

                            default:
                                ref_nuc = ref_seq[nuc_pos - 1]; 
                                alt_nuc = std::string(1, read_seq[seq_idx]);
                                curr_sub_len++;   
                                seq_idx++;
                        }
                        
                        auto smr_itr = site_MutReads_map.find(nuc_pos);
                        //Current position NOT in site_MutReads_map
                        if (smr_itr == site_MutReads_map.end()) {
                            std::vector<std::tuple<std::string, std::string, std::vector<size_t>>> newEntry;
                            newEntry.emplace_back(ref_nuc, alt_nuc, std::vector<size_t>(1, rd_idx));
                            site_MutReads_map[nuc_pos] = newEntry;
                        }
                        else {
                            bool found = false;
                            //Add rd_idx to alt_nuc if PRESENT
                            for (auto& mr_tuple: smr_itr->second) {
                                if ((std::get<0>(mr_tuple) == ref_nuc) && (std::get<1>(mr_tuple) == alt_nuc)) {
                                    auto& rd_idx_vector = std::get<2>(mr_tuple); 
                                    rd_idx_vector.emplace_back(rd_idx);
                                    found = true;
                                    break;
                                }
                            }
                            //Else add alt_nuc to smr_itr->second
                            if (!found)
                                smr_itr->second.emplace_back(ref_nuc, alt_nuc, std::vector<size_t>(1, rd_idx));
                        }
                    }
                }
            } 
        }
    }
}
