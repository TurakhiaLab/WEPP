#include "experiment.hpp"

po::variables_map parse_place_read_command(po::parsed_options parsed) {
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
    ("write-vcf,v", po::value<std::string>()->default_value(""),
     "Output VCF file representing selected subtree. Default is full tree")
    ("ref-fasta,f", po::value<std::string>()->default_value(""),
     "Input Fasta file representing reference sequence")
    ("distribution,d", po::value<std::string>()->default_value(""),
     "Give the distribution of samples, comma delimited.")
    ("read-length,r", po::value<int>()->default_value(100),
     "Give the read length of samples")
    ("haplotype-samples,w", po::value<int>()->default_value(10),
     "Give the number of haplotype samples")
    ("read-error,e", po::value<std::string>()->default_value(""),
     "Give the error in sequence reads")
    ("sequence-depth,s", po::value<int>()->default_value(1),
     "Give the sequenceing depth of samples")
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

void simulate_and_place_reads (po::parsed_options parsed) {
    bool old_vcf = true;
    //main argument for the complex extract command
    //uses included code from multiple modules
    //specifically, the modules select, describe, convert, and filter support this command
    po::variables_map vm = parse_place_read_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string ref_mat_filename = vm["ref-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();

    std::string cmd_lineage = vm["lineage"].as<std::string>();
    std::string cmd_distribution = vm["distribution"].as<std::string>();
    std::string cmd_read_error = vm["read-error"].as<std::string>();
    int read_length = vm["read-length"].as<int>();
    int sample_size = vm["haplotype-samples"].as<int>();
    int sequence_depth = vm["sequence-depth"].as<int>();
    
    std::vector<std::string> in_lineage;
    std::stringstream lin_str(cmd_lineage);
    std::string str;
    while (std::getline(lin_str,str,',')) {
        in_lineage.emplace_back(str);
    }

    std::vector<float> in_distribution;
    std::stringstream dist_str(cmd_distribution);
    while (std::getline(dist_str,str,',')) {
        in_distribution.emplace_back(std::stof(str));
    }

    std::vector<float> read_error;
    std::stringstream re_str(cmd_read_error);
    while (std::getline(re_str,str,',')) {
        read_error.emplace_back(std::stof(str));
    }

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    std::string vcf_filename_samples = dir_prefix + vm["write-vcf"].as<std::string>() + "_samples.vcf";
    std::string vcf_filename_reads = dir_prefix + vm["write-vcf"].as<std::string>() + "_reads.vcf";
    std::string vcf_filename_reads_freyja = dir_prefix + vm["write-vcf"].as<std::string>() + "_reads_freyja.vcf";
    std::string depth_filename_reads_freyja = dir_prefix + vm["write-vcf"].as<std::string>() + "_reads_freyja.depth";
    std::string mismatch_matrix_file = dir_prefix + vm["write-vcf"].as<std::string>() + "_mismatch_matrix.csv";
    std::string barcode_file = dir_prefix + vm["write-vcf"].as<std::string>() + "_barcode.csv";
    std::string read_abundance_vcf = dir_prefix + vm["write-vcf"].as<std::string>() + "_abundance.vcf";
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

    fprintf(stderr, "Ref Seq size: %ld\n", ref_seq.size()); 
    //std::cout << read_error[0] << "\n";

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

    fprintf(stderr,"Ref Seq Length: %ld\n", ref_seq.size());
    fprintf(stderr,"Leaves in tree: %ld\n\n", T.get_num_leaves());


    //////////////////////////////////////////////////////// Main Code HERE
    timer.Start();
    
    std::vector<std::vector<Mutation_Annotated_Tree::Node*>> ancestors, all_lineages;
    std::map<int, std::vector<struct ances_sample_list*>*> sample_map;
    std::unordered_map<int, std::vector<Mutation_Annotated_Tree::Node*>*> back_mut_map;
    std::vector<MAT::Node*> dfs, lineage_list, lineage_selected;
    
    //Depth first expansion to get all nodes in the tree and 
    // comparison with given lineage to get all nodes of the required lineage 
    dfs = T.depth_first_expansion(); 
    
    if (!old_vcf) {
        for (auto lineage: in_lineage) {
            for (auto n: dfs) {
                if (n->clade_annotations[1] == lineage) {
                    std::queue<Mutation_Annotated_Tree::Node*> remaining_nodes;
                    remaining_nodes.push(n);
                    while (remaining_nodes.size() > 0) {
                        Mutation_Annotated_Tree::Node* curr_node = remaining_nodes.front();
                        remaining_nodes.pop();
                        if ((curr_node->clade_annotations[1] == "") || (curr_node->clade_annotations[1] == lineage)) {
                            if (curr_node->children.size() == 0)
                                lineage_list.emplace_back(curr_node);
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
            all_lineages.emplace_back(lineage_list);
            lineage_list.clear();
        }

        clock_t time = clock();
        srand(int(time));
        
        //Random selection of required samples from a lineage
        std::vector<std::vector<Mutation_Annotated_Tree::Node*>>::iterator lineage_ptr = all_lineages.begin(); 
        for (auto dist: in_distribution) {
            for (int i = 0; i < ceil(dist * sample_size); i++) {
                int rand_val = int(rand() % lineage_ptr->size());
                lineage_selected.emplace_back((*lineage_ptr)[rand_val]);
            }
            lineage_ptr++;
        }

        for (auto sample: lineage_selected) {
            auto clade = get_clade(T, sample);
            printf("Sample: %s, Clade: %s\n", sample->identifier.c_str(), clade.c_str());
        }

        fprintf(stderr, "\n%ld Samples Selected in %ld msec \n\n", lineage_selected.size(), timer.Stop());

        int num_reads = sequence_depth * ((int(ref_seq.size()) / read_length) + ((int(ref_seq.size()) % read_length) != 0));
        fprintf(stderr, "Num reads per sample: %d\n", num_reads);
    }
   

    std::unordered_map<size_t, struct read_info*> read_map;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score_map;
    std::vector<std::string> vcf_samples;

    //read samples.vcf to check how close are peaks to samples 
    read_sample_vcf(vcf_samples, vcf_filename_samples);
    //Get the reads.vcf data
    read_vcf(read_map, vcf_filename_reads);
    
    //Core algorithm 
    analyze_reads(T, T_ref, read_map, node_score_map, vcf_samples, barcode_file, read_abundance_vcf);
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
            //Checking for header
            size_t field_offset = 9;
            if (words[1] == "POS") {
                header_found = true;
                //Leave certain fields based on our VCF format
                for (size_t j = field_offset; j < words.size(); j++) {
                    struct read_info* rp = new struct read_info;
                    rp->read = words[j];
                    //Get start-end positions and depth of the Read
                    std::regex rgx(".*READ_(\\w+)_(\\w+)_(\\w+)");
                    std::smatch match;
                    if (std::regex_search(words[j], match, rgx)) {
                        rp->start = std::stoi(match[1]);
                        rp->end = std::stoi(match[2]);
                        rp->depth = std::stoi(match[3]);
                    }
                    //Consider the entire genome
                    else {
                        rp->start = 1;
                        rp->end = 29903;
                        rp->depth = 1;
                    }
                    read_map.insert({(j-field_offset), rp});
                }
            }
            else if (header_found) {
                std::vector<std::string> alleles;
                alleles.clear();
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
                            if (mut_present) {
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
                mut_read_map.clear();
            }
        }
    }

    fprintf(stderr,"%s parsed in %ld sec\n\n", vcf_filename_reads.c_str(), (timer.Stop() / 1000));   
}

//Parsimonious placement search for reads
int place_reads(const MAT::Tree &T, const struct read_info* rp, const std::vector<MAT::Node*> &check_nodes, tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings) {
    auto dfs = T.depth_first_expansion();
    std::stack<struct parsimony> parsimony_stack;
    struct min_parsimony min_par;
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
                //Check in Mutation position found in parsimony of parent node
                auto curr_node_par_itr = curr_node_par_mut.begin();
                while (curr_node_par_itr != curr_node_par_mut.end()) {
                    if (curr_node_par_itr->position == node_mut.position) {
                        //mut_nuc matches => remove mutation from parsimony
                        if ((curr_node_par_itr->mut_nuc == node_mut.mut_nuc) || (curr_node_par_itr->mut_nuc == 0b1111)) {
                           curr_node_par_mut.erase(curr_node_par_itr);
                           common_node_mut.emplace_back(*curr_node_par_itr);
                        }
                        //Update the par_nuc in parsimony if only position matches
                        else {
                            curr_node_par_itr->par_nuc = node_mut.mut_nuc;
                        }
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
        //Place as a sibling if common_node_mut is not empty and NOT root
        if ((common_node_mut.size()) && (!curr_node->is_root())) {
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
        else if (!node_mappings.find(curr_node)->second.front()->is_leaf()) {
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
                    block_num_nodes += node_mappings.find(dfs[n_idx])->second.size();
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
                    for (auto const& node: node_mappings.find(dfs[n_idx])->second) {
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
                    for (auto const& node: node_mappings.find(dfs[n_idx])->second) {
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

//Main peak search algorithm
void analyze_reads(const MAT::Tree &T, const MAT::Tree &T_ref, const std::unordered_map<size_t, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, const std::vector<std::string> &vcf_samples, const std::string &barcode_file, const std::string &read_mutation_depth_vcf) {
    timer.Start();
    int top_n = 25, m_dist_thresh = 2, neighbor_dist_thresh = 7, neighbor_peaks_thresh = 100, tree_range = 600, tree_increment = 400; 
    std::vector<MAT::Node*> dfs, peak_nodes, curr_peak_nodes, prohibited_nodes, neighbor_nodes, curr_neighbor_nodes;
    std::vector<size_t> remaining_reads(read_map.size());
    
    //Depth first expansion to get all nodes in the tree and 
    dfs = T.depth_first_expansion();
    
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

    //Iterate till no reads left
    while ((int)remaining_reads.size() > 0) {
        printf("\n");
        fprintf(stderr, "\n");
        
        //Calculating node score for remaining reads
        static tbb::affinity_partitioner ap;
        using my_mutex_t = tbb::queuing_mutex;
        my_mutex_t my_mutex;
        
        //MAPPING reads to nodes
        std::vector<size_t>remaining_read_ids = remaining_reads;
        for (int start = 0; start < 30000; start += tree_increment) {
            int end = start + tree_range;
            std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> node_mappings;
            MAT::Tree range_Tree;
            std::vector<size_t> curr_read_ids;
            
            //Extract all the remaining_read_ids within range
            tbb::parallel_for(tbb::blocked_range<size_t>(0, remaining_read_ids.size()),
                [&](tbb::blocked_range<size_t> k) {
                    for (size_t i = k.begin(); i < k.end(); ++i) {
                        auto rm_idx = remaining_read_ids[i];
                        auto read_id = read_map.find(rm_idx)->second;	
                        if ((read_id->start >= start) && (read_id->end <= end)) { 
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            curr_read_ids.emplace_back(rm_idx);
                        }
                    }
                },
            ap);
            //Sort curr_read_ids
            tbb::parallel_sort(curr_read_ids.begin(), curr_read_ids.end());
            //Remove curr_read_ids from remaining_read_ids
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
                get_range_Tree(T.root, start, end, node_mappings, range_Tree);

            //Tree Search
            tbb::parallel_for(tbb::blocked_range<size_t>(0, curr_read_ids.size()),
                [&](tbb::blocked_range<size_t> k) {
                    for (size_t i = k.begin(); i < k.end(); ++i) {
                        size_t rm_idx = curr_read_ids[i];
                        auto read_id = read_map.find(rm_idx)->second;	
                        place_reads(range_Tree, read_id, curr_peak_nodes, node_score_map, node_mappings);
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
            fprintf(stderr,"Smaller MAT Code NOT working properly!!!");
        
        //REMOVE prohibited_nodes from node_score_map
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

        //SORT node_scores
        tbb::concurrent_vector<std::pair<MAT::Node*, double>> node_score_vector;
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
             [&T](const auto& a, const auto& b) {
                return compare_node_score(T, a, b);
        });

        //GET top_n top_score values from node_score_vector into top_n_node_scores
        std::vector<std::pair<MAT::Node*, double>> top_n_node_scores;
        top_n_node_scores.reserve(std::min(top_n, (int)node_score_vector.size()));
        auto top_score = node_score_vector.begin()->second;
        for (int i = 0; i < std::min(top_n, (int)node_score_vector.size()); ++i) {
            auto n_s = node_score_vector[i];
            //Only consider top_n nodes that have score equal to top node
            if (abs(top_score - n_s.second) < 1e-9)
                top_n_node_scores.emplace_back(n_s);
            else 
                break;
        }
        node_score_vector.clear();

        //REMOVE peak from top_n_node_scores that are in each other's neighborhood
        //Finding top_n_node_scores in neighborhood
        std::vector<MAT::Node*> top_n_node_scores_remove_nodes;
        for (int idx = 0; idx < (int)top_n_node_scores.size(); idx++) {
            auto ref_n_s = top_n_node_scores[idx];
            int num_check_nodes = (int)top_n_node_scores.size() - idx - 1;
            if (std::find(top_n_node_scores_remove_nodes.begin(), top_n_node_scores_remove_nodes.end(), ref_n_s.first) != top_n_node_scores_remove_nodes.end())
                continue;
            tbb::parallel_for(tbb::blocked_range<int>(0, num_check_nodes),
                [&](tbb::blocked_range<int> k) {
                    for (int i = k.begin(); i < k.end(); ++i) {
                        auto curr_n_s = top_n_node_scores[i + idx + 1];
                        //Only consider curr_n_s if mutation_distance > m_dist_thresh
                        int m_dist = mutation_distance(T, T, ref_n_s.first, curr_n_s.first);
                        if (m_dist <= m_dist_thresh) {
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            top_n_node_scores_remove_nodes.emplace_back(curr_n_s.first);
                        }
                    }
                },
            ap);
        }
        //Removing nieghborhood nodes from top_n_node_scores
        auto top_n_itr = top_n_node_scores.begin();
        while (top_n_itr != top_n_node_scores.end()) {
            if (std::find(top_n_node_scores_remove_nodes.begin(), top_n_node_scores_remove_nodes.end(), top_n_itr->first) != top_n_node_scores_remove_nodes.end())
                top_n_itr = top_n_node_scores.erase(top_n_itr);
            else {
                curr_peak_nodes.emplace_back(top_n_itr->first);
                auto curr_clade = get_clade(T, top_n_itr->first);
                printf("PEAK: %s, Score: %f, Clade:%s, reads: %d\n",top_n_itr->first->identifier.c_str(), top_n_itr->second, curr_clade.c_str(), (int)remaining_reads.size());
                fprintf(stderr, "PEAK: %s, Score: %f, Clade:%s, reads: %d\n",top_n_itr->first->identifier.c_str(), top_n_itr->second, curr_clade.c_str(), (int)remaining_reads.size());
                top_n_itr++;
            }
        } 
        top_n_node_scores_remove_nodes.clear();
        top_n_node_scores.clear();
        
        //ADD neighbors of curr_peak_nodes to neighbor_nodes
        curr_neighbor_nodes = update_neighbor_nodes(T, curr_peak_nodes, peak_nodes, node_score_map, neighbor_nodes, neighbor_dist_thresh, neighbor_peaks_thresh);
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
        update_prohibited_nodes(T, curr_peak_nodes, prohibited_nodes, m_dist_thresh);
        
        //REMOVE reads mapped to curr_peak_nodes
        std::vector<size_t> remove_reads;
        remaining_read_ids = remaining_reads;
        for (int start = 0; start < 30000; start += tree_increment) {
            int end = start + tree_range;
            std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> node_mappings;
            MAT::Tree range_Tree;
            std::vector<size_t> curr_read_ids;
            
            //Extract all the remaining_read_ids within range
            tbb::parallel_for(tbb::blocked_range<size_t>(0, remaining_read_ids.size()),
                [&](tbb::blocked_range<size_t> k) {
                    for (size_t i = k.begin(); i < k.end(); ++i) {
                        auto rm_idx = remaining_read_ids[i];
                        auto read_id = read_map.find(rm_idx)->second;	
                        if ((read_id->start >= start) && (read_id->end <= end)) { 
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            curr_read_ids.emplace_back(rm_idx);
                        }
                    }
                },
            ap);
            //Sort curr_read_ids
            tbb::parallel_sort(curr_read_ids.begin(), curr_read_ids.end());
            //Remove curr_read_ids from remaining_read_ids
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
                get_range_Tree(T.root, start, end, node_mappings, range_Tree);

            //Tree Search
            tbb::parallel_for(tbb::blocked_range<size_t>(0, curr_read_ids.size()),
                [&](tbb::blocked_range<size_t> k) {
                    for (size_t i = k.begin(); i < k.end(); ++i) {
                        size_t rm_idx = curr_read_ids[i];
                        auto read_id = read_map.find(rm_idx)->second;	
                        int read_present = place_reads(range_Tree, read_id, curr_peak_nodes, node_score_map, node_mappings);
                        //If Read contains current node -> REMOVE
                        if (read_present) {
                            my_mutex_t::scoped_lock my_lock{my_mutex};
                            remove_reads.emplace_back(rm_idx);
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
            fprintf(stderr,"Smaller MAT Code NOT working properly!!!");

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
            int curr_dist = mutation_distance(T, T_ref, pn, sample_node);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = pn;
            }
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld msec\n\n", timer.Stop());
    
    //REMOVE nodes already present in peak_nodes from neighbor_nodes
    for (const auto& node: peak_nodes) {
        auto nn_itr = neighbor_nodes.begin();
        while (nn_itr != neighbor_nodes.end()) {
            if (*nn_itr == node) {
                neighbor_nodes.erase(nn_itr);
                break;
            }
            nn_itr++;
        }
    }
    //ADD neighbor_nodes to peak_nodes
    peak_nodes.reserve(peak_nodes.size() + neighbor_nodes.size());
    peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    neighbor_nodes.clear();
    
    generate_regression_abundance_data(T, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf);
}

//Add neighboring nodes to prohibited nodes
void update_prohibited_nodes(const MAT::Tree &T, const std::vector<MAT::Node*> &curr_peak_nodes, std::vector<MAT::Node*> &prohibited_nodes, const int& m_dist_thresh) {
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
std::vector<MAT::Node*> update_neighbor_nodes(const MAT::Tree &T, const std::vector<MAT::Node*> &curr_peak_nodes, const std::vector<MAT::Node*> &peak_nodes, const tbb::concurrent_hash_map<MAT::Node*, double> &node_score_map, std::vector<MAT::Node*> &neighbor_nodes, const int& neighbor_dist_thresh, const int& neighbor_peaks_thresh) {
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
                    int m_dist = mutation_distance(T, T, n, peak);
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
                    int m_dist = mutation_distance(T, T, present_node, peak);
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
                                int mut_dist = mutation_distance(T, T, hap, present_node);
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
                                            int mut_dist = mutation_distance(T, T, hap, present_node);
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
                                            int mut_dist = mutation_distance(T, T, hap, present_node);
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
                        return compare_node_score(T, a, b);
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
    tbb::parallel_sort(node1_mutations.begin(), node1_mutations.end(), compare_mutations);
    tbb::parallel_sort(node2_mutations.begin(), node2_mutations.end(), compare_mutations);
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
            //If mutation not found then don't check rest of peaks
            if (cmp_mut_itr == peak_itr->second.end())
                break;
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
    tbb::parallel_sort(peak_mut_list.begin(), peak_mut_list.end(), compare_mutations);
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
        tbb::parallel_sort(mut_idx_list.begin(), mut_idx_list.end());
        
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
        
    fprintf(stderr,"Barcode and VCF file writing took %ld sec\n\n", (timer.Stop() / 1000));
}

//Comparing mutations for sorting a vector 
bool compare_mutations(const MAT::Mutation &a, const MAT::Mutation &b) {
    return a.position < b.position;
}

//Comparing different node_scores  
bool compare_node_score (const MAT::Tree &T, const std::pair<MAT::Node*, double>& a, const std::pair<MAT::Node*, double>& b) {
    //Compare based on double values
    if (abs(a.second - b.second) > 1e-9)
        return a.second > b.second;
    else {
        size_t a_leaves = get_num_leaves(T, a.first);
        size_t b_leaves = get_num_leaves(T, b.first);
        //Compare based on number of leaf nodes
        if (a_leaves != b_leaves)
          return a_leaves > b_leaves;
        //Compare alphabetically
        else
          return a.first->identifier > b.first->identifier;
   }
};

//Get number of leaves
size_t get_num_leaves(const MAT::Tree &T, MAT::Node* node) {
    std::vector<MAT::Node*> leaves;
    std::queue<MAT::Node*> remaining_nodes;
    remaining_nodes.push(node);
    while (remaining_nodes.size() > 0) {
        MAT::Node* curr_node = remaining_nodes.front();
        remaining_nodes.pop();
        if (curr_node->children.size() == 0)
            leaves.emplace_back(curr_node);
        else {
            for (auto c: curr_node->children)
                remaining_nodes.push(c);
        }
    }
    return leaves.size();
}

//Get MAT within range
void get_range_Tree(MAT::Node* ref_root, const int &start, const int &end, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings, MAT::Tree &T) {
    //Have a queue of current_node (from ref_Tree)and parent_node (from new_Tree) pair
    std::queue<std::pair<MAT::Node*, MAT::Node*>> remaining_nodes;
    auto new_node = T.create_node("DUMMY", -1.0, 0);
    node_mappings[new_node] = std::vector<MAT::Node*>();
    remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(ref_root, new_node));
    
    //Add the new_node to the node_mappings
    while(remaining_nodes.size() > 0) {
        auto r_curr_node = remaining_nodes.front().first;
        auto n_parent_node = remaining_nodes.front().second;
        remaining_nodes.pop();
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
            for (const auto& mut: range_mutations)
                new_node->mutations.emplace_back(mut);
            range_mutations.clear();
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