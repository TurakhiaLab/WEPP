#include "experiment.hpp"

po::variables_map parse_place_read_command(po::parsed_options parsed) {
    uint32_t num_cores = tbb::task_scheduler_init::default_num_threads();
    std::string num_threads_message = "Number of threads to use when possible [DEFAULT uses all available cores, " + std::to_string(num_cores) + " detected on this machine]";
    po::variables_map vm;
    po::options_description conv_desc("place_read options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
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
    bool old_vcf = false; // If true, use older VCF

    //main argument for the complex extract command
    //uses included code from multiple modules
    //specifically, the modules select, describe, convert, and filter support this command
    po::variables_map vm = parse_place_read_command(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
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

    fprintf(stderr, "\nLoading input MAT file %s.\n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T;
    if (input_mat_filename.find(".pb\0") != std::string::npos) {
        T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T.uncondense_leaves();
    } else if (input_mat_filename.find(".json\0") != std::string::npos) {
        T = load_mat_from_json(input_mat_filename);
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
    dfs = T.depth_first_expansion(T.root); 
    
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
            //printf("Sample Name: %s\n", (*lineage_ptr)[rand_val]->identifier.c_str());
        }
        lineage_ptr++;
    }

    for (auto sample: lineage_selected) {
        auto clade = get_clade(T, sample);
        printf("Sample: %s, Clade: %s\n", sample->identifier.c_str(), clade.c_str());
    }

    fprintf(stderr, "\n%ld Samples Selected in %ld msec \n\n", lineage_selected.size(), timer.Stop());

    timer.Start();
    const std::vector<int8_t> nuc_array{1, 2, 4, 8};
    std::vector<int> random_vec, read_positions;
    std::vector<struct pos_misread*> misread_pos;
    std::vector<std::vector<int>> sample_reads_vector; 
    std::vector<struct sample_read_pair*> sample_read_list;
    std::map<int,std::vector<struct sample_read_pair*>*> read_map;
    int num_reads = sequence_depth * ((int(ref_seq.size()) / read_length) + ((int(ref_seq.size()) % read_length) != 0));
    fprintf(stderr, "Num reads per sample: %d\n", num_reads);

    //Creating Reads by randomly selecting a seed position
    for (auto node: lineage_selected) {
        ancestors.emplace_back(T.rsearch(node->identifier, true));
        random_vec.clear();
        clock_t time = clock();
        srand(int(time));
        
        for (int i = 0; i < num_reads; i++) {
            int rand_val = int(rand() % int( (int(ref_seq.size()) - ((read_length / 2) + (read_length % 2))) - ((read_length / 2) + (read_length % 2))) ) + ((read_length / 2) + (read_length % 2));
            auto it = std::find(random_vec.begin(), random_vec.end(), rand_val);
            
            while (it != random_vec.end()) {
                rand_val = int(rand() % int( (int(ref_seq.size()) - ((read_length / 2) + (read_length % 2))) - ((read_length / 2) + (read_length % 2))) ) + ((read_length / 2) + (read_length % 2));
                it = std::find(random_vec.begin(), random_vec.end(), rand_val);
            }
            random_vec.emplace_back(rand_val);
            
            struct sample_read_pair *sample_read = new struct sample_read_pair;
            sample_read->sample = node;
            //sample_read->read = node->identifier + "_" + boost::lexical_cast<std::string>(rand_val) + "_READ_" + boost::lexical_cast<std::string>(i); 
            int start_coord = (rand_val - ((read_length / 2) + (read_length % 2)));
            int end_coord = (rand_val + (read_length / 2)) - 1;
            if (start_coord < 0)
                start_coord = 0; 
            if (end_coord >= (int(ref_seq.size())))
                end_coord = (int(ref_seq.size())) - 1;
            sample_read->read = node->identifier + "_READ_" + boost::lexical_cast<std::string>(start_coord) + "_" + boost::lexical_cast<std::string>(end_coord); 
            sample_read_list.emplace_back(sample_read);

            read_positions.clear();
            time = clock();
            srand(int(time));
            for (int pos = (rand_val - ((read_length / 2) + (read_length % 2))); pos < (rand_val + (read_length / 2)); pos ++) {
                if ((pos >= 0) && (pos < (int(ref_seq.size())))) {
                    read_positions.emplace_back(pos);
                    double rndDouble = (double)rand() / RAND_MAX;
                    if (rndDouble < (read_error[0])) {             
                        struct pos_misread *misreadpos = new struct pos_misread;
                        misreadpos->pos = pos;
                        misreadpos->read = sample_read->read;
                        misreadpos->used = false;
                        misread_pos.emplace_back(misreadpos);
                    }
                }
            }
            sample_reads_vector.emplace_back(read_positions);  
        } 
    }

    fprintf(stderr, "Misread pos = %ld\n\n", misread_pos.size()); 
    //for (auto misreads: misread_pos) 
    //    std::cout << "Read: " << misreads->read << ", Pos: " << misreads->pos << "\n";

    timer.Start();

    // for each sample, we need to make a vector of Mutation objects to store all the mutations that are present in the sample
    // we can iterate over the vector of reads and for each read, we can iterate over the vector of misread positions and check if the position is present in the read
    // if it is, then we add the mutation to the vector of mutations for the sample

    std::vector<std::vector<MAT::Mutation>> node_mutations;
    for (int i = 0; i < (int)dfs.size(); i++) {
        std::vector<MAT::Mutation> node_mut;
        node_mutations.emplace_back(node_mut);
    }

    for (int i = 0; i < (int)sample_read_list.size(); i++) {
        auto sample_read = sample_read_list[i];
        for (int j = 0; j < (int)misread_pos.size(); j++) {
            auto misread = misread_pos[j];
            if (misread->read == sample_read->read) {
                auto node = T.get_node(sample_read->sample->identifier);
                auto node_mut = node_mutations[node->dfs_idx];
                MAT::Mutation m;
                m.chrom = ref_header;
                m.position = misread->pos;
                m.ref_nuc = MAT::get_nuc_id(ref_seq[misread->pos]);
                m.par_nuc = m.ref_nuc;
                m.mut_nuc = MAT::get_nuc_id('N');
                node_mut.emplace_back(m);
                node_mutations[node->dfs_idx] = node_mut;
            }
        }
    }

    // for each node, the back mutations should be removed (i.e. those mutations where the base of the node is the same as the base of the reference)
    // we can iterate over the vector of mutations and remove those where the base is the same as the reference

    for (int i = 0; i < (int)dfs.size(); i++) {
        auto node_mut = node_mutations[i];
        auto itr = node_mut.begin();
        while (itr != node_mut.end()) {
            if (itr->ref_nuc == itr->mut_nuc) {
                itr = node_mut.erase(itr);
            } else {
                itr++;
            }
        }
        node_mutations[i] = node_mut;
    }

    fprintf(stderr, "Node Mutations in %ld msec \n\n", timer.Stop());

    // timer.Start();
 
    // fprintf(stderr, "Both Maps Filled in %ld msec \n\n", timer.Stop());

    
    // Printing and storing in the VCF
    timer.Start();

    // Now for each read, the mutations need to be written to the VCF
    // use the vector of Mutation objects we created earlier to do this
    // a 1 in the VCF indicates the presence of a mutation and a 0 indicates the absence of a mutation

    // for each read, we need to iterate over the vector of mutations and check if the mutation is present in the read
    // if it is, then we write a 1, else we write a 0

    // Open the VCF file for writing
    boost::filesystem::ofstream vcf_fileHandler(vcf_filename_reads);

    for (int i = 0; i < (int)sample_read_list.size(); i++) {
        auto sample_read = sample_read_list[i];
        auto read_mutations = node_mutations[T.get_node(sample_read->sample->identifier)->dfs_idx];
        std::vector<int> read_mutations_vector;
        for (int j = 0; j < (int)read_mutations.size(); j++) {
            auto read_mut = read_mutations[j];
            bool found = false;
            for (int k = 0; k < (int)sample_reads_vector[i].size(); k++) {
                if (read_mut.position == sample_reads_vector[i][k]) {
                    found = true;
                    break;
                }
            }
            if (found) {
                read_mutations_vector.emplace_back(1);
            } else {
                read_mutations_vector.emplace_back(0);
            }
        }

    }
}
}

void read_sample_vcf(std::vector<std::string> &vcf_samples, const std::string vcf_filename_samples) {
    boost::filesystem::ifstream fileHandler(vcf_filename_samples);
    std::string s;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            if (words[1] == "POS") {
                for (int j=9; j < (int)words.size(); j++)
                    vcf_samples.emplace_back(words[j]);
            }
        }
    }
}

void read_vcf(uint32_t num_threads, const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, std::unordered_map<int, struct read_info*> &read_map, const std::string vcf_filename_reads) {
    // Boost library used to stream the contents of the input VCF file
    // Store the header information from VCF
    timer.Start();
    std::vector<int> missing_idx;
    std::vector<struct read_info*> read_ids;
    std::string s;
    boost::filesystem::ifstream fileHandler(vcf_filename_reads);
    bool header_found = false;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            if (words[1] == "POS") {
                header_found = true;
                for (int j=9; j < (int)words.size(); j++) {
                    struct read_info * rp = new struct read_info;
                    rp->read = words[j];
                    std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                    std::smatch match;
                    if (std::regex_search(words[j], match, rgx)) {
                        rp->start = std::stoi(match[1]);
                        rp->end = std::stoi(match[2]);
                    }
                    read_ids.emplace_back(rp);
                    missing_idx.emplace_back(j);
                }
            }
            else if (header_found) {
                std::vector<std::string> alleles;
                alleles.clear();
                MAT::string_split(words[4], ',', alleles);
                int k = 0;
                while (k < (int)missing_idx.size()) {
                    size_t j = missing_idx[k];
                    auto iter = read_ids.begin();
                    std::advance(iter, k);
                    if (iter != read_ids.end()) {
                        read_map.insert({k, (*iter)});
                        MAT::Mutation m;
                        m.chrom = words[0];
                        m.position = std::stoi(words[1]);
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

    fprintf(stderr,"read VCF parsed in %ld sec\n\n", (timer.Stop() / 1000));   
}


int place_reads(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, struct read_info* rp, const MAT::Node* check_node, tbb::concurrent_hash_map<MAT::Node*, double> &node_score, const std::vector<MAT::Node*> &peak_nodes, std::vector<int> &mismatch_vector, const int par_score_lim) {
    std::stack<struct parsimony> parsimony_stack;
    struct min_parsimony min_par;
    int peak_count = 0;
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
                        if(par_node_mut.mut_nuc == node_mut.mut_nuc) {
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
                        //Mutation found in read, add to common_node_mut
                        if (read_mut.mut_nuc == node_mut.mut_nuc)
                            common_node_mut.emplace_back(read_mut);
                        //read_mut is 'N', add to common_node_mut
                        else if (read_mut.mut_nuc == 0b1111)
                            common_node_mut.emplace_back(read_mut);
                        //If mutation is different
                        else {
                            //For root, handle the mutation as child of current node. 
                            if (!i) {
                                struct MAT::Mutation new_mut;
                                new_mut.position = read_mut.position;
                                new_mut.ref_nuc = read_mut.ref_nuc;
                                new_mut.par_nuc = node_mut.mut_nuc;
                                new_mut.mut_nuc = read_mut.mut_nuc;
                                curr_node_par_mut.emplace_back(new_mut);
                                //Placing it in common_mut so don't add this mut again
                                common_node_mut.emplace_back(new_mut);
                            }
                            //Otherwise, handle it as sibling of current node, i.e. add as uniq mutation
                            //Reverse par_nuc and mut_nuc as it will again get flipped in uniq_curr_node_mut
                            else {
                                struct MAT::Mutation new_mut;
                                new_mut.position = read_mut.position;
                                new_mut.ref_nuc = read_mut.ref_nuc;
                                new_mut.par_nuc = read_mut.mut_nuc;
                                new_mut.mut_nuc = node_mut.mut_nuc;
                                uniq_curr_node_mut.emplace_back(new_mut);
                            }
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
                auto itr = common_node_mut.begin();
                while (itr != common_node_mut.end()) {
                    if (itr->position == read_mut.position) {
                        if (itr->mut_nuc != read_mut.mut_nuc)
                            std::cout << "common_node mut does not match read_mut!!!" << "\n";
                        present = true;
                        break;
                    }
                    itr++;
                }
                if (present)
                    continue;

                //Else add it to curr_node_mut if mut_nuc != 'N'
                if (read_mut.mut_nuc != 0b1111)
                    curr_node_par_mut.emplace_back(read_mut);
            }
        }

        //Updating the parsimony score of peak nodes
        auto curr_node = dfs[i];
        for (int j = 0; j < (int)peak_nodes.size(); j++) {
            if (peak_nodes[j] == curr_node) {
                mismatch_vector[j] = (int)curr_node_par_mut.size();
                peak_count++;
                if (peak_count == (int)peak_nodes.size())
                    return 0;
                break;
            }
        }

        //Checking min_parsimony
        int new_min_par = -1; 
        // If best_par_score is empty and curr_par_score >= limit -> CHANGE
        if ((!(min_par.par_list.size())) && ((int)curr_node_par_mut.size() >= par_score_lim))
            new_min_par = 1;
        // If curr_par_score < best_par_score and curr_par_score >= limit -> CHANGE
        else if ((curr_node_par_mut.size() < min_par.par_list[0].size()) && ((int)curr_node_par_mut.size() >= par_score_lim))
            new_min_par = 1;
        // If cur_par_score == best_par_score -> APPEND
        else if (curr_node_par_mut.size() == min_par.par_list[0].size())
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

        //Updating parsimony to be stored as a child
        for (auto uniq_mut: uniq_curr_node_mut) {
            //Just reverse mut_nuc and par_nuc and add it to parsimony
            int8_t temp = uniq_mut.par_nuc;
            uniq_mut.par_nuc = uniq_mut.mut_nuc;
            uniq_mut.mut_nuc = temp;
            curr_node_par_mut.emplace_back(uniq_mut);
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
        //VERIFY only if searching for nodes with zero parsimony score
        if (!par_score_lim) {
            //Check to ensure every read has its corresponding sample as its most parsimonious position
            std::string target = rp->read;
            size_t pos = target.find("_READ");
            target.erase(pos);
            bool found = false;
            int idx, i;
            for (i = 0; i < (int)min_par.idx_list.size(); i++) {
                idx = min_par.idx_list[i];
                if (dfs[idx]->identifier == target) {
                    found = true;
                    break;
                }
            }
            if ((!found) || (min_par.par_list[0].size())) {
                if (found) {
                    for (auto mut: min_par.par_list[i])
                        fprintf(stderr, "Parsimony Pos: %d, mut: %c \n", mut.position, MAT::get_nuc(mut.mut_nuc));
                    fprintf(stderr, "Sample: %s \n", dfs[idx]->identifier.c_str());
                    for (auto mut: rp->mutations)
                        fprintf(stderr, "Read mut Pos: %d, mut: %c \n", mut.position, MAT::get_nuc(mut.mut_nuc));
                }
                else {
                    fprintf(stderr, "Sample not Found !!! \n");
                    auto clade = get_clade(T, T.get_node(target));
                    fprintf(stderr, "Target: %s, Clade: %s\n", target.c_str(), clade.c_str());
                    fprintf(stderr, "mut pos: ");
                    for (auto anc: T.rsearch(target, true)) { //Checking all ancestors of a node to get clade 
                        for (auto mut: anc->mutations)
                            fprintf(stderr, "%c%d%c, ", MAT::get_nuc(mut.ref_nuc), mut.position, MAT::get_nuc(mut.mut_nuc));
                    }
                    fprintf(stderr, "\n");
                }
                fprintf(stderr, "Read: %s, read mutations: %ld, Parsimony score = %ld, parsimonious positions: %ld\n\n", rp->read.c_str(), rp->mutations.size(), min_par.par_list[0].size(), min_par.par_list.size());
                std::cout << "\n";
            }
        }
        //Keeping tab on weighted read score being mapped to each parsimonious node
        min_par.par_list.clear();
        // In absence of specific clades, add scores to nodes
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


void analyze_reads(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, const std::unordered_map<int, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score, const std::vector<std::string> &vcf_samples, const std::string &mismatch_matrix_file, const std::string &barcode_file, const std::string &read_abundance_vcf) {
    timer.Start();
    int top_n = 25, m_dist_thresh = 2, remaining_read_thresh = (int)read_map.size() * 0.0025;
    std::vector<std::pair<MAT::Node*, double>> top_n_node_score(top_n);
    std::vector<MAT::Node*> peak_nodes_dummy, peak_nodes;
    std::vector<int> remaining_reads, mismatch_vector_dummy;
    
    //GREEDY ALGORITHM for getting Lineages
    //Initializing remaining_reads
    auto itr = read_map.begin();
    while (itr != read_map.end()) {
        remaining_reads.emplace_back(itr->first);
        itr++;
    }
    while ((int)remaining_reads.size() > remaining_read_thresh) {
        //Calculating node score for remaining reads
        static tbb::affinity_partitioner ap;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    int rm_idx = remaining_reads[i];
                    auto read_id = read_map.find(rm_idx)->second;
                    place_reads(T, dfs, read_id, NULL, node_score, peak_nodes_dummy, mismatch_vector_dummy, 0);
                }
            },
        ap);
        //Sorting the node scores
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

        //Find top node not seen before
        std::vector<bool> peak_vec(top_n_node_score.size(), true);
        auto top_n_itr = top_n_node_score.begin();
        //Find top_node that is NOT present in neighborhood of peak_nodes
        while (top_n_itr != top_n_node_score.end()) {
            bool present = check_peaks_neighbourhood(top_n_itr->first, peak_nodes, m_dist_thresh);
            if (!present)
                break;
            peak_vec[top_n_itr - top_n_node_score.begin()] = false;
            top_n_itr++;
        }
        auto top_score = top_n_itr->second;
        
        //Find the reads not mapped to the best nodes
        while (top_n_itr != top_n_node_score.end()) {
            //Peak Finding
            auto curr_node = top_n_itr->first;
            //Only add unique nodes in peak_nodes with score == top_score
            if (abs(top_score - top_n_itr->second) < 1e-9) {
                //Check if curr_node does not lie in neighborhood of any node seen in current iterarion or previous
                bool present = check_peaks_neighbourhood(curr_node, peak_nodes, m_dist_thresh);
                if ((peak_vec[top_n_itr - top_n_node_score.begin()]) && (!present)) {
                    peak_nodes.emplace_back(curr_node);
                    //Remove mapped reads from remaining_reads
                    auto clade = get_clade(T, curr_node);
                    printf("PEAK Node: %s, Clade: %s\n", curr_node->identifier.c_str(), clade.c_str());
                    std::vector<int> remove_reads;
                    using my_mutex_t = tbb::queuing_mutex;
                    my_mutex_t my_mutex;
                    static tbb::affinity_partitioner ap;
                    //Parallel_for loop for each remaining read
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
                        [&](tbb::blocked_range<size_t> k) {
                            for (size_t i = k.begin(); i < k.end(); ++i) {
                                int rm_idx = remaining_reads[i];
                                auto read_id = read_map.find(rm_idx)->second;
                                int read_present = place_reads(T, dfs, read_id, curr_node, node_score, peak_nodes_dummy, mismatch_vector_dummy, 0);
                                node_score.clear();
                                //Read contains current node -> REMOVE
                                if (read_present) {
                                    my_mutex_t::scoped_lock my_lock{my_mutex};
                                    remove_reads.emplace_back(rm_idx);
                                }
                            }
                        },
                    ap);
                    //Erase from remove_reads
                    for (auto rm_idx: remove_reads) {
                        auto itr = remaining_reads.begin();
                        while (itr != remaining_reads.end()) {
                            if ((*itr) == rm_idx) {
                                remaining_reads.erase(itr);
                                break;
                            }
                            itr++;
                        }
                    }
                    remove_reads.clear();
                }
                //If present in neighbourhood, move to next node
                else {
                    peak_vec[top_n_itr - top_n_node_score.begin()] = false;
                    top_n_itr++;
                    continue;
                }
            }
            else if (abs(top_score - top_n_itr->second) > 1e-9)
                break;
            
            //Remove nearby peaks from further analysis in this iteration
            auto top_n_peak_cmp_itr = top_n_itr + 1;
            while (top_n_peak_cmp_itr != top_n_node_score.end()) {
                auto cmp_node = top_n_peak_cmp_itr->first;
                // Don't check if score < top_score or if node lies in neighborhood
                if (abs(top_score - top_n_peak_cmp_itr->second) > 1e-9)
                    break;
                else if (!peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()]) {
                    top_n_peak_cmp_itr++;
                    continue;
                }
                //Don't consider peaks within mutation distance limit 
                int m_dist = mutation_distance(curr_node, cmp_node);
                if (m_dist > m_dist_thresh)
                    peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = true;
                else
                    peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = false;
                top_n_peak_cmp_itr++;
            }
            top_n_itr++;
        }        
        top_n_node_score.clear();
        top_n_node_score.resize(top_n);
        peak_vec.clear();
        std::cout << "\n";
    }
    std::cout << "\nRemaining reads: " << remaining_reads.size() << "\n";
    remaining_reads.clear();
    fprintf(stderr,"Peak search took %ld min\n\n", (timer.Stop() / (60 * 1000)));
    
    //Verify Recovery of Input Samples
    timer.Start();
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto ref_clade = get_clade(T, T.get_node(sample));
        auto best_clade = ref_clade;
        auto best_node = T.get_node(sample);
        for (auto n_s: peak_nodes) {
            int curr_dist = mutation_distance(n_s, T.get_node(sample));
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_clade = get_clade(T, n_s);
                best_node = n_s;
            }
        }
        printf("Node: %s, Closest_node: %s, Ref_clade: %s, closest_clade: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), ref_clade.c_str(), best_clade.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld msec\n\n", timer.Stop());

    generate_regression_abundance_data(T, peak_nodes, read_map, barcode_file, read_abundance_vcf);
    //generate_EM_data(T, dfs, read_map, peak_nodes, mismatch_matrix_file);
}


//Function to check if current node is in neighbourhood of peak nodes
bool check_peaks_neighbourhood (const MAT::Node* N, const std::vector<MAT::Node*> &peak_nodes, const int m_dist_thresh) {
    for (auto pn: peak_nodes) {
        //Return false if the Node lies within mutation distance limit 
        int m_dist = mutation_distance(pn, N);
        if (m_dist <= m_dist_thresh)
            return true;
    }
    return false;
}

//Function to calculation distance between two nodes
int mutation_distance(const MAT::Node* N1, const MAT::Node* N2) {
    if (N1 == N2)
        return 0;
    std::vector<std::pair<MAT::Mutation, bool>> mutation_list;
    std::vector<MAT::Mutation> removed_muts, dont_remove_muts;
    while (N1->level > N2->level) {
        update_unique_mutations(N1, mutation_list, removed_muts, dont_remove_muts, 0);
        N1 = N1->parent;
    } 
    while (N1->level < N2->level) {
        update_unique_mutations(N2, mutation_list, removed_muts, dont_remove_muts, 1);
        N2 = N2->parent;
    }
    while (N1->parent != N2->parent) {
        update_unique_mutations(N1, mutation_list, removed_muts, dont_remove_muts, 0);
        update_unique_mutations(N2, mutation_list, removed_muts, dont_remove_muts, 1);
        N1 = N1->parent;
        N2 = N2->parent;
    }
    update_unique_mutations(N1, mutation_list, removed_muts, dont_remove_muts, 0);
    update_unique_mutations(N2, mutation_list, removed_muts, dont_remove_muts, 1);
    return (int)mutation_list.size();
}

void update_unique_mutations(const MAT::Node* N, std::vector<std::pair<MAT::Mutation, bool>> &mutation_list, std::vector<MAT::Mutation> &removed_muts, std::vector<MAT::Mutation> &dont_remove_muts, bool is_N2) {
    //Iterating through mutations of Node N
    for (auto mut: N->mutations) {
        //If present in removed (common mutation) then don't check for this mutation
        auto rm_itr = removed_muts.begin();
        while (rm_itr != removed_muts.end()) {
            if (rm_itr->position == mut.position)
                break;
            rm_itr++;
        }
        if (rm_itr != removed_muts.end())
            continue;
        //Iterating through mutations already stored in mutation_list
        auto mut_itr = mutation_list.begin();
        while (mut_itr != mutation_list.end()) {
            if ((mut_itr->first.position == mut.position) && (mut_itr->first.mut_nuc == mut.mut_nuc)) {
                //Don't remove if different mutation is found near leaf node
                auto dr_itr = dont_remove_muts.begin();
                while (dr_itr != dont_remove_muts.end()) {
                    if (dr_itr->position == mut.position)
                        break;
                    dr_itr++;
                }
                if (dr_itr == dont_remove_muts.end()) {
                    mutation_list.erase(mut_itr);
                    removed_muts.emplace_back(mut_itr->first);
                }
                break;
            }
            else if (mut_itr->first.position == mut.position) {
                //If different mutation is near leaf of other node then can't remove. Should belong to other node
                if (mut_itr->second != is_N2) {
                    auto dr_itr = dont_remove_muts.begin();
                    while (dr_itr != dont_remove_muts.end()) {
                        if (dr_itr->position == mut.position)
                            break;
                        dr_itr++;
                    }
                    if (dr_itr == dont_remove_muts.end())
                        dont_remove_muts.emplace_back(mut_itr->first);
                }
                break;
            }
            mut_itr++;
        }
        if (mut_itr == mutation_list.end())
            mutation_list.emplace_back(std::pair(mut, is_N2));
    }
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

//Generate EM based estimate algorithm data
void generate_EM_data(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, const std::unordered_map<int, struct read_info*> &read_map, const std::vector<MAT::Node*> &peak_nodes, const std::string &mismatch_matrix_file) {
    //Write MISMATCH File for Abundance Estimation using EM algorithm
    timer.Start();
    std::unordered_map<int, std::vector<int>> mismatch_matrix;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    std::ofstream outfile_mismatch_matrix(mismatch_matrix_file, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_mismatch_matrix;
    if (mismatch_matrix_file.find(".gz\0") != std::string::npos) {
            outbuf_mismatch_matrix.push(boost::iostreams::gzip_compressor());
    }
    outbuf_mismatch_matrix.push(outfile_mismatch_matrix);
    std::ostream mismatch_file(&outbuf_mismatch_matrix);

    //Parallel_for loop for each remaining read
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_map.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t i = k.begin(); i < k.end(); ++i) {
                std::vector<int> mismatch_vector(peak_nodes.size(), -1);
                auto read_id = read_map.find(i)->second;
                place_reads(T, dfs, read_id, NULL, node_score, peak_nodes, mismatch_vector, 0);
                node_score.clear();
                my_mutex_t::scoped_lock my_lock{my_mutex};
                mismatch_matrix.insert({i, mismatch_vector});
            }
        },
    ap);
    
    auto m_itr = mismatch_matrix.begin();
    while (m_itr != mismatch_matrix.end()) {
        auto m_element = m_itr->second;
        for (int i = 0; i < (int)m_element.size(); i++) {
            if (i)
                mismatch_file << ",";
            mismatch_file << m_element[i];
        }
        mismatch_file << "\n";
        m_itr++;
    }
    boost::iostreams::close(outbuf_mismatch_matrix);
    outfile_mismatch_matrix.close();
    fprintf(stderr,"Mismatch matrix file writing took %ld min\n\n", (timer.Stop() / (60*1000)));
}

//Generate regression based estimate algorithm data
void generate_regression_abundance_data(const MAT::Tree &T, const std::vector<MAT::Node*> &peak_nodes, const std::unordered_map<int, struct read_info*> &read_map, const std::string &barcode_file, const std::string &read_abundance_vcf) {
    timer.Start();

    //Get Mutations of Peaks for lineage estimation
    std::unordered_map<MAT::Node*, std::vector<MAT::Mutation>> peak_mut_map;
    //Create a map of peak node mutations 
    for (auto pk: peak_nodes) {
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
        peak_mut_map.insert({pk, mut_list});
        mut_list.clear();
    }

    //Remove the mutations that are common to all peak nodes
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
    std::ofstream outfile_vcf(read_abundance_vcf, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_barcode, outbuf_vcf;
    if (barcode_file.find(".gz\0") != std::string::npos) {
            outbuf_barcode.push(boost::iostreams::gzip_compressor());
    }
    if (read_abundance_vcf.find(".gz\0") != std::string::npos) {
            outbuf_vcf.push(boost::iostreams::gzip_compressor());
    }
    outbuf_barcode.push(outfile_barcode);
    outbuf_vcf.push(outfile_vcf);
    std::ostream barcode(&outbuf_barcode);
    std::ostream vcf(&outbuf_vcf);

    //Storing unique mutations in a vector
    std::vector<MAT::Mutation> peak_mut_list;
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

    std::sort(peak_mut_list.begin(), peak_mut_list.end(), compare_mutations);
    //Writing the header mutations
    std::string barcode_print, vcf_print;
    for (auto mut: peak_mut_list) {
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
        float af = (float)match_reads / (float)total_reads;
        vcf_print += "NC_045512v2\t" + std::to_string(mut.position) + "\t" + MAT::get_nuc(mut.par_nuc) + std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc) + "\t" + MAT::get_nuc(mut.par_nuc) + "\t" + MAT::get_nuc(mut.mut_nuc) + "\t.\t.\tAF=";
        vcf_print += std::to_string(af) + "\n";
    } 
    barcode << barcode_print;
    vcf << vcf_print;
    barcode_print.clear();
    vcf_print.clear();

    //Writing peak mutations
    peak_mut_itr = peak_mut_map.begin();
    while (peak_mut_itr != peak_mut_map.end()) {
        barcode_print += "\n";
        barcode_print += peak_mut_itr->first->identifier.c_str();
        std::vector<int> mut_idx_list;
        //Iterating through mutations of this peak
        auto mut_itr = peak_mut_itr->second.begin();
        while (mut_itr != peak_mut_itr->second.end()) {
            //Search for mutation in the list
            auto pm_itr = peak_mut_list.begin();
            while (pm_itr != peak_mut_list.end()) {
                if ((pm_itr->position == mut_itr->position) && (pm_itr->mut_nuc == mut_itr->mut_nuc))
                    break;
                pm_itr++;                
            }
            mut_idx_list.emplace_back((int)(pm_itr - peak_mut_list.begin()));
            mut_itr++;
        }
        //Sort the mutation indexes
        std::sort(mut_idx_list.begin(), mut_idx_list.end());
        int idx = 0;
        for (auto m_idx: mut_idx_list) {
            while (idx < m_idx) {
                barcode_print += ",0.0";
                idx++;
            }
            if (idx == m_idx) {
                barcode_print += ",1.0";
                idx++;
            }
        }
        while (idx < (int)peak_mut_list.size()) {
            barcode_print += ",0.0";
            idx++;
        }
        mut_idx_list.clear();
        barcode << barcode_print;
        barcode_print.clear();
        peak_mut_itr++;
    }
        
    fprintf(stderr,"Barcode file writing took %ld sec\n\n", (timer.Stop() / 1000));
}


//Comparing mutations for sorting a vector 
bool compare_mutations(const MAT::Mutation &a, const MAT::Mutation &b) {
    return a.position < b.position;
}