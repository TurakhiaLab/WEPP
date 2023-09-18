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
                if (start_coord < 1)
                    start_coord = 1; 
                if (end_coord >= (int(ref_seq.size())))
                    end_coord = (int(ref_seq.size()));
                sample_read->read = node->identifier + "_READ_" + boost::lexical_cast<std::string>(start_coord) + "_" + boost::lexical_cast<std::string>(end_coord); 
                sample_read_list.emplace_back(sample_read);

                read_positions.clear();
                time = clock();
                srand(int(time));
                for (int pos = (rand_val - ((read_length / 2) + (read_length % 2))); pos < (rand_val + (read_length / 2)); pos ++) {
                    if ((pos > 0) && (pos < (int(ref_seq.size()+1)))) {
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
        
        
        //Creating read map to place the reads acc to mut positions
        std::vector<struct sample_read_pair*>::iterator read_name_ptr;
        read_name_ptr = sample_read_list.begin();
        for (auto reads: sample_reads_vector) {
            for (auto pos: reads) {
                if (read_map.find(pos) == read_map.end()) {
                    std::vector<struct sample_read_pair*> *sr_list = new std::vector<struct sample_read_pair*>;
                    sr_list->emplace_back(*read_name_ptr);
                    read_map.insert({pos, sr_list});
                }
                else {
                    std::vector<struct sample_read_pair*> *sr_list;
                    sr_list = read_map[pos];
                    auto itr = std::find(sr_list->begin(), sr_list->end(), *read_name_ptr);
                    if (itr == sr_list->end())
                        sr_list->emplace_back(*read_name_ptr);
                }
            }
            read_name_ptr++;
        }


        // Inserting selected samples in the Map
        for (auto anc: ancestors) {
            for (auto node: anc) {
                for (auto mut: node->mutations){
                    bool back_mutation = true;
                    if (back_mut_map.find(mut.position) == back_mut_map.end()) {
                        if (mut.ref_nuc != mut.mut_nuc)
                            back_mutation =  false;
                    }
                    // No Back Mutation
                    if (!back_mutation) {
                        if(sample_map.find(mut.position) == sample_map.end()){
                            std::vector<struct ances_sample_list*> *anc_sample_list = new std::vector<struct ances_sample_list*>;
                            std::vector<Mutation_Annotated_Tree::Node*> *samples = new std::vector<Mutation_Annotated_Tree::Node*>;
                            struct ances_sample_list *anc_nodes = new ances_sample_list;

                            anc_nodes->ancestor_node = node;
                            samples->emplace_back(anc[0]);
                            anc_nodes->sample_nodes = samples;
                            anc_sample_list->emplace_back(anc_nodes);
                            sample_map.insert({mut.position, anc_sample_list});
                        }
                        else {
                            std::vector<struct ances_sample_list*> *anc_sample_list; 
                            std::vector<struct ances_sample_list*>::iterator ptr; 
                            anc_sample_list = sample_map[mut.position];

                            for (ptr = anc_sample_list->begin(); ptr < anc_sample_list->end(); ptr++) 
                                if ((*ptr)->ancestor_node == node)
                                    break;

                            if (ptr == anc_sample_list->end()) {
                                struct ances_sample_list *anc_nodes = new ances_sample_list;
                                std::vector<Mutation_Annotated_Tree::Node*> *samples = new std::vector<Mutation_Annotated_Tree::Node*>;
                                anc_nodes->ancestor_node = node;
                                samples->emplace_back(anc[0]);
                                anc_nodes->sample_nodes = samples;
                                anc_sample_list->emplace_back(anc_nodes);
                            }
                            else {
                                bool present = false;
                                for (auto sample: *(*ptr)->sample_nodes)
                                    if (sample == anc[0]){
                                        present = true;
                                        break;
                                    }
                                if (!present) {
                                    (*ptr)->sample_nodes->emplace_back(anc[0]);
                                }
                            }
                        }
                    }   
                    // First Back Mutation
                    else if (back_mut_map.find(mut.position) == back_mut_map.end()) {
                        if(sample_map.find(mut.position) != sample_map.end()){
                            std::vector<struct ances_sample_list*> *anc_sample_list; 
                            anc_sample_list = sample_map[mut.position];
                            for (auto ances: *anc_sample_list) {
                                // New Back mutation is not ancestor of ancestor node from sample_map 
                                if (T.is_ancestor(ances->ancestor_node->identifier, node->identifier)) {
                                    for (auto leaf: *ances->sample_nodes) {
                                        if (T.is_ancestor(node->identifier, leaf->identifier) || (node->identifier == leaf->identifier)) {
                                            ances->sample_nodes->erase(std::remove(ances->sample_nodes->begin(), ances->sample_nodes->end(), leaf), ances->sample_nodes->end());
                                        }
                                    }
                                }
                                // New Back Mutation node is either ancestor or sibling of previous node in map, i.e., BM is irrelevant
                                // Do nothing as iteration is for current node only
                                if (ances->sample_nodes->size() < 1)
                                    anc_sample_list->erase(std::remove(anc_sample_list->begin(), anc_sample_list->end(), ances), anc_sample_list->end());
                            }
                        }

                        std::vector<Mutation_Annotated_Tree::Node*> *back_mut_list = new std::vector<Mutation_Annotated_Tree::Node*>;
                        back_mut_list->emplace_back(node);
                        back_mut_map.insert({mut.position, back_mut_list});
                    } 
                    // Back Mutation present at that position
                    else {
                        auto back_mut_list = back_mut_map[mut.position];
                        int no_BM = 0;
                        //New node without Back Mutation
                        if (mut.ref_nuc != mut.mut_nuc) {
                            for (auto bm_node: *back_mut_list) {
                                // If new Node is ancestor of Back Mutation Node and sample does not belong to BM 
                                if ((T.is_ancestor(node->identifier, bm_node->identifier)) && (!( (bm_node->identifier == anc[0]->identifier) || (T.is_ancestor(bm_node->identifier, anc[0]->identifier) )))) {
                                    no_BM += 1;
                                }
                                // If Back Mutation node is ancestor of new node
                                else if (T.is_ancestor(bm_node->identifier, node->identifier)) {
                                    no_BM += 1;
                                }
                                // If Back Mutation node is sibling of new node
                                else if (!( (T.is_ancestor(bm_node->identifier, node->identifier)) || (T.is_ancestor(node->identifier, bm_node->identifier)) )) {
                                    no_BM += 1;
                                }
                            }
                            if (no_BM == int(back_mut_list->size())) {
                                if(sample_map.find(mut.position) == sample_map.end()){
                                    std::vector<struct ances_sample_list*> *anc_sample_list = new std::vector<struct ances_sample_list*>;
                                    std::vector<Mutation_Annotated_Tree::Node*> *samples = new std::vector<Mutation_Annotated_Tree::Node*>;
                                    struct ances_sample_list *anc_nodes = new ances_sample_list;

                                    anc_nodes->ancestor_node = node;
                                    samples->emplace_back(anc[0]);
                                    anc_nodes->sample_nodes = samples;
                                    anc_sample_list->emplace_back(anc_nodes);
                                    sample_map.insert({mut.position, anc_sample_list});
                                }
                                else {
                                    std::vector<struct ances_sample_list*> *anc_sample_list; 
                                    std::vector<struct ances_sample_list*>::iterator ptr; 
                                    anc_sample_list = sample_map[mut.position];

                                    for (ptr = anc_sample_list->begin(); ptr < anc_sample_list->end(); ptr++) 
                                        if ((*ptr)->ancestor_node == node)
                                            break;

                                    if (ptr == anc_sample_list->end()) {
                                        struct ances_sample_list *anc_nodes = new ances_sample_list;
                                        std::vector<Mutation_Annotated_Tree::Node*> *samples = new std::vector<Mutation_Annotated_Tree::Node*>;
                                        anc_nodes->ancestor_node = node;
                                        samples->emplace_back(anc[0]);
                                        anc_nodes->sample_nodes = samples;
                                        anc_sample_list->emplace_back(anc_nodes);
                                    }
                                    else {
                                        bool present = false;
                                        for (auto sample: *(*ptr)->sample_nodes)
                                            if (sample == anc[0]){
                                                present = true;
                                                break;
                                            }
                                        if (!present) {
                                            (*ptr)->sample_nodes->emplace_back(anc[0]);
                                        }
                                    }
                                }
                            }
                        }
                        //New node with Back Mutation
                        else {
                            std::vector<struct ances_sample_list*> *anc_sample_list; 
                            if (sample_map.find(mut.position) != sample_map.end()) {
                                anc_sample_list = sample_map[mut.position];
                                for (auto ances: *anc_sample_list) {
                                    // New Back mutation is not ancestor of node from sample_map
                                    if (T.is_ancestor(ances->ancestor_node->identifier, node->identifier)) {
                                        for (auto leaf: *ances->sample_nodes) {
                                            if (T.is_ancestor(node->identifier, leaf->identifier) || (node->identifier == leaf->identifier))
                                                ances->sample_nodes->erase(std::remove(ances->sample_nodes->begin(), ances->sample_nodes->end(), leaf), ances->sample_nodes->end());
                                        }
                                        if (ances->sample_nodes->size() < 1)
                                            anc_sample_list->erase(std::remove(anc_sample_list->begin(), anc_sample_list->end(), ances), anc_sample_list->end());
                                    }
                                    // New Back Mutation is either ancestor or sibling of node in map i.e., BM is irrelevant
                                    // Do nothing as iteration is for current node only
                                }
                            }
                            back_mut_list->emplace_back(node);
                            back_mut_map.insert({mut.position, back_mut_list});
                        }
                    }
                }
            }
        }

        ancestors.clear();
        fprintf(stderr, "Both Maps Filled in %ld msec \n\n", timer.Stop());

        timer.Start();
        //Printing and storing in VCF
        std::ofstream outfile_samples(vcf_filename_samples, std::ios::out | std::ios::binary);
        std::ofstream outfile_reads(vcf_filename_reads, std::ios::out | std::ios::binary);
        std::ofstream outfile_reads_freyja(vcf_filename_reads_freyja, std::ios::out | std::ios::binary);
        std::ofstream outfile_reads_freyja_depth(depth_filename_reads_freyja, std::ios::out | std::ios::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_samples, outbuf_reads, outbuf_reads_freyja, outbuf_reads_freyja_depth;
        if (vcf_filename_samples.find(".gz\0") != std::string::npos) {
            outbuf_samples.push(boost::iostreams::gzip_compressor());
        }
        if (vcf_filename_reads.find(".gz\0") != std::string::npos) {
            outbuf_reads.push(boost::iostreams::gzip_compressor());
        }
        if (vcf_filename_reads_freyja.find(".gz\0") != std::string::npos) {
            outbuf_reads_freyja.push(boost::iostreams::gzip_compressor());
        }
        if (depth_filename_reads_freyja.find(".gz\0") != std::string::npos) {
            outbuf_reads_freyja_depth.push(boost::iostreams::gzip_compressor());
        }
        outbuf_samples.push(outfile_samples);
        outbuf_reads.push(outfile_reads);
        outbuf_reads_freyja.push(outfile_reads_freyja);
        outbuf_reads_freyja_depth.push(outfile_reads_freyja_depth);
        std::ostream vcf_file_samples(&outbuf_samples);
        std::ostream vcf_file_reads(&outbuf_reads);
        std::ostream vcf_file_reads_freyja(&outbuf_reads_freyja);
        std::ostream depth_file_reads_freyja(&outbuf_reads_freyja_depth);

        vcf_file_samples << "##fileformat=VCFv4.2\n";
        vcf_file_samples << "##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n";
        vcf_file_samples << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        
        vcf_file_reads << "##fileformat=VCFv4.2\n";
        vcf_file_reads << "##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n";
        vcf_file_reads << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        
        vcf_file_reads_freyja << "##fileformat=VCFv4.2\n";
        vcf_file_reads_freyja << "##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n";
        vcf_file_reads_freyja << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
        
        for (auto node: lineage_selected) {
            vcf_file_samples << "\t" << node->identifier;
        }

        for (auto sample: sample_read_list) {
            vcf_file_reads << "\t" << sample->read;
        }
        
        vcf_file_samples << "\n";
        vcf_file_reads << "\n";
        vcf_file_reads_freyja << "\n";
        
        std::vector<int> lineage_present(lineage_selected.size());
        std::vector<int> read_present(sample_read_list.size());
        std::vector<int> misreads_eligible;
        std::string vcf_file_read_holder, vcf_file_read_match, vcf_file_read_freyja_holder, vcf_file_read_freyja_match, depth_file_read_freyja_holder;
        std::vector<struct pos_misread*>::iterator misread_ptr;
        std::vector<struct sample_ances_list*> sample_ancestors;
        auto previous_pos = sample_map.begin()->first;
        int map_count = 0, allele_pos = 0;
        
        
        //Populating VCF based on mutating positions in the sample map 
        for (auto map: sample_map){
            bool encountered_sample_ancestor = false, encountered_read = false;
            std::string mut_nuc_list_samples, mut_nuc_list_reads;
            std::vector<struct read_mut_pair*> misread_names;
            std::vector<int8_t>nuc_used, nuc_available;
            char ref_nuc_name_samples = {};
            char ref_nuc_name_reads = {};
            mut_nuc_list_samples.clear(); 
            mut_nuc_list_reads.clear();
            misreads_eligible.clear();
            sample_ancestors.clear(); 
            vcf_file_read_holder.clear();
            vcf_file_read_match.clear();
            vcf_file_read_freyja_holder.clear();
            vcf_file_read_freyja_match.clear();
            fill(lineage_present.begin(), lineage_present.end(), 0);
            fill(read_present.begin(), read_present.end(), 0);
            std::vector<int>read_match_count;
    
            //Tackling Misreads not present as mutations in our sample map
            misread_ptr = misread_pos.begin();
            while (misread_ptr != misread_pos.end()) {
                auto m_e_itr = std::find(misreads_eligible.begin(), misreads_eligible.end(), (*misread_ptr)->pos);
                if ( (!(*misread_ptr)->used) && (m_e_itr == misreads_eligible.end()) && ( ( (map_count == (int)sample_map.size()-1) && ((*misread_ptr)->pos > map.first) ) || ( (map.first == sample_map.begin()->first) && ((*misread_ptr)->pos < map.first) ) || ( (map.first != sample_map.begin()->first) && ((*misread_ptr)->pos < map.first) && ((*misread_ptr)->pos > previous_pos) ) ) ) {
                    misreads_eligible.emplace_back((*misread_ptr)->pos);    
                }
                misread_ptr++;
            }
            std::sort(misreads_eligible.begin(), misreads_eligible.end());
            int last_pos = -1;
            misread_names.clear();

            for (auto misread: misreads_eligible) {
                if (last_pos != -1) {
                    fill(read_present.begin(), read_present.end(), 0);
                    vcf_file_read_holder.append("\t");
                    vcf_file_read_holder.push_back(ref_nuc_name_reads);
                    vcf_file_read_holder.append("\t");
                    
                    vcf_file_read_freyja_holder.append("\t");
                    vcf_file_read_freyja_holder.push_back(ref_nuc_name_reads);
                    vcf_file_read_freyja_holder.append("\t");
                    read_match_count.clear();
                    for (int i = 0; i < int(mut_nuc_list_reads.size()); i++) {
                        if (i) {
                            vcf_file_read_holder += ",";
                            vcf_file_read_freyja_holder += ",";
                        }
                        vcf_file_read_holder += mut_nuc_list_reads[i];
                        vcf_file_read_freyja_holder += mut_nuc_list_reads[i];
                        read_match_count.emplace_back(0);
                    }
                    vcf_file_read_holder += "\t.\t.\t.\t.\t";
                    vcf_file_read_freyja_holder += "\t.\t.\t";

                    std::vector<struct sample_read_pair*>::iterator s_r_ptr;
                    s_r_ptr = sample_read_list.begin();
                    for (auto read_mut: misread_names) {
                        std::string read_name = read_mut->read_name;
                        while (s_r_ptr != sample_read_list.end()) {
                            if ((*s_r_ptr)->read == read_name) {
                                read_present[s_r_ptr - sample_read_list.begin()] = mut_nuc_list_reads.find(read_mut->mut_nuc) + 1;
                                break;
                            }
                            s_r_ptr++;
                        } 
                    }

                    for (auto read: read_present) {
                        vcf_file_read_holder.append(std::to_string(read) + "\t");
                        if (read)
                            read_match_count[read-1] += 1;
                    }
                    vcf_file_read_holder += "\n";
                    
                    int total_count = 0;
                    s_r_ptr = sample_read_list.begin();
                    int start, end;
                    while (s_r_ptr != sample_read_list.end()) {
                        std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                        std::smatch match;
                        if (std::regex_search((*s_r_ptr)->read, match, rgx)) {
                            start = std::stoi(match[1]);
                            end = std::stoi(match[2]);
                            if ((misread >= start) && (misread <= end))
                                total_count++;
                        }
                        s_r_ptr++;
                    }
                    vcf_file_read_freyja_holder += "AF=";
                    for (int i = 0; i < int(read_match_count.size()); i++) {
                        if (i)
                            vcf_file_read_freyja_holder += ",";
                        float af = (float)read_match_count[i] / (float)total_count;
                        vcf_file_read_freyja_holder.append(std::to_string(af));
                    }
                    vcf_file_read_freyja_holder += "\n";
                    
                    ref_nuc_name_reads={};
                    misread_names.clear();
                    mut_nuc_list_reads.clear();
                }

                misread_ptr = misread_pos.begin();
                while (misread_ptr != misread_pos.end()) {
                    if ((misread == (*misread_ptr)->pos) && (!(*misread_ptr)->used)) {
                        if (misread != last_pos) {
                            nuc_available.clear();
                            nuc_used.clear();
                            nuc_used.emplace_back(MAT::get_nuc_id(ref_seq[(*misread_ptr)->pos-1]));
                            std::vector<int8_t> misread_nuc;
                            struct read_mut_pair *read_mut = new struct read_mut_pair;

                            for (auto nuc_in_use: nuc_used) {
                                for (auto nuc: nuc_array){
                                    if (nuc_in_use != nuc) {
                                        nuc_available.emplace_back(nuc);
                                    }
                                }
                            }
                            std::sample(
                                nuc_available.begin(),
                                nuc_available.end(),
                                std::back_inserter(misread_nuc),
                                1,
                                std::mt19937{std::random_device{}()}
                            );
                            ref_nuc_name_reads = ref_seq[(*misread_ptr)->pos-1];
                            read_mut->read_name = (*misread_ptr)->read;
                            read_mut->mut_nuc = MAT::get_nuc(misread_nuc[0]);
                            misread_names.emplace_back(read_mut);
                            
                            vcf_file_read_holder += "NC_045512v2\t" + std::to_string((*misread_ptr)->pos) + "\t";
                            vcf_file_read_holder += ref_nuc_name_reads + std::to_string((*misread_ptr)->pos) + MAT::get_nuc(misread_nuc[0]);
                            vcf_file_read_freyja_holder += "NC_045512v2\t" + std::to_string((*misread_ptr)->pos) + "\t";
                            vcf_file_read_freyja_holder += ref_nuc_name_reads + std::to_string((*misread_ptr)->pos) + MAT::get_nuc(misread_nuc[0]);
                            mut_nuc_list_reads.push_back(MAT::get_nuc(misread_nuc[0])); 
                            nuc_used.emplace_back(misread_nuc[0]);
                            misread_nuc.clear();
                        }
                        else {
                            std::vector<int8_t> misread_nuc;
                            struct read_mut_pair *read_mut = new struct read_mut_pair;
                            bool unique = true;
                            std::sample(
                                nuc_available.begin(),
                                nuc_available.end(),
                                std::back_inserter(misread_nuc),
                                1,
                                std::mt19937{std::random_device{}()}
                            );
                            for (auto nuc_not_avail: nuc_used) {
                                if (misread_nuc[0] == nuc_not_avail)
                                    unique = false;
                            }
                            if (unique) {
                                vcf_file_read_holder += ",";
                                vcf_file_read_holder += ref_nuc_name_reads + std::to_string((*misread_ptr)->pos) + MAT::get_nuc(misread_nuc[0]); 
                                vcf_file_read_freyja_holder += ",";
                                vcf_file_read_freyja_holder += ref_nuc_name_reads + std::to_string((*misread_ptr)->pos) + MAT::get_nuc(misread_nuc[0]); 
                                mut_nuc_list_reads.push_back(MAT::get_nuc(misread_nuc[0]));  
                                nuc_used.emplace_back(misread_nuc[0]);
                            }
                            read_mut->read_name = (*misread_ptr)->read;
                            read_mut->mut_nuc = MAT::get_nuc(misread_nuc[0]);
                            misread_names.emplace_back(read_mut);
                            misread_nuc.clear();
                        }
                        (*misread_ptr)->used = true;
                        last_pos = misread;
                    }   
                    misread_ptr++;
                }
            }

            if (misreads_eligible.size()) {
                vcf_file_read_holder.append("\t");
                vcf_file_read_holder.push_back(ref_nuc_name_reads);
                vcf_file_read_holder.append("\t");
                
                vcf_file_read_freyja_holder.append("\t");
                vcf_file_read_freyja_holder.push_back(ref_nuc_name_reads);
                vcf_file_read_freyja_holder.append("\t");
                read_match_count.clear();
                for (int i = 0; i < int(mut_nuc_list_reads.size()); i++) {
                    if (i) {
                        vcf_file_read_holder += ",";
                        vcf_file_read_freyja_holder += ",";
                    }
                    vcf_file_read_holder += mut_nuc_list_reads[i];
                    vcf_file_read_freyja_holder += mut_nuc_list_reads[i];
                    read_match_count.emplace_back(0);
                }
                vcf_file_read_holder += "\t.\t.\t.\t.\t";
                vcf_file_read_freyja_holder += "\t.\t.\t";
                
                fill(read_present.begin(), read_present.end(), 0);
                std::vector<struct sample_read_pair*>::iterator s_r_ptr;
                s_r_ptr = sample_read_list.begin();
                for (auto read_mut: misread_names) {
                    std::string read_name = read_mut->read_name;
                    while (s_r_ptr != sample_read_list.end()) {
                        if ((*s_r_ptr)->read == read_name) {
                            read_present[s_r_ptr - sample_read_list.begin()] = mut_nuc_list_reads.find(read_mut->mut_nuc) + 1;
                            break;
                        }
                        s_r_ptr++;
                    } 
                }
                
                for (auto read: read_present) {
                    vcf_file_read_holder += std::to_string(read) + "\t";
                    if (read)
                        read_match_count[read-1] += 1;
                }
                vcf_file_read_holder += "\n";
                
                int total_count = 0;
                s_r_ptr = sample_read_list.begin();
                int start, end;
                while (s_r_ptr != sample_read_list.end()) {
                    std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                    std::smatch match;
                    if (std::regex_search((*s_r_ptr)->read, match, rgx)) {
                        start = std::stoi(match[1]);
                        end = std::stoi(match[2]);
                        if ((last_pos >= start) && (last_pos <= end))
                            total_count++;
                    }
                    s_r_ptr++;
                }
                vcf_file_read_freyja_holder += "AF=";
                for (int i = 0; i < int(read_match_count.size()); i++) {
                    if (i)
                        vcf_file_read_freyja_holder += ",";
                    float af = (float)read_match_count[i] / (float)total_count;
                    vcf_file_read_freyja_holder.append(std::to_string(af));
                }
                vcf_file_read_freyja_holder += "\n";
                
                if (map_count != (int)sample_map.size()-1) {
                    vcf_file_reads << vcf_file_read_holder;
                    vcf_file_reads_freyja << vcf_file_read_freyja_holder;
                }
            }


            //Normal VCF Writing for mutations in samples
            mut_nuc_list_reads.clear();
            ref_nuc_name_reads = {};
            misread_names.clear();
            fill(read_present.begin(), read_present.end(), 0);

            for (auto anc_smp: *(map.second)) {
                auto anc = anc_smp->ancestor_node;
                for (auto mut_anc: anc->mutations) {
                    if (mut_anc.position == map.first) {
                        ref_nuc_name_samples = MAT::get_nuc(mut_anc.ref_nuc);
                        if (!encountered_sample_ancestor) {
                            vcf_file_samples << "NC_045512v2\t" << map.first << "\t";
                            vcf_file_samples << ref_nuc_name_samples << map.first << MAT::get_nuc(mut_anc.mut_nuc);
                            encountered_sample_ancestor = true;                     
                            mut_nuc_list_samples.push_back(MAT::get_nuc(mut_anc.mut_nuc));  
                        }
                        else if (mut_nuc_list_samples.find(MAT::get_nuc(mut_anc.mut_nuc)) == std::string::npos) {
                            vcf_file_samples << "," << ref_nuc_name_samples << map.first << MAT::get_nuc(mut_anc.mut_nuc);
                            mut_nuc_list_samples.push_back(MAT::get_nuc(mut_anc.mut_nuc));                     
                        }

                        //For read vcf
                        if (read_map.find(map.first) != read_map.end()) {
                            for (auto sample_reads: *read_map[map.first]) {
                                for (auto sample: *anc_smp->sample_nodes){
                                    if (sample_reads->sample == sample) {
                                        ref_nuc_name_reads = MAT::get_nuc(mut_anc.ref_nuc);
                                        
                                        //Checking misread at this position
                                        bool mis_read_pos = false;
                                        nuc_available.clear();
                                        auto misread_ptr = misread_pos.begin();
                                        while (misread_ptr != misread_pos.end()) {
                                            if ((map.first == (*misread_ptr)->pos) && (!(*misread_ptr)->used) && (sample_reads->read == (*misread_ptr)->read) ) {
                                                std::vector<int8_t> misread_nuc;
                                                for (auto nuc: nuc_array){
                                                if (mut_anc.mut_nuc != nuc)
                                                    nuc_available.emplace_back(nuc);
                                                }
                                                std::sample(
                                                    nuc_available.begin(),
                                                    nuc_available.end(),
                                                    std::back_inserter(misread_nuc),
                                                    1,
                                                    std::mt19937{std::random_device{}()}
                                                );
                                                struct read_mut_pair *read_mut = new struct read_mut_pair;
                                                read_mut->read_name = (*misread_ptr)->read;
                                                read_mut->mut_nuc = MAT::get_nuc(misread_nuc[0]);
                                                misread_names.emplace_back(read_mut);
                                                if (!encountered_read) {
                                                    vcf_file_read_match += "NC_045512v2\t" + std::to_string(map.first) + "\t";
                                                    vcf_file_read_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(misread_nuc[0]);
                                                    vcf_file_read_freyja_match += "NC_045512v2\t" + std::to_string(map.first) + "\t";
                                                    vcf_file_read_freyja_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(misread_nuc[0]);
                                                    mut_nuc_list_reads.push_back(MAT::get_nuc(misread_nuc[0]));  
                                                    encountered_read = true;
                                                }
                                                else if (mut_nuc_list_reads.find(MAT::get_nuc(misread_nuc[0])) == std::string::npos) {
                                                    vcf_file_read_match += ",";
                                                    vcf_file_read_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(misread_nuc[0]);
                                                    vcf_file_read_freyja_match += ",";
                                                    vcf_file_read_freyja_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(misread_nuc[0]);
                                                    mut_nuc_list_reads.push_back(MAT::get_nuc(misread_nuc[0]));  
                                                }
                                                misread_nuc.clear();
                                                (*misread_ptr)->used = true;
                                                mis_read_pos = true;
                                            }
                                            misread_ptr++;
                                        }
                                        
                                        if (!mis_read_pos) {
                                            if (!encountered_read) {
                                                vcf_file_read_match += "NC_045512v2\t" + std::to_string(map.first) + "\t";
                                                vcf_file_read_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(mut_anc.mut_nuc);
                                                vcf_file_read_freyja_match += "NC_045512v2\t" + std::to_string(map.first) + "\t";
                                                vcf_file_read_freyja_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(mut_anc.mut_nuc);
                                                mut_nuc_list_reads.push_back(MAT::get_nuc(mut_anc.mut_nuc));  
                                                encountered_read = true;                     
                                            }
                                            else if (mut_nuc_list_reads.find(MAT::get_nuc(mut_anc.mut_nuc)) == std::string::npos) {
                                                vcf_file_read_match += ",";
                                                vcf_file_read_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(mut_anc.mut_nuc); 
                                                vcf_file_read_freyja_match += ",";
                                                vcf_file_read_freyja_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(mut_anc.mut_nuc); 
                                                mut_nuc_list_reads.push_back(MAT::get_nuc(mut_anc.mut_nuc));  
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }

            vcf_file_samples << "\t" <<  ref_nuc_name_samples << "\t";
            for (int i = 0; i < int(mut_nuc_list_samples.size()); i++) {
                if (i) {
                    vcf_file_samples << ",";
                }
                vcf_file_samples << mut_nuc_list_samples[i];
            }
            vcf_file_samples << "\t.\t.\t.\t.\t";
            
            read_match_count.clear();
            if (encountered_read) {
                vcf_file_reads << vcf_file_read_match;
                vcf_file_reads << "\t" <<  ref_nuc_name_reads << "\t";
                vcf_file_reads_freyja << vcf_file_read_freyja_match;
                vcf_file_reads_freyja << "\t" <<  ref_nuc_name_reads << "\t";
                for (int i = 0; i < int(mut_nuc_list_reads.size()); i++) {
                    if (i) {
                        vcf_file_reads << ",";
                        vcf_file_reads_freyja << ",";
                    }
                    vcf_file_reads << mut_nuc_list_reads[i];
                    vcf_file_reads_freyja << mut_nuc_list_reads[i];
                    read_match_count.emplace_back(0);
                }
                vcf_file_reads << "\t.\t.\t.\t.\t";
                vcf_file_reads_freyja << "\t.\t.\t";
            }

            std::vector<struct ances_sample_list*> *anc_sample_list; 
            anc_sample_list = map.second; 
            
            for (auto anc_sample: *anc_sample_list) {
                for (auto sample: *anc_sample->sample_nodes) {
                    //This doesn't work for some reason
                    //auto itr = std::find(lineage_selected.begin(), lineage_selected.end(), sample);
                    auto itr = lineage_selected.begin();
                    while (itr != lineage_selected.end()) {
                        itr = std::find(itr, lineage_selected.end(), sample);
                        if (itr != lineage_selected.end()) {
                            // If that lineage is not assigned something
                            if (!lineage_present[itr - lineage_selected.begin()]) {
                                for (auto mut_anc: anc_sample->ancestor_node->mutations) {
                                    if (mut_anc.position == map.first) {
                                        struct sample_ances_list *sample_ancestor_list = new struct sample_ances_list;
                                        sample_ancestor_list->sample = sample;
                                        std::vector<Mutation_Annotated_Tree::Node*> * ances_list = new std::vector<Mutation_Annotated_Tree::Node*>;
                                        ances_list->emplace_back(anc_sample->ancestor_node);
                                        sample_ancestor_list->ances = ances_list;
                                        sample_ancestors.emplace_back(sample_ancestor_list);
                                        
                                        if (mut_nuc_list_samples.find(MAT::get_nuc(mut_anc.mut_nuc)) != std::string::npos){
                                            //Find the correct mutation in sample
                                            lineage_present[itr - lineage_selected.begin()] = mut_nuc_list_samples.find(MAT::get_nuc(mut_anc.mut_nuc)) + 1; 
                                            
                                            //Find the same sample in read list
                                            if (read_map.find(map.first) != read_map.end()) {
                                                for (auto sample_reads: *read_map[map.first]) {
                                                    if (sample_reads->sample == sample) {
                                                        std::vector<struct sample_read_pair*>::iterator ptr;
                                                        ptr = sample_read_list.begin();
                                                        while (ptr != sample_read_list.end()) {
                                                            if ((*ptr)->read == sample_reads->read) {
                                                                bool misread_present = false;
                                                                if (misread_names.size()) {
                                                                    /////NEED to change
                                                                    for (auto read_mut: misread_names) {
                                                                        if ((*ptr)->read == read_mut->read_name) {
                                                                            read_present[ptr - sample_read_list.begin()] = mut_nuc_list_reads.find(read_mut->mut_nuc) + 1;
                                                                            misread_present = true;
                                                                            break;
                                                                        }
                                                                    }
                                                                }
                                                                if (!(misread_present))
                                                                    read_present[ptr - sample_read_list.begin()] = mut_nuc_list_reads.find(MAT::get_nuc(mut_anc.mut_nuc)) + 1;
                                                                //break;
                                                            }
                                                            ptr++;
                                                        } 
                                                    }
                                                }
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                            // Sample already has a mutation assigned in VCF
                            else {
                                for (auto mut_anc: anc_sample->ancestor_node->mutations) {
                                    if (mut_anc.position == map.first) {
                                        std::vector<struct sample_ances_list*>::iterator pointer = sample_ancestors.begin();
                                        while (pointer != sample_ancestors.end()) {
                                            if ((*pointer)->sample == sample)
                                                break;
                                            pointer++;
                                        }
                                        bool ances_found = false;
                                        for (auto anc: *(*pointer)->ances) {
                                            if (anc == anc_sample->ancestor_node) {
                                                ances_found = true;
                                                break;
                                            }
                                        }
                                        if (!ances_found) {
                                            bool is_daughter = true;
                                            for (auto anc: *(*pointer)->ances) {
                                                if (T.is_ancestor(anc_sample->ancestor_node->identifier, anc->identifier)) {
                                                is_daughter = false;
                                                break;
                                                }    
                                            }
                                            if (is_daughter) {
                                                if (mut_nuc_list_samples.find(MAT::get_nuc(mut_anc.mut_nuc)) != std::string::npos){
                                                    //Find the correct mutation in sample
                                                    lineage_present[itr - lineage_selected.begin()] = mut_nuc_list_samples.find(MAT::get_nuc(mut_anc.mut_nuc)) + 1; 

                                                    //Find the same sample in read list
                                                    if (read_map.find(map.first) != read_map.end()) {
                                                        for (auto sample_reads: *read_map[map.first]) {
                                                            if (sample_reads->sample == sample) {
                                                                std::vector<struct sample_read_pair*>::iterator ptr;
                                                                ptr = sample_read_list.begin();
                                                                while (ptr != sample_read_list.end()) {
                                                                    if ((*ptr)->read == sample_reads->read) {
                                                                        bool misread_present = false;
                                                                        if (misread_names.size()) {
                                                                        ///////NEED to change
                                                                            for (auto read_mut: misread_names) {
                                                                                if ((*ptr)->read == read_mut->read_name) {
                                                                                    read_present[ptr - sample_read_list.begin()] = mut_nuc_list_reads.find(read_mut->mut_nuc) + 1;
                                                                                    misread_present = true;
                                                                                    break;
                                                                                }
                                                                            }
                                                                        }
                                                                        if (!(misread_present))
                                                                            read_present[ptr - sample_read_list.begin()] = mut_nuc_list_reads.find(MAT::get_nuc(mut_anc.mut_nuc)) + 1;
                                                                        //break;
                                                                    }
                                                                    ptr++;
                                                                } 
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            (*pointer)->ances->emplace_back(anc_sample->ancestor_node);
                                        }
                                    }
                                }
                            }
                            itr++;
                        }
                    }
                }
            }
            
            //Print in VCF from lineage_present
            for (auto hap: lineage_present)
                vcf_file_samples << hap <<"\t";
            vcf_file_samples << "\n";
            
            if (encountered_read) {
                for (auto read: read_present) {
                    vcf_file_reads << read <<"\t";
                    if (read)
                        read_match_count[read-1] += 1;
                }
                vcf_file_reads << "\n";
            }

            int total_count = 0;
            int start, end;
            using my_mutex_t = tbb::queuing_mutex;
            my_mutex_t my_mutex;
            static tbb::affinity_partitioner ap;
            tbb::parallel_for( tbb::blocked_range<size_t>(0, sample_read_list.size()),
                [&](tbb::blocked_range<size_t> k) {
                    for (size_t s = k.begin(); s < k.end(); ++s) {
                        auto sample_read = sample_read_list[s];
                        std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                        std::smatch match;
                        my_mutex_t::scoped_lock my_lock;
                        my_lock.acquire(my_mutex);
                        if (std::regex_search(sample_read->read, match, rgx)) {
                            start = std::stoi(match[1]);
                            end = std::stoi(match[2]);
                            if ((map.first >= start) && (map.first <= end)) {
                                    total_count++;
                            }
                        }
                        my_lock.release();
                    }
                },
            ap);
            for (int i = 0; i < int(read_match_count.size()); i++) {
                if (!i)
                    vcf_file_reads_freyja << "AF=";
                if (i)
                    vcf_file_reads_freyja << ",";
                float af = (float)read_match_count[i] / (float)total_count;
                vcf_file_reads_freyja << std::to_string(af);
                if (i == int(read_match_count.size() - 1))
                    vcf_file_reads_freyja << "\n";
            }

            while (allele_pos <= map.first) {
                int read_count = 0;
                int start, end;
                using my_mutex_t = tbb::queuing_mutex;
                my_mutex_t my_mutex;
                static tbb::affinity_partitioner ap;
                tbb::parallel_for( tbb::blocked_range<size_t>(0, sample_read_list.size()),
                    [&](tbb::blocked_range<size_t> k) {
                        for (size_t s = k.begin(); s < k.end(); ++s) {
                            auto sample_read = sample_read_list[s];
                            std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                            std::smatch match;
                            my_mutex_t::scoped_lock my_lock;
                            my_lock.acquire(my_mutex);
                            if (std::regex_search(sample_read->read, match, rgx)) {
                                start = std::stoi(match[1]);
                                end = std::stoi(match[2]);
                                if ((allele_pos >= start) && (allele_pos <= end)) {
                                        read_count++;
                                }
                            }
                            my_lock.release();
                        }
                    },
                ap);
                depth_file_read_freyja_holder += "NC_045512v2\t" + std::to_string(allele_pos+1) + "\t" + ref_seq[allele_pos-1] + "\t" + std::to_string(read_count) + "\n";
                allele_pos++;
            }


            if ((map_count == ((int)sample_map.size()-1)) && (misreads_eligible.size())) {
                vcf_file_reads << vcf_file_read_holder;
                vcf_file_reads_freyja << vcf_file_read_freyja_holder;
            }
            if (map_count == ((int)sample_map.size()-1)) {
                while (allele_pos < (int)ref_seq.size()) {
                    int read_count = 0;
                    int start, end;
                    using my_mutex_t = tbb::queuing_mutex;
                    my_mutex_t my_mutex;
                    static tbb::affinity_partitioner ap;
                    tbb::parallel_for( tbb::blocked_range<size_t>(0, sample_read_list.size()),
                        [&](tbb::blocked_range<size_t> k) {
                            for (size_t s = k.begin(); s < k.end(); ++s) {
                                auto sample_read = sample_read_list[s];
                                std::regex rgx(".*READ_(\\w+)_(\\w+).*");
                                std::smatch match;
                                my_mutex_t::scoped_lock my_lock;
                                my_lock.acquire(my_mutex);
                                if (std::regex_search(sample_read->read, match, rgx)) {
                                    start = std::stoi(match[1]);
                                    end = std::stoi(match[2]);
                                    if ((allele_pos >= start) && (allele_pos <= end)) {
                                            read_count++;
                                    }
                                }
                                my_lock.release();
                            }
                        },
                    ap);
                    depth_file_read_freyja_holder += "NC_045512v2\t" + std::to_string(allele_pos+1) + "\t" + ref_seq[allele_pos-1] + "\t" + std::to_string(read_count) + "\n";
                    allele_pos++;
                }
                depth_file_reads_freyja << depth_file_read_freyja_holder;
            }
            previous_pos = map.first;
            map_count ++;
        }

        fprintf(stderr, "VCFs Written in %ld min\n\n", (timer.Stop() / (60 * 1000)));

        boost::iostreams::close(outbuf_samples);
        boost::iostreams::close(outbuf_reads);
        boost::iostreams::close(outbuf_reads_freyja);
        boost::iostreams::close(outbuf_reads_freyja_depth);
        outfile_samples.close();
        outfile_reads.close();
        outfile_reads_freyja.close();
        outfile_reads_freyja_depth.close();
    }
   

    std::unordered_map<size_t, struct read_info*> read_map;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    std::vector<std::string> vcf_samples;
    //read samples.vcf to check how close are peaks to samples 
    read_sample_vcf(vcf_samples, vcf_filename_samples);
    //Get the reads.vcf data
    read_vcf(read_map, vcf_filename_reads);
    //Core algorithm 
    analyze_reads(T, T_ref, read_map, node_score, vcf_samples, barcode_file, read_abundance_vcf);
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
    std::vector<size_t> remaining_reads(read_map.size());
    
    //Depth first expansion to get all nodes in the tree and 
    dfs = T.depth_first_expansion(); 
    
    //GREEDY ALGORITHM for peak nodes
    //    1. Place remaining reads on tree
    //    2. Remove nodes not to be considered (present in Prohibited nodes)
    //    3. Find top scoring nodes
    //    4. Only consider highest scoring node(s) not seen before
    //    5. Remove reads mapped to these nodes
    //    6. Don't consider top scoring nodes in the neighborhood of other top nodes
    //    7. Add neighbors of selected top scoring nodes to Prohibited nodes 
    
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

        //SORT node_scores
        std::vector<std::pair<MAT::Node*, double>> node_score_vector, top_n_node_score;
        node_score_vector.reserve(node_score.size());
        // Convert the concurrent hash map into a vector in parallel
        tbb::parallel_for_each(node_score.begin(), node_score.end(), [&node_score_vector](const auto& element) {
            node_score_vector.emplace_back(element);
        });
        node_score.clear();
        // Sort the node_score_vector filled above
        tbb::parallel_sort(node_score_vector.begin(), node_score_vector.end(), 
             [&T](const auto& a, const auto& b) {
                return compare_node_score(T, a, b);
        });

        // Take top_n values from nod_score_vector into top_n_node_score
        int node_count = 0;
        auto top_score = node_score_vector.begin()->second;
        for (const auto& n_s: node_score_vector) {
            if (abs(top_score - n_s.second) < 1e-9) {
                top_n_node_score.emplace_back(n_s);
                node_count++;
                if (node_count == top_n)
                    break;
            }
            else 
                break;
        }
        
        //Only consider peak nodes not seen before
        std::vector<bool> peak_vec(top_n_node_score.size(), true);
        //Find the reads not mapped to the best nodes
        auto top_n_itr = top_n_node_score.begin();
        while (top_n_itr != top_n_node_score.end()) {
            auto curr_node = top_n_itr->first;
            //Check if current node is not in neighborhood of peaks seen till now
            if (peak_vec[top_n_itr - top_n_node_score.begin()]) {
                curr_peak_nodes.emplace_back(curr_node);
                auto curr_clade = get_clade(T, curr_node);
                printf("PEAK: %s, Score: %f, Clade:%s\n",curr_node->identifier.c_str(), top_n_itr->second, curr_clade.c_str());
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
                //Erase from remove_reads
                std::sort(remove_reads.begin(), remove_reads.end());
                auto rmr_itr = remaining_reads.begin();
                auto rr_itr = remove_reads.begin();
                while ((rmr_itr != remaining_reads.end()) && (rr_itr != remove_reads.end())) {
                    if (*rmr_itr == *rr_itr) {
                        rr_itr++;
                        rmr_itr = remaining_reads.erase(rmr_itr);
                    }
                    else
                        rmr_itr++;
                }
                remove_reads.clear();
            }
            //If present in neighbourhood of current peaks, move to next node
            else {
                top_n_itr++;
                continue;
            }
            
            //Remove nearby peaks from further analysis
            auto top_n_peak_cmp_itr = top_n_itr;
            std::advance(top_n_peak_cmp_itr, 1);
            while (top_n_peak_cmp_itr != top_n_node_score.end()) {
                auto cmp_node = top_n_peak_cmp_itr->first;
                //Check next peak if this peak is already in neighborhood of some other peak
                if (!peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()]) {
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
        //  1. Find the farthest ancestor of every peak node within m_dist_thresh
        //  2. Recursively only analyze its children within m_dist_thresh
        //  3. Only include unique neighborhood nodes to prohibited_nodes
        using my_mutex_t = tbb::queuing_mutex;
        my_mutex_t my_mutex;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, curr_peak_nodes.size()),
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
                            my_mutex_t::scoped_lock my_lock;
                            my_lock.acquire(my_mutex);
                            if (!std::binary_search(prohibited_nodes.begin(), prohibited_nodes.end(), present_node))
                                prohibited_nodes.emplace_back(present_node);
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

    add_neighbor_peaks(T, peak_nodes, neighbor_dist_thresh, neighbor_peaks_thresh);
    generate_regression_abundance_data(T, peak_nodes, read_map, barcode_file, read_mutation_depth_vcf);
}

//Add neighboring nodes to peaks
void add_neighbor_peaks(const MAT::Tree &T, std::vector<MAT::Node*> &peak_nodes, const int neighbor_dist_thresh, const int neighbor_peaks_thresh) {
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
                
                //Keep looking for neighbors until neighbors added == neighbor_dist_thresh AND neighbors keep getting added 
                while ((neighbor_added < neighbor_peaks_thresh) && (feasible_branch_length == branch_length_limit)) {
                    // 1. Search children nodes within branch_length_limit
                    branch_length_limit++;
                    auto curr_node = peak;
                    int curr_branch_length = 0;
                    feasible_branch_length = child_nodes_addition(T, my_mutex, peak_nodes, neighbor_nodes, neighbor_added, prev_visited_nodes, peak, curr_node, curr_branch_length, branch_length_limit, neighbor_peaks_thresh, neighbor_dist_thresh);
                    
                    // 2. Search parent nodes within branch_length_limit
                    while ((curr_branch_length < branch_length_limit) && (neighbor_added < neighbor_peaks_thresh)) {
                        if (curr_node->parent == NULL)
                            break;
                        curr_node = curr_node->parent;
                        curr_branch_length++;
                        int temp_feasible_branch_length = child_nodes_addition(T, my_mutex, peak_nodes, neighbor_nodes, neighbor_added, prev_visited_nodes, peak, curr_node, curr_branch_length, branch_length_limit, neighbor_peaks_thresh, neighbor_dist_thresh);
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
int child_nodes_addition(const MAT::Tree &T, my_mutex_t& my_mutex, const std::vector<MAT::Node*>& peak_nodes, std::vector<MAT::Node*>& neighbor_nodes, int& neighbor_added, std::vector<MAT::Node*>& prev_visited_nodes, MAT::Node* peak, MAT::Node* init_node, const int& init_branch_length, const int& branch_length_thresh, const int& neighbor_peaks_thresh, const int& neighbor_dist_thresh) {
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
            int m_dist = mutation_distance(T, T, curr_node_branch.node, peak);
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
        
    fprintf(stderr,"Barcode and VCF file writing took %ld sec\n\n", (timer.Stop() / 1000));
}

//Comparing mutations for sorting a vector 
bool compare_mutations(const MAT::Mutation &a, const MAT::Mutation &b) {
    return a.position < b.position;
}

//Comparing different node_scores  
bool compare_node_score (const MAT::Tree &T, const std::pair<MAT::Node*, double>& a, const std::pair<MAT::Node*, double>& b) {
    //Compare based on double values
    if (abs(a.second - b.second) > 1e-9) {
        return a.second > b.second;
    }
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
        for (auto c: curr_node->children) {
            remaining_nodes.push(c);
        }
    }
    return leaves.size();
}