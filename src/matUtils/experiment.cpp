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
    bool old_vcf = false;
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

    
    fprintf(stderr,"\nRef Seq Length: %ld\n\n", ref_seq.size());
    fprintf(stderr,"Leaves in tree: %ld\n\n", T.get_num_leaves());
    //std::cout << "\nRead Length: " << read_length << ", Range: " << 100 - ((read_length / 2) + (read_length % 2)) << "-" 
    //<< 100 + (read_length / 2) - 1 << std::endl;
    //std::cout << "\n";



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
                        nuc_used.emplace_back(MAT::get_nuc_id(ref_seq[(*misread_ptr)->pos]));
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
                        ref_nuc_name_reads = ref_seq[(*misread_ptr)->pos];
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
            depth_file_read_freyja_holder += "NC_045512v2\t" + std::to_string(allele_pos+1) + "\t" + ref_seq[allele_pos] + "\t" + std::to_string(read_count) + "\n";
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
                depth_file_read_freyja_holder += "NC_045512v2\t" + std::to_string(allele_pos+1) + "\t" + ref_seq[allele_pos] + "\t" + std::to_string(read_count) + "\n";
                allele_pos++;
            }
            depth_file_reads_freyja << depth_file_read_freyja_holder;
        }
        previous_pos = map.first;
        map_count ++;
    }

    fprintf(stderr, "VCFs Written in %ld sec \n\n", (timer.Stop() / 1000));

    boost::iostreams::close(outbuf_samples);
    boost::iostreams::close(outbuf_reads);
    boost::iostreams::close(outbuf_reads_freyja);
    boost::iostreams::close(outbuf_reads_freyja_depth);
    outfile_samples.close();
    outfile_reads.close();
    outfile_reads_freyja.close();
    outfile_reads_freyja_depth.close();
    }

    std::vector<struct read_info *> read_ids;
    
    timer.Start();
    read_vcf(vcf_filename_reads, read_ids);
    fprintf(stderr,"Reads_VCF parsed in %ld msec\n\n", timer.Stop());

    tbb::concurrent_hash_map<MAT::Node*, score_read> node_score;
    tbb::concurrent_hash_map<size_t, struct min_parsimony> read_min_parsimony;
    place_reads(T, dfs, read_ids, node_score, read_min_parsimony, vcf_filename_reads);
    analyze_reads(T, dfs, read_ids, node_score, read_min_parsimony);
}


void place_reads(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, const std::vector<struct read_info*> &read_ids, tbb::concurrent_hash_map<MAT::Node*, score_read> &node_score, tbb::concurrent_hash_map<size_t, struct min_parsimony> &read_min_parsimony, std::string vcf_filename_reads) {
    fprintf(stderr, "Total nodes: %ld, Reads: %ld\n\n", dfs.size(), read_ids.size());
    timer.Start();

    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_ids.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t r = k.begin(); r < k.end(); ++r) {
                auto rp = read_ids[r];
                std::stack<struct parsimony> parsimony_stack;
                struct min_parsimony min_par;
                while (!parsimony_stack.empty())
                    parsimony_stack.pop();
                min_par.idx_list.clear();
                min_par.par_list.clear();
                min_par.is_sibling_list.clear();
            
                for (size_t i = 0; i < dfs.size(); i++) {
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

                    //Checking min_parsimony
                    int new_min_par = -1; 
                    if (!(min_par.par_list.size()))
                        new_min_par = 1;
                    else if (curr_node_par_mut.size() < min_par.par_list[0].size())
                        new_min_par = 1;
                    else if (curr_node_par_mut.size() == min_par.par_list[0].size())
                        new_min_par = 0;
                    
                    if (new_min_par == 1) {
                        min_par.idx_list.clear();
                        min_par.par_list.clear();
                        min_par.is_sibling_list.clear();
                        min_par.idx_list.emplace_back(i);
                        min_par.par_list.emplace_back(curr_node_par_mut);
                        if (uniq_curr_node_mut.size())
                            min_par.is_sibling_list.emplace_back(true);
                        else
                            min_par.is_sibling_list.emplace_back(false);
                    }
                    else if (new_min_par == 0) {
                        min_par.idx_list.emplace_back(i);
                        min_par.par_list.emplace_back(curr_node_par_mut);
                        if (uniq_curr_node_mut.size())
                            min_par.is_sibling_list.emplace_back(true);
                        else
                            min_par.is_sibling_list.emplace_back(false);
                    }

                    //Updating parsimony to be stored a a child
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

                //Keeping tab on weighted read score being mapped to each parsimonious node
                for (auto n_idx: min_par.idx_list) {
                    auto node = dfs[n_idx];
                    //Give score to parsimonious node
                    struct score_read sc;
                    sc.score = (1.0 / pow(min_par.idx_list.size(), 2.0));
                    sc.reads.emplace_back(r);
                    tbb::concurrent_hash_map<MAT::Node*, score_read>::accessor ac;
                    auto created = node_score.insert(ac, std::make_pair(node, sc));
                    if (!created) {
                        ac->second.score += (1.0 / pow(min_par.idx_list.size(), 2.0));
                        ac->second.reads.emplace_back(r);
                    }
                    ac.release();
                }

                //Only keeping for the check
                tbb::concurrent_hash_map<size_t, struct min_parsimony>::accessor ac;
                read_min_parsimony.insert(ac, r);
                ac->second = min_par;
                ac.release();
            }
        },
        ap);

    fprintf(stderr,"Reads placed in %ld sec\n\n", (timer.Stop() / 1000));
    
}


void analyze_reads(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, const std::vector<struct read_info*> &read_ids, tbb::concurrent_hash_map<MAT::Node*, score_read> &node_score, tbb::concurrent_hash_map<size_t, struct min_parsimony> &read_min_parsimony) {
    timer.Start();
    //Check to ensure every read has its corresponding sample as its most parsimonious position
    unsigned long long avg = 0;
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
    static tbb::affinity_partitioner ap;
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_ids.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t r = k.begin(); r < k.end(); ++r) {
                auto rp = read_ids[r];
                tbb::concurrent_hash_map<size_t, struct min_parsimony>::const_accessor k_ac;
                read_min_parsimony.find(k_ac, r);
                auto min_par = k_ac->second;
                k_ac.release();
                std::string target = rp->read;
                size_t pos = target.find("_READ");
                target.erase(pos);
                bool found = false;
                size_t idx, i;
                for (i = 0; i < min_par.idx_list.size(); i++) {
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
                    
                    fprintf(stderr, "Read: %s, read mutations: %ld, Parsimony score = %ld, parsimonious positions: %ld, sample_idx = %ld \n\n", rp->read.c_str(), rp->mutations.size(), min_par.par_list[0].size(), min_par.par_list.size(), idx);
                    std::cout << "\n";

                }

                my_mutex_t::scoped_lock my_lock{my_mutex};
                    avg += (min_par.idx_list.size() / read_ids.size());
            }
        },
            ap);

    fprintf(stderr, "Avg par pos: %lld found in %ld sec\n", avg, (timer.Stop() / 1000));

   //GREEDY ALGORITHM for getting handful of clades
   //Get top scoring nodes
    timer.Start();
    int top_n = 25;
    size_t b_dist_thresh = 100;
    double score_thresh = 0.8;
    std::vector<std::pair<MAT::Node*, score_read>> top_n_node_score(top_n);
    std::vector<size_t> remaining_reads;
    std::vector<std::string> selected_clades;
    
    auto r_min_itr = read_min_parsimony.begin();
    while (r_min_itr != read_min_parsimony.end()) {
        remaining_reads.emplace_back(r_min_itr->first);
        r_min_itr++;
    }
    
    while (!remaining_reads.empty()) {
        //Sorting the node scores
        std::partial_sort_copy(node_score.begin(),
                            node_score.end(),
                            top_n_node_score.begin(),
                            top_n_node_score.end(),
                            [](std::pair<const MAT::Node*, score_read> const& l,
                               std::pair<const MAT::Node*, score_read> const& r)
                            {
                                return l.second.score > r.second.score;
                            });
        node_score.clear();

        //Find the reads not mapped to the best nodes
        auto top_n_itr = top_n_node_score.begin();
        while (top_n_itr != top_n_node_score.end()) {
            auto curr_clade = get_clade(T, top_n_itr->first);
            auto c_check_itr = std::find(selected_clades.begin(), selected_clades.end(), curr_clade);
            //Allow top scorer or others with equal score but different clade
            if ( (abs(top_n_itr->second.score - top_n_node_score[0].second.score) < 1e-9) && (c_check_itr == selected_clades.end()) ) {
                //Store the clades corresponding to top scoring node
                selected_clades.emplace_back(curr_clade);
                //Remove mapped reads from remaining_reads
                auto rem_r_itr = remaining_reads.begin();
                while (rem_r_itr != remaining_reads.end()) {
                    auto itr = std::find(top_n_itr->second.reads.begin(), top_n_itr->second.reads.end(), *rem_r_itr);
                    if (itr != top_n_itr->second.reads.end())
                        rem_r_itr = remaining_reads.erase(rem_r_itr);
                    else rem_r_itr++;
                }
            }
            else if (abs(top_n_itr->second.score - top_n_node_score[0].second.score) > 1e-9)
                break;
            top_n_itr++;
        }        
        top_n_node_score.clear();
        top_n_node_score.resize(top_n);

        //Calculating node score for remaining reads
        tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto r = remaining_reads[i];
                    tbb::concurrent_hash_map<size_t, struct min_parsimony>::const_accessor k_ac;
                    read_min_parsimony.find(k_ac, r);
                    auto min_par = k_ac->second;
                    k_ac.release();
                    for (auto n_idx: min_par.idx_list) {
                        auto node = dfs[n_idx];
                        //Give score to parsimonious node
                        struct score_read sc;
                        sc.score = (1.0 / pow(min_par.idx_list.size(), 2.0));
                        sc.reads.emplace_back(r);
                        tbb::concurrent_hash_map<MAT::Node*, score_read>::accessor ac;
                        auto created = node_score.insert(ac, std::make_pair(node, sc));
                        if (!created) {
                            ac->second.score += (1.0 / pow(min_par.idx_list.size(), 2.0));
                            ac->second.reads.emplace_back(r);
                        }
                        ac.release();
                    }
                }
            },
        ap);
    }


    //Print the clades otained from greedy algorithm
    for (size_t i = 0; i < selected_clades.size(); i++) {
        printf("CLADE Selected: %s\n", selected_clades[i].c_str());
    }

    //Rescoring the nodes only belonging to selected lineages
    node_score.clear();        
    tbb::parallel_for( tbb::blocked_range<size_t>(0, read_ids.size()),
        [&](tbb::blocked_range<size_t> k) {
            for (size_t r = k.begin(); r < k.end(); ++r) {
                tbb::concurrent_hash_map<size_t, struct min_parsimony>::const_accessor k_ac;
                read_min_parsimony.find(k_ac, r);
                auto min_par = k_ac->second;
                k_ac.release();
                for (auto n_idx: min_par.idx_list) {
                    auto node = dfs[n_idx];
                    auto clade = get_clade(T, node);
                    auto clade_itr = std::find(selected_clades.begin(), selected_clades.end(), clade);
                    //On finding the node among the selected clades, update its score
                    if (clade_itr != selected_clades.end()) {
                        struct score_read sc;
                        sc.score = (1.0 / pow(min_par.idx_list.size(), 2.0));
                        sc.reads.emplace_back(r);
                        tbb::concurrent_hash_map<MAT::Node*, score_read>::accessor ac;
                        auto created = node_score.insert(ac, std::make_pair(node, sc));
                        if (!created) {
                            ac->second.score += (1.0 / pow(min_par.idx_list.size(), 2.0));
                            ac->second.reads.emplace_back(r);
                        }
                        ac.release();
                    }

                }
            }
        },
     ap);
    
    std::cout << "\n Top Nodes after rescoring among limited clades\n\n";

    //Greedy Peak Detection
    std::vector<std::pair<MAT::Node*, score_read>> peak_node_score;
    r_min_itr = read_min_parsimony.begin();
    while (r_min_itr != read_min_parsimony.end()) {
        remaining_reads.emplace_back(r_min_itr->first);
        r_min_itr++;
    }
    while (!remaining_reads.empty()) {
        //Sorting the node scores
        top_n_node_score.clear();
        top_n_node_score.resize(top_n);
        std::partial_sort_copy(node_score.begin(),
                            node_score.end(),
                            top_n_node_score.begin(),
                            top_n_node_score.end(),
                            [](std::pair<const MAT::Node*, score_read> const& l,
                               std::pair<const MAT::Node*, score_read> const& r)
                            {
                                return l.second.score > r.second.score;
                            });
        node_score.clear();

        //Get a leaf top node not seen before 
        std::vector<bool> peak_vec(top_n_node_score.size(), true);
        auto top_n_itr = top_n_node_score.begin();
        while (top_n_itr != top_n_node_score.end()) {
            //Find unseen top node
            auto p_itr = peak_node_score.begin();
            while (p_itr != peak_node_score.end()) {
                if (p_itr->first == top_n_itr->first)
                    break;
                p_itr++;
            }
            //Check if leaf
            if ((top_n_itr->first->is_leaf()) && (p_itr == peak_node_score.end())) 
                break;
            // Skip for peak consideration if node not leaf or seen in pervious iteration
            peak_vec[top_n_itr - top_n_node_score.begin()] = false;
            top_n_itr++;
        }
        double top_score = top_n_itr->second.score;
        
        //Find the peak nodes and remove reads mapped to peak nodes
        while (top_n_itr != top_n_node_score.end()) {
            //Peak Finding
            auto curr_node = top_n_itr->first;
            if (!(curr_node->is_leaf()))
                peak_vec[top_n_itr - top_n_node_score.begin()] = false;
            if (!peak_vec[top_n_itr - top_n_node_score.begin()]) {
                top_n_itr++;
                continue;
            }
            auto top_n_peak_cmp_itr = top_n_itr + 1;
            while (top_n_peak_cmp_itr != top_n_node_score.end()) {
                auto cmp_node = top_n_peak_cmp_itr->first;
                if (!(cmp_node->is_leaf()))
                    peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = false;
                if (!peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()]) {
                    top_n_peak_cmp_itr++;
                    continue;
                }
                //Calculate branch distance if scores of both nodes are not equal
                if (abs(top_n_itr->second.score - top_n_peak_cmp_itr->second.score) < 1e-9 )
                    peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = true;
                else {
                    size_t b_dist = branch_distance(curr_node, cmp_node);
                    if (b_dist > b_dist_thresh)
                        peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = true;
                    else
                        peak_vec[top_n_peak_cmp_itr - top_n_node_score.begin()] = false;
                }
                top_n_peak_cmp_itr++;
            }
            //Only add unique leaf nodes in peak score with score >= 0.1 top_score
            auto p_itr = peak_node_score.begin();
            while (p_itr != peak_node_score.end()) {
                if (p_itr->first == curr_node)
                    break;
                p_itr++;
            }
            if ((p_itr == peak_node_score.end()) && (top_n_itr->second.score >= (score_thresh * top_score))) {
                peak_node_score.emplace_back(*top_n_itr);
                //Remove reads mapped to this node from remaining_reads
                auto rem_r_itr = remaining_reads.begin();
                while (rem_r_itr != remaining_reads.end()) {
                    auto itr = std::find(top_n_itr->second.reads.begin(), top_n_itr->second.reads.end(), *rem_r_itr);
                    if (itr != top_n_itr->second.reads.end())
                        rem_r_itr = remaining_reads.erase(rem_r_itr);
                    else rem_r_itr++;
                }
            }
            top_n_itr++;
        }        

        //Calculating node score for remaining reads
        tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    auto r = remaining_reads[i];
                    tbb::concurrent_hash_map<size_t, struct min_parsimony>::const_accessor k_ac;
                    read_min_parsimony.find(k_ac, r);
                    auto min_par = k_ac->second;
                    k_ac.release();
                    for (auto n_idx: min_par.idx_list) {
                        auto node = dfs[n_idx];
                        //Give score to parsimonious node
                        struct score_read sc;
                        sc.score = (1.0 / pow(min_par.idx_list.size(), 2.0));
                        sc.reads.emplace_back(r);
                        tbb::concurrent_hash_map<MAT::Node*, score_read>::accessor ac;
                        auto created = node_score.insert(ac, std::make_pair(node, sc));
                        if (!created) {
                            ac->second.score += (1.0 / pow(min_par.idx_list.size(), 2.0));
                            ac->second.reads.emplace_back(r);
                        }
                        ac.release();
                    }
                }
            },
        ap);

    }

    for (auto n_s: peak_node_score) {
        auto clade = get_clade(T, n_s.first);
        printf("PEAK score = %f, read_count = %lu, Node: %s, Clade: %s\n", n_s.second.score, n_s.second.reads.size(), n_s.first->identifier.c_str(), clade.c_str());
    }


    fprintf(stderr,"\nAnalysis algorithms took %ld sec\n\n", (timer.Stop() / 1000));
}


void read_vcf(std::string vcf_filename_reads, std::vector<struct read_info*> &read_ids) {
    // Boost library used to stream the contents of the input VCF file in
    // uncompressed or compressed .gz format
    std::ifstream infile(vcf_filename_reads, std::ios_base::in | std::ios_base::binary);
    if (!infile) {
        fprintf(stderr, "ERROR: Could not open the VCF file: %s!\n", vcf_filename_reads.c_str());
        exit(1);
    }
    boost::iostreams::filtering_istream instream;
    try {
        if (vcf_filename_reads.find(".gz\0") != std::string::npos) {
            instream.push(boost::iostreams::gzip_decompressor());
        }
        instream.push(infile);
    } catch(const boost::iostreams::gzip_error& e) {
        std::cout << e.what() << '\n';
    }
    bool header_found = false;
    std::vector<size_t> missing_idx;
    std::string s;
    while (instream.peek() != EOF) {
        std::getline(instream, s);
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if ((not header_found) && (words.size() > 1)) {
            if (words[1] == "POS") {
                for (size_t j=9; j < words.size(); j++) {
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
                header_found = true;
                //fprintf(stderr, "Reads: %ld\n", read_ids.size());
            }
        } else if (header_found) {
            if (words.size() != 9+read_ids.size()) {
                fprintf(stderr, "ERROR! Incorrect VCF format. Expected %zu columns but got %zu.\n", 9+read_ids.size(), words.size());
                exit(1);
            }
            std::vector<std::string> alleles;
            alleles.clear();
            MAT::string_split(words[4], ',', alleles);
            for (size_t k = 0; k < missing_idx.size(); k++) {
                size_t j = missing_idx[k];
                auto iter = read_ids.begin();
                std::advance(iter, k);
                if (iter != read_ids.end()) {
                    MAT::Mutation m;
                    m.chrom = words[0];
                    m.position = std::stoi(words[1]);
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
                    //if ((m.mut_nuc & (m.mut_nuc-1)) !=0) {
                    //    (*iter)num_ambiguous++;
                    //}
                }
            }
        }
    }
}


//Function to calculation distance between two nodes
size_t branch_distance(MAT::Node* N1, MAT::Node* N2) {
    if (N1 == N2)
        return 0;
    size_t distance = 0;
    MAT::Node node1 = *N1, node2 = *N2;
    while (node1.level > node2.level) {
        node1 = *node1.parent;
        distance++;
    } 
    while (node1.level < node2.level) {
        node2 = *node2.parent;
        distance++; 
    }
    while (node1.parent != node2.parent) {
        node1 = *node1.parent;
        node2 = *node2.parent;
        distance += 2;
    }
    distance += 2;
    return distance;
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