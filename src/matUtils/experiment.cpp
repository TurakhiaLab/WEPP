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
    bool old_vcf = true;
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
    std::vector<MAT::Node*> traversal, lineage_list, lineage_selected;
    
    //Depth first expansion to get all nodes in the tree and 
    // comparison with given lineage to get all nodes of the required lineage 
    traversal = T.depth_first_expansion(T.root); 
    
    if (!old_vcf) {
    for (auto lineage: in_lineage) {
        bool list_start = false; 
        for (auto n: traversal) {
            for (auto str: n->clade_annotations) {
                if (str == lineage) {
                    std::queue<Mutation_Annotated_Tree::Node*> remaining_nodes;
                    remaining_nodes.push(n);
                    while (remaining_nodes.size() > 0) {
                        Mutation_Annotated_Tree::Node* curr_node = remaining_nodes.front();
                        remaining_nodes.pop();
                        for (auto curr_mut: curr_node->clade_annotations) {
                            if ((curr_mut == "") || curr_mut == lineage) {
                                if (curr_node->children.size() == 0)
                                    lineage_list.emplace_back(curr_node);
                                else {
                                    for (auto c: curr_node->children) {
                                        remaining_nodes.push(c);
                                    }
                                }
                                break;
                            }
                        }
                    }
                    list_start = true;
                    break;
                }
            }
            if (list_start)
                break;
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
            printf("Sample Name: %s\n", (*lineage_ptr)[rand_val]->identifier.c_str());
        }
        lineage_ptr++;
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
                if (back_mut_map.find(mut.position) == back_mut_map.end()) 
                    if (mut.ref_nuc != mut.mut_nuc)
                        back_mutation =  false;
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
                            if (!(T.is_ancestor(node->identifier, ances->ancestor_node->identifier))) {
                                for (auto leaf: *ances->sample_nodes) {
                                    if (T.is_ancestor(node->identifier, leaf->identifier) || (node->identifier == leaf->identifier)) {
                                        ances->sample_nodes->erase(std::remove(ances->sample_nodes->begin(), ances->sample_nodes->end(), leaf), ances->sample_nodes->end());
                                    }
                                }
                            }
                            // New Back Mutation node is ancestor of previous node in map, i.e., BM is irrelevant
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
                        for (auto ances: *back_mut_list) {
                            // If new Node is ancestor of Back Mutation Node and sample does not belong to BM 
                            if ((T.is_ancestor(node->identifier, ances->identifier)) && (!( (ances->identifier == anc[0]->identifier) || (T.is_ancestor(ances->identifier, anc[0]->identifier) )))) {
                                no_BM += 1;
                            }
                            // If Back Mutation node is ancestor of new node
                            else if (T.is_ancestor(ances->identifier, node->identifier)) {
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
                    // New node with Back Mutation
                    else {
                        std::vector<struct ances_sample_list*> *anc_sample_list; 
                        if (sample_map.find(mut.position) != sample_map.end()) {
                            anc_sample_list = sample_map[mut.position];
                            for (auto ances: *anc_sample_list) {
                                // New Back mutation is not ancestor of node from sample_map
                                if (!(T.is_ancestor(node->identifier, ances->ancestor_node->identifier))) {
                                    for (auto leaf: *ances->sample_nodes) {
                                        if (T.is_ancestor(node->identifier, leaf->identifier) || (node->identifier == leaf->identifier))
                                            ances->sample_nodes->erase(std::remove(ances->sample_nodes->begin(), ances->sample_nodes->end(), leaf), ances->sample_nodes->end());
                                    }
                                    if (ances->sample_nodes->size() < 1)
                                        anc_sample_list->erase(std::remove(anc_sample_list->begin(), anc_sample_list->end(), ances), anc_sample_list->end());
                                }
                                // New Back Mutation is ancestor of node in map i.e., BM is irrelevant
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
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_samples, outbuf_reads;
    if (vcf_filename_samples.find(".gz\0") != std::string::npos) {
        outbuf_samples.push(boost::iostreams::gzip_compressor());
    }
    if (vcf_filename_reads.find(".gz\0") != std::string::npos) {
        outbuf_reads.push(boost::iostreams::gzip_compressor());
    }
    outbuf_samples.push(outfile_samples);
    outbuf_reads.push(outfile_reads);
    std::ostream vcf_file_samples(&outbuf_samples);
    std::ostream vcf_file_reads(&outbuf_reads);

    vcf_file_samples << "##fileformat=VCFv4.2\n";
    vcf_file_samples << "##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n";
    vcf_file_samples << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    vcf_file_reads << "##fileformat=VCFv4.2\n";
    vcf_file_reads << "##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n";
    vcf_file_reads << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    
    for (auto node: lineage_selected) {
        vcf_file_samples << "\t" << node->identifier;
    }

    for (auto sample: sample_read_list) {
        vcf_file_reads << "\t" << sample->read;
    }
    
    vcf_file_samples << "\n";
    vcf_file_reads << "\n";
    
    std::vector<int> lineage_present(lineage_selected.size());
    std::vector<int> read_present(sample_read_list.size());
    std::vector<int> misreads_eligible;
    std::string vcf_file_read_holder, vcf_file_read_match;
    std::vector<struct pos_misread*>::iterator misread_ptr;
    std::vector<struct sample_ances_list*> sample_ancestors;
    auto previous_pos = sample_map.begin()->first;
    long unsigned int map_count = 0;
    
    
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
        fill(lineage_present.begin(), lineage_present.end(), 0);
        fill(read_present.begin(), read_present.end(), 0);
 
        //Tackling Misreads not present as mutations in our sample map
        misread_ptr = misread_pos.begin();
        while (misread_ptr != misread_pos.end()) {
            auto m_e_itr = std::find(misreads_eligible.begin(), misreads_eligible.end(), (*misread_ptr)->pos);
            if ( (!(*misread_ptr)->used) && (m_e_itr == misreads_eligible.end()) && ( ( (map_count == sample_map.size()-1) && ((*misread_ptr)->pos > map.first) ) || ( (map.first == sample_map.begin()->first) && ((*misread_ptr)->pos < map.first) ) || ( (map.first != sample_map.begin()->first) && ((*misread_ptr)->pos < map.first) && ((*misread_ptr)->pos > previous_pos) ) ) ) {
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
                for (int i = 0; i < int(mut_nuc_list_reads.size()); i++) {
                    if (i) {
                        vcf_file_read_holder += ",";
                    }
                    vcf_file_read_holder += mut_nuc_list_reads[i];
                }
                vcf_file_read_holder += "\t.\t.\t.\t.\t";

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

                for (auto read: read_present)
                    vcf_file_read_holder.append(std::to_string(read) + "\t");
                vcf_file_read_holder += "\n";
                
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
                            mut_nuc_list_reads.push_back(MAT::get_nuc(misread_nuc[0]));  
                            nuc_used.emplace_back(misread_nuc[0]);
                        }
                        read_mut->read_name = (*misread_ptr)->read;
                        read_mut->mut_nuc = MAT::get_nuc(misread_nuc[0]);
                        misread_names.emplace_back(read_mut);
                        misread_nuc.clear();
                    }
                    (*misread_ptr)->used = true;
                    //std::cout << "Misread pos: " << (*misread_ptr)->pos << ", Sample: " << (*misread_ptr)->read <<  "\n";   
                    last_pos = misread;
                }   
                misread_ptr++;
            }
        }

        if (misreads_eligible.size()) {
            vcf_file_read_holder.append("\t");
            vcf_file_read_holder.push_back(ref_nuc_name_reads);
            vcf_file_read_holder.append("\t");
            for (int i = 0; i < int(mut_nuc_list_reads.size()); i++) {
                if (i) {
                    vcf_file_read_holder += ",";
                }
                vcf_file_read_holder += mut_nuc_list_reads[i];
            }
            vcf_file_read_holder += "\t.\t.\t.\t.\t";

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
            for (auto read: read_present)
                vcf_file_read_holder += std::to_string(read) + "\t";
            vcf_file_read_holder += "\n";
            
            
            if (map_count != sample_map.size()-1)
                vcf_file_reads << vcf_file_read_holder;
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
                                                mut_nuc_list_reads.push_back(MAT::get_nuc(misread_nuc[0]));  
                                                encountered_read = true;
                                            }
                                            else if (mut_nuc_list_reads.find(MAT::get_nuc(misread_nuc[0])) == std::string::npos) {
                                                vcf_file_read_match += ",";
                                                vcf_file_read_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(misread_nuc[0]);
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
                                            mut_nuc_list_reads.push_back(MAT::get_nuc(mut_anc.mut_nuc));  
                                            encountered_read = true;                     
                                        }
                                        else if (mut_nuc_list_reads.find(MAT::get_nuc(mut_anc.mut_nuc)) == std::string::npos) {
                                            vcf_file_read_match += ",";
                                            vcf_file_read_match += ref_nuc_name_reads + std::to_string(map.first) + MAT::get_nuc(mut_anc.mut_nuc); 
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
        
        if (encountered_read) {
            vcf_file_reads << vcf_file_read_match;
            vcf_file_reads << "\t" <<  ref_nuc_name_reads << "\t";
            for (int i = 0; i < int(mut_nuc_list_reads.size()); i++) {
                if (i) {
                    vcf_file_reads << ",";
                }
                vcf_file_reads << mut_nuc_list_reads[i];
            }
            vcf_file_reads << "\t.\t.\t.\t.\t";
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
            for (auto read: read_present)
                vcf_file_reads << read <<"\t";
            vcf_file_reads << "\n";
        }

        if ((map_count == sample_map.size()-1) && (misreads_eligible.size()))
            vcf_file_reads << vcf_file_read_holder;
        
        previous_pos = map.first;
        map_count ++;
        //fprintf(stderr, "Map position %ld written in both VCF \n", map_count-1);
    }

    fprintf(stderr, "VCFs Written in %ld msec \n\n", timer.Stop());

    boost::iostreams::close(outbuf_samples);
    boost::iostreams::close(outbuf_reads);
    outfile_samples.close();
    outfile_reads.close();
    }

    std::vector<struct read_info *> read_ids;
    
    timer.Start();
    read_vcf(vcf_filename_reads, read_ids);
    fprintf(stderr,"Reads_VCF parsed in %ld msec\n\n", timer.Stop());

    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    place_reads(traversal, read_ids, node_score, vcf_filename_reads);
    analyze_reads(T, traversal, node_score);
}


void place_reads(const std::vector<MAT::Node*> &dfs, const std::vector<struct read_info*> &read_ids, tbb::concurrent_hash_map<MAT::Node*, double> &node_score, std::string vcf_filename_reads) {
    fprintf(stderr, "Total nodes: %ld, Reads: %ld\n\n", dfs.size(), read_ids.size());
    timer.Start();

    tbb::concurrent_hash_map<size_t, struct min_parsimony> read_min_parsimony;
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

                                    //std::string target = rp->read;
                                    //size_t pos = target.find("_READ");
                                    //target.erase(pos);
                                    //bool there = false; 
                                    //if (dfs[i]->identifier == target) {
                                    //    for (auto mut: dfs[i]->mutations) {
                                    //        if ((mut.position >= rp->start) && (mut.position <= rp->end)) {
                                    //            fprintf(stderr, "sample mut: %d", mut.position);
                                    //            there = true;    
                                    //        }
                                    //    }
                                    //    if ((rp->mutations.empty()) && (there)) {
                                    //        fprintf(stderr, "Start: %d, End: %d, Read Mutation gone!! \n", rp->start, rp->end);
                                    //    }
                                    //}
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

                            //Else add it to curr_node_mut
                            curr_node_par_mut.emplace_back(read_mut);
                        }
                    }


                    //std::string target = rp->read;
                    //size_t pos = target.find("_READ");
                    //target.erase(pos);
                    //if (dfs[i]->identifier == target) {
                    //    fprintf(stderr, "Node: %s\n", dfs[i]->identifier.c_str());
                    //    for (auto m: curr_node_par_mut)
                    //        fprintf(stderr, "curr par mut: %d, mut_nuc: %c \n", m.position, MAT::get_nuc(m.mut_nuc));
                    //    for (auto m: rp->mutations)
                    //        fprintf(stderr, "read mut pos: %d, mut_nuc: %c\n", m.position, MAT::get_nuc(m.mut_nuc));
                    //}

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
                }

                for (auto n_idx: min_par.idx_list) {
                    auto node = dfs[n_idx];
                    tbb::concurrent_hash_map<MAT::Node*, double>::accessor ac;
                    auto created = node_score.insert(ac, std::make_pair(node, (1.0 / min_par.idx_list.size())));
                    if (!created)
                        ac->second += (1.0 / min_par.idx_list.size());
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
    
    //Check to ensure every read has its corresponding sample as its most parsimonious position
    unsigned long long avg = 0;
    using my_mutex_t = tbb::queuing_mutex;
    my_mutex_t my_mutex;
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
                        fprintf(stderr, "Target: %s\n", target.c_str());
                    }
                    fprintf(stderr, "Read: %s, mutations: %ld, Parsimony score = %ld, parsimonious positions: %ld, sample_idx = %ld \n", rp->read.c_str(), rp->mutations.size(), min_par.par_list[0].size(), min_par.par_list.size(), idx);
                    std::cout << "\n"; 
                }

                my_mutex_t::scoped_lock my_lock{my_mutex};
                    avg += (min_par.idx_list.size() / read_ids.size());
            }
        },
            ap);

    fprintf(stderr, "Avg par pos: %lld\n", avg); 
}


void analyze_reads(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, tbb::concurrent_hash_map<MAT::Node*, double> &node_score) {
    timer.Start();
    
    std::unordered_map<std::string, std::pair<double, size_t>> clade_score;
    //Print the node scores
    std::vector<std::pair<MAT::Node*, double>> top_n(50);
    std::partial_sort_copy(node_score.begin(),
                        node_score.end(),
                        top_n.begin(),
                        top_n.end(),
                        [](std::pair<const MAT::Node*, double> const& l,
                           std::pair<const MAT::Node*, double> const& r)
                        {
                            return l.second > r.second;
                        });
    
    for (auto top: top_n) {
        for (auto anc: T.rsearch(top.first->identifier, true)) { //Checking all ancestors of a node to get clade 
            bool found = false;
            //Don't consider the clades that start with numbers {Different clade naming}
            auto clade = anc->clade_annotations[1];
            if(clade != "") {
            printf("score = %f, Node: %s, Clade: %s\n", top.second, top.first->identifier.c_str(), clade.c_str());
                found = true;
            }
            if (found)
                break;
        }
    }
    
    //Making clade map
    for (auto itr = node_score.begin(); itr != node_score.end(); itr++) {
        for (auto anc: T.rsearch(itr->first->identifier, true)) { //Checking all ancestors of a node to get clade 
            bool found = false;
            //Don't consider the clades that start with numbers {Different clade naming}
            auto clade = anc->clade_annotations[1];
            if(clade != "") {
                if (clade_score.find(clade) == clade_score.end())
                    clade_score.insert({clade, {itr->second, 1}});
                else
                    clade_score[clade] = {clade_score[clade].first + itr->second, clade_score[clade].second + 1}; 
                found = true;
            }
            if (found)
                break;
        }
    }

    std::cout << "\n\n";
    std::vector<std::pair<std::string, std::pair<double, size_t>>> top_c(5);
    std::partial_sort_copy(clade_score.begin(),
                       clade_score.end(),
                       top_c.begin(),
                       top_c.end(),
                       [](std::pair<const std::string, std::pair<double, size_t>> const& l,
                          std::pair<const std::string, std::pair<double, size_t>> const& r)
                       {
                           return (l.second.first / l.second.second) > (r.second.first / r.second.second);
                       });

    for (auto top: top_c) {
        printf("Lineage: %s, Node_sum: %f, num_nodes: %zu, Avg: %f\n", top.first.c_str(), top.second.first, top.second.second, (top.second.first / top.second.second));
    }


    fprintf(stderr,"\nTop Nodes found in %ld sec\n\n", (timer.Stop() / 1000));
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
