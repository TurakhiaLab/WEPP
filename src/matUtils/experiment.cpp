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
    bool old_vcf = false;
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
                end_coord = (int(ref_seq.size())) - 1;
            sample_read->read = node->identifier + "_READ_" + boost::lexical_cast<std::string>(start_coord) + "_" + boost::lexical_cast<std::string>(end_coord); 
            sample_read_list.emplace_back(sample_read);

            read_positions.clear();
            time = clock();
            srand(int(time));
            for (int pos = (rand_val - ((read_length / 2) + (read_length % 2))); pos < (rand_val + (read_length / 2)); pos ++) {
                if ((pos > 0) && (pos < (int(ref_seq.size())))) {
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

    // for each sample, we need to make a vector of Mutation objects to store all the mutations that are present in the sample
    // note that each sample has the same length as the reference sequence
    timer.Start();
    std::vector<std::vector<MAT::Mutation>> sample_mutations_vector;
    for (auto node: lineage_selected) {
        std::vector<MAT::Mutation> sample_mutations;
        auto path = T.rsearch(node->identifier, true);
        for (auto p: path) {
            for (auto m: p->mutations) {
                sample_mutations.emplace_back(m);
            }
        }
        sample_mutations_vector.emplace_back(sample_mutations);
    }
    fprintf(stderr, "Sample Mutations in %ld msec \n\n", timer.Stop());

    // Similarly, make a vector of Mutation objects to store all mutations present in the reads
    // Note that each read has the same length as the read length, and this is much smaller than the reference sequence
    // each read also has a start and end position on the reference sequence
    timer.Start();
    std::vector<std::vector<MAT::Mutation>> read_mutations_vector;
    for (auto sample_read: sample_read_list) {
        std::vector<MAT::Mutation> read_mutations;
        auto path = T.rsearch(sample_read->sample->identifier, true);
        for (auto p: path) {
            for (auto m: p->mutations) {
                // only consider mutations that are present in the read, i.e. the mutation position is between the start and end position of the read
                int start_pos = std::stoi(sample_read->read.substr(sample_read->read.find("_READ_") + 6, sample_read->read.find("_", sample_read->read.find("_READ_") + 6) - (sample_read->read.find("_READ_") + 6)));
                int end_pos = std::stoi(sample_read->read.substr(sample_read->read.find("_", sample_read->read.find("_READ_") + 6) + 1, sample_read->read.size() - (sample_read->read.find("_", sample_read->read.find("_READ_") + 6) + 1)));
                if ((m.position >= start_pos) && (m.position <= end_pos)) {
                    read_mutations.emplace_back(m);
                }
            }
        }
        read_mutations_vector.emplace_back(read_mutations);
    }

    fprintf(stderr, "Read Mutations in %ld msec \n\n", timer.Stop());


    // Removing back mutations
    // Iterate over the vector of mutations for each sample (created earlier) and remove the back mutations
    // e.g. a mutation from A to C and then back to A is a back mutation
    // we only want to keep the mutations that are not back mutations
    timer.Start();
    fprintf(stderr, "Number of samples: %ld\n\n", sample_mutations_vector.size());
    for (int i = 0; i < (int)sample_mutations_vector.size(); i++) {
        std::vector<MAT::Mutation> sample_mutations = sample_mutations_vector[i];
        std::vector<std::int64_t> back_mutation_indices;
        for (int j = 0; j < (int)sample_mutations.size(); j++) {
            for (int k = 0; k < (int)sample_mutations.size(); k++) {
                auto m = sample_mutations[j];
                auto m2 = sample_mutations[k];
                if ((std::find(back_mutation_indices.begin(), back_mutation_indices.end(), j) != back_mutation_indices.end()) || (std::find(back_mutation_indices.begin(), back_mutation_indices.end(), k) != back_mutation_indices.end())) {
                    continue;
                }
                if ((m.chrom == m2.chrom) && (m.position == m2.position) && (m.par_nuc == m2.mut_nuc) && (m.mut_nuc == m2.par_nuc)) {
                    //fprintf(stderr, "Back Mutation: Chrom: %s, Pos: %d, Ref: %d, Mut: %d\n", m.chrom.c_str(), m.position, m.ref_nuc, m.mut_nuc );
                    // add the index of the back mutation to the vector of back mutation indices
                    back_mutation_indices.emplace_back(j);
                    back_mutation_indices.emplace_back(k);
                    break;
                }
            }
        }
        std::vector<MAT::Mutation> new_sample_mutations;
        for (int j = 0; j < (int)sample_mutations.size(); j++) {
            if (std::find(back_mutation_indices.begin(), back_mutation_indices.end(), j) == back_mutation_indices.end()) {
                new_sample_mutations.emplace_back(sample_mutations[j]);
            }
        }
        sample_mutations_vector[i] = new_sample_mutations;
    }
    
    fprintf(stderr, "Sample Back Mutations removed in %ld msec \n\n", timer.Stop());

        // Remove back mutations from the read mutations

    timer.Start();

    // print the number of reads
    fprintf(stderr, "Number of reads: %ld\n\n", read_mutations_vector.size());

    for (int i = 0; i < (int)read_mutations_vector.size(); i++) {
        std::vector<MAT::Mutation> read_mutations = read_mutations_vector[i];
        std::vector<std::int64_t> back_mutation_indices;

        for (int j = 0; j < (int)read_mutations.size(); j++) {
            for (int k = 0; k < (int)read_mutations.size(); k++) {
                auto m = read_mutations[j];
                auto m2 = read_mutations[k];
                if ((std::find(back_mutation_indices.begin(), back_mutation_indices.end(), j) != back_mutation_indices.end()) || (std::find(back_mutation_indices.begin(), back_mutation_indices.end(), k) != back_mutation_indices.end())) {
                    continue;
                }
                if ((m.chrom == m2.chrom) && (m.position == m2.position) && (m.par_nuc == m2.mut_nuc) && (m.mut_nuc == m2.par_nuc)) {
                    //fprintf(stderr, "Back Mutation: Chrom: %s, Pos: %d, Ref: %d, Mut: %d\n", m.chrom.c_str(), m.position, m.ref_nuc, m.mut_nuc );
                    // add the index of the back mutation to the vector of back mutation indices
                    back_mutation_indices.emplace_back(j);
                    back_mutation_indices.emplace_back(k);
                    break;
                }
            }
        }
        std::vector<MAT::Mutation> new_read_mutations;
        for (int j = 0; j < (int)read_mutations.size(); j++) {
            if (std::find(back_mutation_indices.begin(), back_mutation_indices.end(), j) == back_mutation_indices.end()) {
                new_read_mutations.emplace_back(read_mutations[j]);
            }
        }
        read_mutations_vector[i] = new_read_mutations;
    }

    

    // Account for read error, by going over the misread positions
    // if the misread position has no mutations, then add the mutation to the read mutations
    // if the misread position has a mutation, then update the most recent mutation at that position to the misread mutation
    for (auto misread: misread_pos) {
        for (int i = 0; i < (int)read_mutations_vector.size(); i++) {
            if (misread->read != sample_read_list[i]->read) {
                continue;
            }
            std::vector<MAT::Mutation> read_mutations = read_mutations_vector[i];
            bool misreadpos = false;
            for (auto m: read_mutations) {
                if ((m.chrom == "NC_045512v2") && (m.position == misread->pos)) {
                    misreadpos = true;
                    break;
                }
            }
            if (!misreadpos) {
                MAT::Mutation new_mutation;
                new_mutation.chrom = "NC_045512v2";
                new_mutation.position = misread->pos;
                new_mutation.ref_nuc = ref_seq[misread->pos];
                new_mutation.mut_nuc = ref_seq[misread->pos];
                while(new_mutation.mut_nuc == ref_seq[misread->pos]) {
                    new_mutation.mut_nuc = MAT::get_nuc(nuc_array[rand() % nuc_array.size()]);
                    
                }
                read_mutations.emplace_back(new_mutation);
                read_mutations_vector[i] = read_mutations;
            } else {
                for (int j = 0; j < (int)read_mutations.size(); j++) {
                    auto m = read_mutations[j];
                    if ((m.chrom == "NC_045512v2") && (m.position == misread->pos)) {
                        read_mutations[j].mut_nuc = ref_seq[misread->pos];
                        while(read_mutations[j].mut_nuc == ref_seq[misread->pos]) {
                            read_mutations[j].mut_nuc = MAT::get_nuc(nuc_array[rand() % nuc_array.size()]);
                        }
                        read_mutations_vector[i] = read_mutations;
                        break;
                    }
                }
            }
        }
    }
    
   // fprintf(stderr, "Read Back Mutations removed in %ld msec \n\n", timer.Stop());

    // Printing and storing in the VCF for samples
    timer.Start();
    // print out the VCF header to the VCF file for samples
    std::ofstream vcf_file_samples_test;
    vcf_file_samples_test.open(vcf_filename_samples);
    vcf_file_samples_test << "##fileformat=VCFv4.2\n";
    vcf_file_samples_test << "##source=matUtils\n";
    vcf_file_samples_test << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (auto sample: lineage_selected) {
        vcf_file_samples_test << "\t" << sample->identifier;
    }
    vcf_file_samples_test << "\n";

    std::vector<std::string> vcf_file_samples_lines_list;

    // make a vector of mutations to store the mutations that are traversed through
    std::vector<MAT::Mutation> mutations_observed;
    // Now, add the sample mutations to the VCF file
    for (int i = 0; i < (int)sample_mutations_vector.size(); i++) {
        std::vector<MAT::Mutation> sample_mutations = sample_mutations_vector[i];
        for (auto m: sample_mutations) {
            // check if mutation has already been traversed through
            std::string current_vcf_line = "";
            bool mutation_already_traversed = false;
            for (auto m2: mutations_observed) {
                if ((m.chrom == m2.chrom) && (m.position == m2.position) && (m.ref_nuc == m2.ref_nuc) && (m.mut_nuc == m2.mut_nuc)) {
                    mutation_already_traversed = true;
                    break;
                }
            }
            if (mutation_already_traversed) {
                continue;
            }
            //vcf_file_samples << m.chrom << "\t" << m.position << "\t" << MAT::get_nuc(m.ref_nuc) << m.position << MAT::get_nuc(m.mut_nuc) << "\t" << MAT::get_nuc(m.ref_nuc) << "\t" << MAT::get_nuc(m.mut_nuc) << "\t.\t.\t.\t.";
            current_vcf_line += m.chrom + "\t" + std::to_string(m.position) + "\t" + MAT::get_nuc(m.ref_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc) + "\t" + MAT::get_nuc(m.ref_nuc) + "\t" + MAT::get_nuc(m.mut_nuc) + "\t.\t.\t.\t.";
            // Loop over all samples and check if the mutation is present in the sample
            for (int j = 0; j < (int)lineage_selected.size(); j++) {
                bool mutation_present = false;
                for (auto m2: sample_mutations_vector[j]) {
                    if ((m.chrom == m2.chrom) && (m.position == m2.position) && (m.ref_nuc == m2.ref_nuc) && (m.mut_nuc == m2.mut_nuc)) {
                        mutation_present = true;
                        break;
                    }
                }
                if (mutation_present) {
                    //vcf_file_samples << "\t1";
                    current_vcf_line += "\t1";
                } else {
                    //vcf_file_samples << "\t0";
                    current_vcf_line += "\t0";
                }
            }
            // Add the mutation to the vector of mutations that are traversed through
            mutations_observed.emplace_back(m);         
            //vcf_file_samples << "\n";
            vcf_file_samples_lines_list.emplace_back(current_vcf_line);
        }
    }

    // sort the lines by position
    std::sort(vcf_file_samples_lines_list.begin(), vcf_file_samples_lines_list.end(), [](const std::string& a, const std::string& b) {
        std::vector<std::string> a_words, b_words;
        MAT::string_split(a, a_words);
        MAT::string_split(b, b_words);
        return std::stoi(a_words[1]) < std::stoi(b_words[1]);
    });

    // remove lines where the ref and alt are the same
    std::vector<std::string> vcf_file_samples_lines_list_filtered;
    for (auto line: vcf_file_samples_lines_list) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words[3] != words[4]) {
            vcf_file_samples_lines_list_filtered.emplace_back(line);
        }
    }
    vcf_file_samples_lines_list = vcf_file_samples_lines_list_filtered;

    std::vector<std::string> vcf_lines_merged;
    for (int i = 0; i < (int)vcf_file_samples_lines_list.size(); i++) {
        std::vector<std::string> words;
        MAT::string_split(vcf_file_samples_lines_list[i], words);
        int pos = std::stoi(words[1]);
        std::string alt = words[4];
        std::string id = words[2];
        std::vector<int> sample_i_counts;
        
        for (int k = 9; k < (int)words.size(); k++) {
            sample_i_counts.emplace_back(std::stoi(words[k]));
        }

        int j;

        for (j = i + 1; j < (int)vcf_file_samples_lines_list.size(); j++) {
            std::vector<std::string> words2;
            MAT::string_split(vcf_file_samples_lines_list[j], words2);
            int pos2 = std::stoi(words2[1]);
            std::string alt2 = words2[4];
            std::string id2 = words2[2];
            if (pos == pos2) {
                alt = alt + "," + alt2;
                id = id + "," + id2;
                for (int k = 9; k < (int)words2.size(); k++) {
                    if (words[k] == "1") {
                        sample_i_counts[k - 9] = 1;
                    }
                    if (words2[k] == "1") {
                        sample_i_counts[k - 9] = j - i + 1;
                    }
                
            }
            } 
            
            else {
                break;
            }
        }

        i = j - 1;

        vcf_lines_merged.emplace_back(words[0] + "\t" + words[1] + "\t" + id + "\t" + words[3] + "\t" + alt + "\t" + words[5] + "\t" + words[6] + "\t" + words[7] + "\t" + words[8]);
        for (int k = 0; k < (int)sample_i_counts.size(); k++) {
            vcf_lines_merged.back() += "\t" + std::to_string(sample_i_counts[k]);
        }
    }
    std::string vcf_lines_sample_all = "";
    for (auto line: vcf_lines_merged) {
        vcf_lines_sample_all += line + "\n";
    }
    vcf_file_samples_test << vcf_lines_sample_all;
    vcf_file_samples_test.close();

    fprintf(stderr, "VCF file for samples written in %ld msec \n\n", timer.Stop());

    // Print and store in the VCF for reads
    timer.Start();

    // print out the VCF header to the VCF file for reads
    std::ofstream vcf_file_reads_test;
    vcf_file_reads_test.open(vcf_filename_reads);
    vcf_file_reads_test << "##fileformat=VCFv4.2\n";
    vcf_file_reads_test << "##source=matUtils\n";
    vcf_file_reads_test << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (auto sample_read: sample_read_list) {
        vcf_file_reads_test << "\t" << sample_read->read;
    }
    vcf_file_reads_test << "\n";
    std::vector<std::string> vcf_file_reads_lines_list;
    // make a vector of mutations to store the mutations that are traversed through
    mutations_observed.clear();
    // Now, add the read mutations to the VCF file
    for (int i = 0; i < (int)read_mutations_vector.size(); i++) {
        std::vector<MAT::Mutation> read_mutations = read_mutations_vector[i];
        for (auto m: read_mutations) {
            std::string current_vcf_line = "";
            // check if mutation has already been traversed through
            bool mutation_already_traversed = false;
            for (auto m2: mutations_observed) {
                if ((m.chrom == m2.chrom) && (m.position == m2.position) && (m.ref_nuc == m2.ref_nuc) && (m.mut_nuc == m2.mut_nuc)) {
                    mutation_already_traversed = true;
                    break;
                }
            }
            if (mutation_already_traversed) {
                continue;
            }
            // vcf_file_reads << m.chrom << "\t" << m.position << "\t" << MAT::get_nuc(m.ref_nuc) << m.position << MAT::get_nuc(m.mut_nuc) << "\t" << MAT::get_nuc(m.ref_nuc) << "\t" << MAT::get_nuc(m.mut_nuc) << "\t.\t.\t.\t.";
            current_vcf_line += m.chrom + "\t" + std::to_string(m.position) + "\t" + MAT::get_nuc(m.ref_nuc) + std::to_string(m.position) + MAT::get_nuc(m.mut_nuc) + "\t" + MAT::get_nuc(m.ref_nuc) + "\t" + MAT::get_nuc(m.mut_nuc) + "\t.\t.\t.\t.";
            // Loop over all reads and check if the mutation is present in the read
            for (int j = 0; j < (int)sample_read_list.size(); j++) {
                bool mutation_present = false;
                for (auto m2: read_mutations_vector[j]) {
                    if ((m.chrom == m2.chrom) && (m.position == m2.position) && (m.ref_nuc == m2.ref_nuc) && (m.mut_nuc == m2.mut_nuc)) {
                        mutation_present = true;
                        break;
                    }
                }
                if (mutation_present) {
                    // vcf_file_reads << "\t1";
                    current_vcf_line += "\t1";
                } else {
                    //vcf_file_reads << "\t0";
                    current_vcf_line += "\t0";
                }
            }

            // Add the mutation to the vector of mutations that are traversed through
            mutations_observed.emplace_back(m);         
            //vcf_file_reads << "\n";
            vcf_file_reads_lines_list.emplace_back(current_vcf_line);
        }
    }

    // sort the lines by position
    std::sort(vcf_file_reads_lines_list.begin(), vcf_file_reads_lines_list.end(), [](const std::string& a, const std::string& b) {
        std::vector<std::string> a_words, b_words;
        MAT::string_split(a, a_words);
        MAT::string_split(b, b_words);
        return std::stoi(a_words[1]) < std::stoi(b_words[1]);
    });

    // remove lines where the ref and alt are the same
    std::vector<std::string> vcf_file_reads_lines_list_filtered;
    for (auto line: vcf_file_reads_lines_list) {
        std::vector<std::string> words;
        MAT::string_split(line, words);
        if (words[3] != words[4]) {
            vcf_file_reads_lines_list_filtered.emplace_back(line);
        }
    }

    vcf_file_reads_lines_list = vcf_file_reads_lines_list_filtered;
    
    vcf_lines_merged.clear();

    for (int i = 0; i < (int)vcf_file_reads_lines_list.size(); i++) {
        std::vector<std::string> words;
        MAT::string_split(vcf_file_reads_lines_list[i], words);
        int pos = std::stoi(words[1]);
        std::string alt = words[4];
        std::string id = words[2];
        std::vector<int> sample_i_counts;
        for (int k = 9; k < (int)words.size(); k++) {
            sample_i_counts.emplace_back(std::stoi(words[k]));
        }
        
        int j;

        for (j = i + 1; j < (int)vcf_file_reads_lines_list.size(); j++) {
            std::vector<std::string> words2;
            MAT::string_split(vcf_file_reads_lines_list[j], words2);
            int pos2 = std::stoi(words2[1]);
            std::string alt2 = words2[4];
            std::string id2 = words2[2];
            if (pos == pos2) {
                alt = alt + "," + alt2;
                id = id + "," + id2;
                for (int k = 9; k < (int)words2.size(); k++) {
                    if (words[k] == "1") {
                        sample_i_counts[k - 9] = 1;
                    }
                    if (words2[k] == "1") {
                        sample_i_counts[k - 9] = j - i + 1;
                    }
                
            }
            } 
            
            else {
                break;
            }
        }

        i = j - 1;

        vcf_lines_merged.emplace_back(words[0] + "\t" + words[1] + "\t" + id + "\t" + words[3] + "\t" + alt + "\t" + words[5] + "\t" + words[6] + "\t" + words[7] + "\t" + words[8]);
        for (int k = 0; k < (int)sample_i_counts.size(); k++) {
            vcf_lines_merged.back() += "\t" + std::to_string(sample_i_counts[k]);
        }
    }

    std::string vcf_lines_read_all = "";
    for (auto line: vcf_lines_merged) {
        vcf_lines_read_all += line + "\n";
    }

    vcf_file_reads_test << vcf_lines_read_all;
    vcf_file_reads_test.close();
    fprintf(stderr, "VCF file for reads written in %ld msec \n\n", timer.Stop());

    // Generating VCF file for reads --> but for Freyja's input
    // Now, for each line in the VCF, we add up the total number of non-zeros and divide by the total number of reads
    // This ratio is stored in another VCF file
    timer.Start();
    std::ofstream vcf_file_reads_freyja_test;
    vcf_file_reads_freyja_test.open(vcf_filename_reads_freyja);
    vcf_file_reads_freyja_test << "##fileformat=VCFv4.2\n";
    vcf_file_reads_freyja_test << "##source=matUtils\n";
    vcf_file_reads_freyja_test << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    // QUAL and FILTER will be set to . and . respectively for each row
    // but INFO will be set to AF=<ratio>, where <ratio> is the ratio of 1's to total number of reads

    // the lines are already sorted by position
    // so we can just iterate over the lines and calculate the ratio for each line

    // group the vcf lines by position, making a vector of vectors
    std::vector<std::vector<std::string>> vcf_lines_grouped;
    std::vector<std::string> vcf_lines_grouped_temp;

    for (int i = 0; i < (int)vcf_file_reads_lines_list.size(); i++) {
        std::vector<std::string> words;
        MAT::string_split(vcf_file_reads_lines_list[i], words);
        int pos = std::stoi(words[1]);
        vcf_lines_grouped_temp.emplace_back(vcf_file_reads_lines_list[i]);
        for (int j = i + 1; j < (int)vcf_file_reads_lines_list.size(); j++) {
            std::vector<std::string> words2;
            MAT::string_split(vcf_file_reads_lines_list[j], words2);
            int pos2 = std::stoi(words2[1]);
            if (pos == pos2) {
                vcf_lines_grouped_temp.emplace_back(vcf_file_reads_lines_list[j]);
                i = j + 1;
            } else {
                break;
            }
        }
        vcf_lines_grouped.emplace_back(vcf_lines_grouped_temp);
        vcf_lines_grouped_temp.clear();
    }

    // Make a vector of strings to store the lines for the VCF file for reads for Freyja
    std::vector<std::string> vcf_lines_freyja_list;

    for (auto group: vcf_lines_grouped) {
        std::vector<std::string> alts;
        std::vector<std::string> ids;
        std::vector<std::string> infos;
        std::string current_vcf_line = "";

        // store chrom, pos, ref, qual, filter for the first line in the group
        std::vector<std::string> words_first;
        MAT::string_split(group[0], words_first);

        for (auto line: group) {
            std::vector<std::string> words;
            MAT::string_split(line, words);

            int num_ones = 0;
            for (int i = 9; i < (int)words.size(); i++) {
                if (words[i] != "0") {
                    num_ones++;
                }
            }
            int num_reads = 0;
            for (auto sample_read: sample_read_list) {
                int start_pos = std::stoi(sample_read->read.substr(sample_read->read.find("_READ_") + 6, sample_read->read.find("_", sample_read->read.find("_READ_") + 6) - (sample_read->read.find("_READ_") + 6)));
                int end_pos = std::stoi(sample_read->read.substr(sample_read->read.find("_", sample_read->read.find("_READ_") + 6) + 1, sample_read->read.size() - (sample_read->read.find("_", sample_read->read.find("_READ_") + 6) + 1)));
                if ((std::stoi(words[1]) >= start_pos) && (std::stoi(words[1]) <= end_pos)) {
                    num_reads++;
                }
            }
            double ratio = (double)num_ones / (double)num_reads;
            infos.emplace_back("AF=" + std::to_string(ratio));
            alts.emplace_back(words[4]);
            ids.emplace_back(words[2]);
        }

        // vcf_file_reads_freyja << words_first[0] << "\t" << words_first[1] << "\t";
        current_vcf_line += words_first[0] + "\t" + words_first[1] + "\t";
        for (int i = 0; i < (int)ids.size(); i++) {
            if (i == 0) {
                current_vcf_line += ids[i];
            } else {
                current_vcf_line += "," + ids[i];
            }
        }
        current_vcf_line += "\t" + words_first[3] + "\t";
        for (int i = 0; i < (int)alts.size(); i++) {
            if (i == 0) {
                current_vcf_line += alts[i];
            } else {
                current_vcf_line += "," + alts[i];
            }
        }
        current_vcf_line += "\t.\t.\t";
        for (int i = 0; i < (int)infos.size(); i++) {
            if (i == 0) {
                current_vcf_line += infos[i];
            } else {
                current_vcf_line += ";" + infos[i];
            }
        }
        current_vcf_line += "\n";
        vcf_lines_freyja_list.emplace_back(current_vcf_line);
    }

    // combine all the lines for VCF freyja into one string, and write all the lines at once
    std::string vcf_lines_freyja_combined;
    for (auto line: vcf_lines_freyja_list) {
        vcf_lines_freyja_combined += line;
    }
    vcf_file_reads_freyja_test << vcf_lines_freyja_combined;
    vcf_file_reads_freyja_test.close();
    fprintf(stderr, "VCF file for reads for Freyja written in %ld msec \n\n", timer.Stop());

    timer.Start();
    std::ofstream vcf_file_reads_freyja_depth;
    vcf_file_reads_freyja_depth.open(depth_filename_reads_freyja);

    // Iterate over the length of the reference sequence
    // For each position, count the number of reads that cover that position
    // and write the position and the number of reads to the depth file

    std::string reference = "NC_045512v2";
    // Make a variable to store all the lines for the depth file, so we can write them all at once
    std::vector<std::string> depth_lines;
    for (int i = 0; i < (int)ref_seq.size(); i++) {
        int read_count = 0;
        for (auto sample_read: sample_read_list) {
            int start_pos = std::stoi(sample_read->read.substr(sample_read->read.find("_READ_") + 6, sample_read->read.find("_", sample_read->read.find("_READ_") + 6) - (sample_read->read.find("_READ_") + 6)));
            int end_pos = std::stoi(sample_read->read.substr(sample_read->read.find("_", sample_read->read.find("_READ_") + 6) + 1, sample_read->read.size() - (sample_read->read.find("_", sample_read->read.find("_READ_") + 6) + 1)));
            if ((i >= start_pos) && (i <= end_pos)) {
                read_count++;
            }
        }
        //vcf_file_reads_freyja_depth << reference << "\t" << i+1 << "\t" << ref_seq[i] << "\t" << read_count << "\n";
        depth_lines.emplace_back(reference + "\t" + std::to_string(i+1) + "\t" + ref_seq[i] + "\t" + std::to_string(read_count));

    }

    // combine the lines into one string, and write all the lines at once
    std::string depth_lines_combined;
    for (auto line: depth_lines) {
        depth_lines_combined += line + "\n";
    }
    vcf_file_reads_freyja_depth << depth_lines_combined;

    vcf_file_reads_freyja_depth.close();
    fprintf(stderr, "Depth file for reads for Freyja written in %ld msec \n\n", timer.Stop());
    }

    //Checking how close are input samples with peaks
    std::unordered_map<int, struct read_info*> read_map;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    std::vector<std::string> vcf_samples;
    read_sample_vcf(vcf_samples, vcf_filename_samples);
    read_vcf(read_map, vcf_filename_reads);
    analyze_reads(T, T_ref, dfs, read_map, node_score, vcf_samples, mismatch_matrix_file, barcode_file, read_abundance_vcf);
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


void read_vcf(std::unordered_map<int, struct read_info*> &read_map, const std::string vcf_filename_reads) {
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

    fprintf(stderr,"%s parsed in %ld sec\n\n", vcf_filename_reads.c_str(), (timer.Stop() / 1000));   
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
        }
        //Place as a child if current node is not a leaf node or ROOT node
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
        }

        //Updating the parsimony score of peak nodes for EM algorithm
        for (int j = 0; j < (int)peak_nodes.size(); j++) {
            if (peak_nodes[j] == curr_node) {
                mismatch_vector[j] = (int)curr_node_par_mut.size();
                peak_count++;
                if (peak_count == (int)peak_nodes.size())
                    return 0;
                break;
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
        //VERIFY only if searching for nodes with zero parsimony score
        if (!par_score_lim) {
            //Check to ensure every read has its corresponding sample as its most parsimonious position
            std::string target = rp->read;
            size_t pos = target.find("_READ");
            target.erase(pos);
            auto target_p = T.get_node(target)->parent->identifier;
            bool found = false;
            int idx, i;
            for (i = 0; i < (int)min_par.idx_list.size(); i++) {
                idx = min_par.idx_list[i];
                if ((dfs[idx]->identifier == target) || (dfs[idx]->identifier == target_p)) {
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

                    fprintf(stderr, "read mut: ");
                    for (auto mut: rp->mutations)
                        fprintf(stderr, "%c%d%c, ", MAT::get_nuc(mut.ref_nuc), mut.position, MAT::get_nuc(mut.mut_nuc));
                    fprintf(stderr, "\n");

                }
                fprintf(stderr, "Read: %s, read mutations: %ld, Parsimony score = %ld, parsimonious positions: %ld\n\n", rp->read.c_str(), rp->mutations.size(), min_par.par_list[0].size(), min_par.par_list.size());
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


void analyze_reads(const MAT::Tree &T, const MAT::Tree &T_ref, const std::vector<MAT::Node*> &dfs, const std::unordered_map<int, struct read_info*> &read_map, tbb::concurrent_hash_map<MAT::Node*, double> &node_score, const std::vector<std::string> &vcf_samples, const std::string &mismatch_matrix_file, const std::string &barcode_file, const std::string &read_abundance_vcf) {
    timer.Start();
    int top_n = 25, m_dist_thresh = 2, neighbor_dist_thresh = 5, neighbor_peaks_thresh = 20;
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
    while ((int)remaining_reads.size() > 0) {
        //Calculating node score for remaining reads
        static tbb::affinity_partitioner ap;
        tbb::parallel_for( tbb::blocked_range<size_t>(0, remaining_reads.size()),
            [&](tbb::blocked_range<size_t> k) {
                for (size_t i = k.begin(); i < k.end(); ++i) {
                    int rm_idx = remaining_reads[i];
                    auto read_id = read_map.find(rm_idx)->second;
                    place_reads(T_ref, dfs, read_id, NULL, node_score, peak_nodes_dummy, mismatch_vector_dummy, 0);
                }
            },
        ap);
        break;
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

        //Find top node not seen before
        std::vector<bool> peak_vec(top_n_node_score.size(), true);
        bool break_loop = false;
        auto top_n_itr = top_n_node_score.begin();
        //Find top_node that is NOT present in neighborhood of peak_nodes
        while (top_n_itr != top_n_node_score.end()) {
            bool present = check_peaks_neighbourhood(T, top_n_itr->first, peak_nodes, m_dist_thresh);
            if (!present)
                break;
            peak_vec[top_n_itr - top_n_node_score.begin()] = false;
            top_n_itr++;
            //If reach node_score.size() limit then quit
            if ((top_n_itr - top_n_node_score.begin()) == (int)node_score.size()) {
                break_loop = true;
                break;
            }
        }
        node_score.clear();
        //Stop the iterative loop if no unique peak if found
        if (break_loop)
            break;

        auto top_score = top_n_itr->second;
        //Find the reads not mapped to the best nodes
        while (top_n_itr != top_n_node_score.end()) {
            //Peak Finding
            auto curr_node = top_n_itr->first;
            //Only add unique nodes in peak_nodes with score == top_score
            if (abs(top_score - top_n_itr->second) < 1e-9) {
                //Check if curr_node does not lie in neighborhood of any node seen in current iterarion or previous
                bool present = check_peaks_neighbourhood(T, curr_node, peak_nodes, m_dist_thresh);
                if ((peak_vec[top_n_itr - top_n_node_score.begin()]) && (!present)) {
                    peak_nodes.emplace_back(curr_node);
                    //Remove mapped reads from remaining_reads
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
                                int read_present = place_reads(T_ref, dfs, read_id, curr_node, node_score, peak_nodes_dummy, mismatch_vector_dummy, 0);
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
                int m_dist = mutation_distance(T, T, curr_node, cmp_node);
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
    }
    remaining_reads.clear();
    printf("\nInital PEAK nodes: %d\n", (int)peak_nodes.size());
    fprintf(stderr,"Peak search took %ld min\n\n", (timer.Stop() / (60 * 1000)));
    
    printf("MUTATION DISTANCE ORIG:\n"); 
    //Verify Recovery of Input Samples
    timer.Start();
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        //auto ref_clade = get_clade(T, T.get_node(sample));
        //auto best_clade = ref_clade;
        auto best_node = T_ref.get_node(sample);
        for (auto pn: peak_nodes) {
            int curr_dist = mutation_distance(T, T_ref, pn, T_ref.get_node(sample));
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                //best_clade = get_clade(T_ref, pn);
                best_node = pn;
            }
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node->identifier.c_str(), min_dist);
    }
    fprintf(stderr,"Mutation Distance Verification took %ld msec\n\n", timer.Stop());

    add_neighbor_peaks(T, peak_nodes, neighbor_dist_thresh, neighbor_peaks_thresh);
    generate_regression_abundance_data(T, peak_nodes, read_map, barcode_file, read_abundance_vcf);
    //generate_EM_data(T, dfs, read_map, peak_nodes, mismatch_matrix_file);
}


//Add neighboring nodes to peaks
void add_neighbor_peaks(const MAT::Tree &T, std::vector<MAT::Node*> &peak_nodes, const int neighbor_dist_thresh, const int neighbor_peaks_thresh) {
    std::vector<MAT::Node*> neighbor_nodes;
    for (const auto& peak: peak_nodes) {
        //Find the farthest ancestor within neighbor_dist_thresh
        auto anc_node = peak;
        for (auto n: T.rsearch(peak->identifier, false)) {
            //Get mutation distance between ancestor and curr_node
            int m_dist = mutation_distance(T, T, n, peak);
            if ((m_dist <= neighbor_dist_thresh) && (!n->is_root()))
                anc_node = n;
            else
                break;
        }
        //Add neighborhood peaks to peak_nodes
        neighbor_nodes.clear();
        if (anc_node != NULL) {
            std::queue<Mutation_Annotated_Tree::Node*> remaining_nodes;
            remaining_nodes.push(anc_node);
            while (remaining_nodes.size() > 0) {
                Mutation_Annotated_Tree::Node* present_node = remaining_nodes.front();
                remaining_nodes.pop();
                //Only add present_node if within neighbor_dist_thresh
                int m_dist = mutation_distance(T, T, present_node, peak);
                if ((m_dist <= neighbor_dist_thresh) && ((int)neighbor_nodes.size() < neighbor_peaks_thresh)) {
                    bool found = false;
                    //Check if a peak with same mutations is not already present
                    for (const auto& hap: peak_nodes) {
                        int mut_dist = mutation_distance(T, T, hap, present_node);
                        if (!mut_dist) {
                            found = true;
                            break;
                        }
                    }
                    //Check if neighbor peak with same mutations already accounted
                    if (!found) {
                        for (const auto& hap: neighbor_nodes) {
                            int mut_dist = mutation_distance(T, T, hap, present_node);
                            if (!mut_dist) {
                                found = true;
                                break;
                            }
                        }
                    }
                    //Only add if the node is unique
                    if (!found)
                        neighbor_nodes.emplace_back(present_node);
                }
                else if ((int)neighbor_nodes.size() == neighbor_peaks_thresh) {
                    while (!remaining_nodes.empty())
                        remaining_nodes.pop();
                    break;
                }
                //Add current node's children in list for checking
                for (auto c: present_node->children)
                    remaining_nodes.push(c);
            }
        }
        peak_nodes.insert(peak_nodes.end(), neighbor_nodes.begin(), neighbor_nodes.end());
    }
}


//Function to check if current node is in neighbourhood of peak nodes
bool check_peaks_neighbourhood (const MAT::Tree &T, const MAT::Node* N, const std::vector<MAT::Node*> &peak_nodes, const int m_dist_thresh) {
    for (auto pn: peak_nodes) {
        //Return false if the Node lies within mutation distance limit 
        int m_dist = mutation_distance(T, T, pn, N);
        if (m_dist <= m_dist_thresh)
            return true;
    }
    return false;
}


//Function to calculation distance between two nodes
int mutation_distance(const MAT::Tree &T1, const MAT::Tree &T2, const MAT::Node* N1, const MAT::Node* N2) {
    std::vector<MAT::Mutation> node1_mutations, node2_mutations;
    for (auto anc: T1.rsearch(N1->identifier, true)) { //Checking all ancestors of a node
        for (auto mut: anc->mutations) {
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
    for (auto anc: T2.rsearch(N2->identifier, true)) { //Checking all ancestors of a node
        for (auto mut: anc->mutations) {
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
    //Remove Back Mutations
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

    std::sort(node1_mutations.begin(), node1_mutations.end(), compare_mutations);
    std::sort(node2_mutations.begin(), node2_mutations.end(), compare_mutations);
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

//Generate EM based estimate algorithm data
void generate_EM_data(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, const std::unordered_map<int, struct read_info*> &read_map, const std::vector<MAT::Node*> &peak_nodes, const std::string &mismatch_matrix_file) {
    //Write MISMATCH File for Abundance Estimation using EM algorithm
    timer.Start();
    std::unordered_map<int, std::vector<int>> mismatch_matrix;
    tbb::concurrent_hash_map<MAT::Node*, double> node_score;
    std::unordered_map<std::string, double> hap_abun_map_dummy;
    tbb::concurrent_hash_map<std::string, double> lineage_score_dummy;
    std::ofstream outfile_mismatch_matrix(mismatch_matrix_file, std::ios::out | std::ios::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_mismatch_matrix;
    if (mismatch_matrix_file.find(".gz\0") != std::string::npos) {
            outbuf_mismatch_matrix.push(boost::iostreams::gzip_compressor());
    }
    outbuf_mismatch_matrix.push(outfile_mismatch_matrix);
    std::ostream mismatch_file(&outbuf_mismatch_matrix);
    std::string file_write;

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

    //Write the header containing peak nodes and corresponding lineages
    for (int i = 0; i < (int)peak_nodes.size(); i++) {
        if (i)
            file_write += "," + peak_nodes[i]->identifier + "_" + get_clade(T, peak_nodes[i]);
        else
            file_write += peak_nodes[i]->identifier + "_" + get_clade(T, peak_nodes[i]);
    } 
    file_write += "\n";
    
    auto m_itr = mismatch_matrix.begin();
    while (m_itr != mismatch_matrix.end()) {
        auto m_element = m_itr->second;
        int zero_count = 0;
        //Only write the vectos where atleast one 1 is present
        for (auto m: m_element) {
            if (!m)
                zero_count++;
            else 
                break;
        }
        if (zero_count < (int)m_element.size()) {
            for (int i = 0; i < (int)m_element.size(); i++) {
                if (i)
                    file_write += "," + std::to_string(m_element[i]);
                else 
                    file_write += std::to_string(m_element[i]);
            }
            file_write += "\n";
        }
        m_itr++;
    }
    mismatch_file << file_write;
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

    //Add unique mutations captured from read_vcf first
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

    std::sort(peak_mut_list.begin(), peak_mut_list.end(), compare_mutations);
    //Writing the header mutations
    std::string barcode_print, vcf_print;
    vcf_print += "##fileformat=VCFv4.2\n##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tDEPTH\n";
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
        float af = 0.0;
        if (total_reads > 0)
            af = (float)match_reads / (float)total_reads;
        vcf_print += "NC_045512v2\t" + std::to_string(mut.position) + "\t" + MAT::get_nuc(mut.par_nuc) + std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc) + "\t" + MAT::get_nuc(mut.par_nuc) + "\t" + MAT::get_nuc(mut.mut_nuc) + "\t.\t.\tAF=";
        vcf_print += std::to_string(af) + "\t" + std::to_string(total_reads) + "\n";
    } 
    barcode << barcode_print;
    vcf << vcf_print;
    barcode_print.clear();
    vcf_print.clear();

    //Writing peak mutations
    peak_mut_itr = peak_mut_map.begin();
    while (peak_mut_itr != peak_mut_map.end()) {
        barcode_print += "\n";
        barcode_print += peak_mut_itr->first->identifier + "_" + get_clade(T, peak_mut_itr->first);
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
                barcode_print += ",0";
                idx++;
            }
            if (idx == m_idx) {
                barcode_print += ",1";
                idx++;
            }
        }
        while (idx < (int)peak_mut_list.size()) {
            barcode_print += ",0";
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