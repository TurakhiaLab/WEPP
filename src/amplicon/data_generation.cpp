#include "data_generation.hpp"

po::variables_map parse_data_gen_command(po::parsed_options parsed) {
    po::variables_map vm;
    po::options_description conv_desc("place_read options");
    conv_desc.add_options()
    ("input-mat,i", po::value<std::string>()->required(),
     "Input mutation-annotated tree file [REQUIRED]")
    ("output-directory,o", po::value<std::string>()->default_value("./"),
     "Write output files to the target directory. Default is current directory.")
    ("write-vcf,v", po::value<std::string>()->default_value(""),
     "Output VCF file start names. Default is full tree")
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

void simulate_reads (po::parsed_options parsed) {
    //main argument for the complex extract command
    //uses included code from multiple modules
    //specifically, the modules select, describe, convert, and filter support this command
    po::variables_map vm = parse_data_gen_command(parsed);
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

    fprintf(stderr, "\nLoading input MAT files %s and %s.\n", input_mat_filename.c_str(), ref_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T;
    if ((input_mat_filename.find(".pb\0") != std::string::npos) || (ref_mat_filename.find(".pb\0") != std::string::npos)) {
        T = MAT::load_mutation_annotated_tree(input_mat_filename);
        T.uncondense_leaves();
    } else if ((input_mat_filename.find(".json\0") != std::string::npos) || (ref_mat_filename.find(".json\0") != std::string::npos)) {
        T = load_mat_from_json(input_mat_filename);
    } else {
        fprintf(stderr, "ERROR: Input file ending not recognized. Must be .json or .pb\n");
        exit(1);
    }
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));

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