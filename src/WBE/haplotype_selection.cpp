#include "wbe.hpp"

void selectHaplotypes (po::parsed_options parsed) {
    //main argument for the complex extract command
    //uses included code from multiple modules
    //specifically, the modules select, describe, convert, and filter support this command
    po::variables_map vm = parseWBEcommand(parsed);
    std::string dir_prefix = vm["output-directory"].as<std::string>();
    
    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    
    std::string input_mat_filename = dir_prefix + vm["input-mat"].as<std::string>();
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    std::string cmd_lineage = vm["lineage"].as<std::string>();
    std::string cmd_distribution = vm["distribution"].as<std::string>();
    int sample_size = vm["haplotype-samples"].as<int>();
    std::string vcf_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.vcf";
    std::string fasta_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.fa";
    
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


    timer.Start();
    //Loading reference genome
    std::ifstream fasta_f(ref_fasta);
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "\nTree loaded in %ld sec \n\n", (timer.Stop() / 1000));

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
        }
        lineage_ptr++;
    }

    for (auto sample: lineage_selected) {
        auto clade = getLineage(T, sample);
        printf("Sample: %s, Clade: %s\n", sample->identifier.c_str(), clade.c_str());
    }

    // for each sample, we need to make a vector of Mutation objects to store all the mutations that are present in the sample
    // note that each sample has the same length as the reference sequence
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

    // Removing back mutations
    // Iterate over the vector of mutations for each sample (created earlier) and remove the back mutations
    // e.g. a mutation from A to C and then back to A is a back mutation
    // we only want to keep the mutations that are not back mutations
    fprintf(stderr, "\nNumber of samples: %ld\n\n", sample_mutations_vector.size());
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

    // Printing and storing in the VCF for samples
    std::ofstream vcf_file_samples, fasta_file_samples;
    vcf_file_samples.open(vcf_filename_samples);
    vcf_file_samples << "##fileformat=VCFv4.2\n";
    vcf_file_samples << "##source=matUtils\n";
    vcf_file_samples << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (auto sample: lineage_selected) {
        vcf_file_samples << "\t" << sample->identifier;
    }
    vcf_file_samples << "\n";

    //Iterate through all the samples for writing the fasta
    fasta_file_samples.open(fasta_filename_samples);
    for (const auto& node: lineage_selected) {
        std::vector<MAT::Mutation> sample_mutations;
        for (const auto& anc_node: T.rsearch(node->identifier, true)) {
            for (const auto& mut: anc_node->mutations) {
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
        //Removing back mutations
        auto sm_itr = sample_mutations.begin();
        while (sm_itr != sample_mutations.end()) {
            if (sm_itr->mut_nuc == sm_itr->ref_nuc)
                sm_itr = sample_mutations.erase(sm_itr);
            else
                sm_itr++;
        }
        tbb::parallel_sort(sample_mutations.begin(), sample_mutations.end(), compareMutations);
        //Writing into fasta
        fasta_file_samples << ">" + node->identifier + "\n";
        sm_itr = sample_mutations.begin();
        for (int i = 1; i <= (int)ref_seq.size(); i++) {
            if (sm_itr != sample_mutations.end()) {
                if (i == sm_itr->position) {
                    fasta_file_samples << MAT::get_nuc(sm_itr->mut_nuc);
                    sm_itr++;
                }
                else
                    fasta_file_samples << ref_seq[i - 1];
            }
            else
                fasta_file_samples << ref_seq[i - 1];
        }
        fasta_file_samples << "\n"; 
        sample_mutations.clear();
    }

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
                    current_vcf_line += "\t1";
                } else {
                    current_vcf_line += "\t0";
                }
            }
            // Add the mutation to the vector of mutations that are traversed through
            mutations_observed.emplace_back(m);         
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
            else
                break;
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
    vcf_file_samples << vcf_lines_sample_all;
    vcf_file_samples.close();

    fprintf(stderr, "VCF and FASTA files for samples written in %ld msec \n\n", timer.Stop());
}