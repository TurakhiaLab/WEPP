#include "wbe.hpp"

void refinePeaks(po::parsed_options parsed) {
    //main argument for the complex extract command
    po::variables_map vm = parseWBEcommand(parsed);
    std::string input_mat_filename = vm["input-mat"].as<std::string>();
    std::string dir_prefix = vm["output-directory"].as<std::string>();

    boost::filesystem::path path(dir_prefix);
    if (!boost::filesystem::exists(path)) {
        fprintf(stderr, "Creating output directory.\n\n");
        boost::filesystem::create_directory(dir_prefix);
    }
    path = boost::filesystem::canonical(dir_prefix);
    dir_prefix = path.generic_string();
    dir_prefix += "/";
    std::string vcf_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.vcf";
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string hap_vcf_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotypes.vcf";
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
    fprintf(stderr, "\nLoading input MAT files %s \n", input_mat_filename.c_str());

    // Load input MAT and uncondense tree
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));
    
    //Checking how close are input samples with peaks
    std::unordered_map<size_t, struct read_info*> read_map, hap_map;
    tbb::concurrent_hash_map<std::string, std::vector<size_t>> hap_read_map;
    std::vector<std::string> vcf_samples;
    std::unordered_map<std::string, double> hap_abun_map;
    
    readSampleVCF(vcf_samples, vcf_filename_samples);
    readVCF(hap_map, hap_vcf_filename, ref_seq.size(), false);
    computeDistance(T, hap_map, vcf_samples);

    ///////////////////////////////
    std::vector<std::string> clade_list = {
"B.1.1.28",
"B.1.429",
"P.1",
"B.1.177.50",
"B.11",
"B.1.160",
"B.1.426",
"R.1",
"B.1.160.7",
"B.1.177.7",
"B.1.243",
"B.61"
    };

    int hap_count = 0;
    std::vector<MAT::Node*> dfs = T.depth_first_expansion(); 
    for (const auto& curr_clade: clade_list) {
        for (auto n: dfs) {
            if (n->clade_annotations[1] == curr_clade) {
                std::queue<MAT::Node*> remaining_nodes;
                remaining_nodes.push(n);
                while (remaining_nodes.size() > 0) {
                    MAT::Node* curr_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    if ((curr_node->clade_annotations[1] == "") || (curr_node->clade_annotations[1] == curr_clade)) {
                        if (curr_node->children.size() == 0)
                            hap_count++;
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

    }
    fprintf(stderr, "HAP count = %d\n", hap_count);
}

//Computes distance between peaks and given samples
void computeDistance(const MAT::Tree &T, const std::unordered_map<size_t, struct read_info*> &hap_map, const std::vector<std::string> &vcf_samples) {
    fprintf(stderr, "Haplotypes: %d\n\n", (int)hap_map.size());
    printf("\nMUTATION DISTANCE NEW:\n");
    //Get Mutations of samples
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = getMutations(T, sample);
        std::string best_node = "";
        auto itr = hap_map.begin();
        while (itr != hap_map.end()) {
            int curr_dist = mutationDistance(sample_mutations, itr->second->mutations);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = itr->second->read;
            }
            itr++;
        }
        printf("Node: %s, Closest_node: %s, mutationDistance: %d\n", sample.c_str(), best_node.c_str(), min_dist);
    }

    //Samples' avg mut distance from other haplotypes of lineage
    printf("\n\nINTRA-LINEAGE AVG MUTATION DISTANCE:\n");
    //Depth first expansion to get all nodes in the tree
    std::vector<MAT::Node*> dfs = T.depth_first_expansion(); 
    for (auto sample: vcf_samples) {
        std::vector<MAT::Node*> hap_list;
        //Find all haplotypes of the lineage
        auto curr_clade = getLineage(T, T.get_node(sample));
        for (auto n: dfs) {
            if (n->clade_annotations[1] == curr_clade) {
                std::queue<MAT::Node*> remaining_nodes;
                remaining_nodes.push(n);
                while (remaining_nodes.size() > 0) {
                    MAT::Node* curr_node = remaining_nodes.front();
                    remaining_nodes.pop();
                    if ((curr_node->clade_annotations[1] == "") || (curr_node->clade_annotations[1] == curr_clade)) {
                        if (curr_node->children.size() == 0)
                            hap_list.emplace_back(curr_node);
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
        //Compute Intra-lineage mutation distance
        double avg_distance = 0.0;
        auto sample_mutations = getMutations(T, sample);
        for (auto hap: hap_list) {
            std::string hap_name = hap->identifier;
            if (hap_name != sample) {
                auto hap_mutations = getMutations(T, hap_name);
                avg_distance += mutationDistance(sample_mutations, hap_mutations);
            }
        }
        avg_distance /= (int)hap_list.size() - 1;
        hap_list.clear();
        printf("Node: %s, avg_mutation_distance: %f\n", sample.c_str(), avg_distance);
    }
}