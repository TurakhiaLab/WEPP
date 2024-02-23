#include "wbe.hpp"

void refinePeaks(po::parsed_options parsed) {
    //main argument for the complex extract command
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
    std::string vcf_filename_samples = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_samples.vcf";
    std::string vcf_filename_reads = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_reads.vcf";
    std::string hap_vcf_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotypes.vcf";
    std::string hap_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_haplotype_abundance.csv";
    std::string freyja_lineage_csv_filename = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_freyja_results.csv";
    std::string condensed_nodes_csv = dir_prefix + vm["output-files-prefix"].as<std::string>() + "_condensed_nodes.csv";
    std::string ref_fasta = dir_prefix + vm["ref-fasta"].as<std::string>();
    uint32_t num_threads = vm["threads"].as<uint32_t>();
    fprintf(stderr, "\nNum Cores: %d\n\n", num_threads);
    tbb::task_scheduler_init init(num_threads);
    
    //Loading reference genome
    timer.Start();
    std::ifstream fasta_f(ref_fasta);
    if (!fasta_f.is_open()) {
        std::cerr << "Error: Unable to open file " << ref_fasta << std::endl;
        exit(1); 
    }
    std::string ref_header;
    std::getline(fasta_f, ref_header);
    std::string temp;
    std::string ref_seq;
    while (fasta_f) {
        std::getline(fasta_f, temp);
        ref_seq += temp;
    }

    // Load input MAT and uncondense tree
    fprintf(stderr, "\nLoading input MAT files %s \n", input_mat_filename.c_str());
    MAT::Tree T = MAT::load_mutation_annotated_tree(input_mat_filename);
    T.uncondense_leaves();
    fprintf(stderr, "Completed in %ld sec \n\n", (timer.Stop() / 1000));

    //Checking how close are input samples with peaks
    std::unordered_map<size_t, struct read_info*> read_map, hap_map;
    tbb::concurrent_hash_map<std::string, std::vector<size_t>> hap_read_map;
    std::vector<std::string> vcf_samples;
    std::unordered_map<std::string, double> hap_abun_map, freyja_lineage_abun_map;
    std::vector<MAT::Node*> peak_nodes;
    std::unordered_map<std::string, std::vector<std::string>> condensed_nodeNames_map;
    
    //Compute mutation distance by reading files returned by python code
    readSampleVCF(vcf_samples, vcf_filename_samples);
    readVCF(hap_map, hap_vcf_filename, ref_seq.size());
    //readCSV(hap_abun_map, hap_csv_filename);
    readCSV(condensed_nodeNames_map, condensed_nodes_csv);
    readCSV(freyja_lineage_abun_map, freyja_lineage_csv_filename);
    readVCF(read_map, vcf_filename_reads, ref_seq.size());

    //UPDATE hap_map if CONDENSED nodes are found
    std::vector<struct read_info*> uncondensed_nodes;
    for (size_t i = 0; i < hap_map.size(); i++) {
        auto hm = hap_map[i];
        if (hm->read.find("CONDENSED") != std::string::npos) {
            auto node_names_list = condensed_nodeNames_map[hm->read];
            for (size_t i = 0 ; i < node_names_list.size(); i++) {
                auto node_name = node_names_list[i];
                if (!i) {
                    hm->read = node_name;
                }
                else {
                    struct read_info* rp = new struct read_info;
                    rp->read = node_name;
                    rp->mutations = hm->mutations;
                    uncondensed_nodes.emplace_back(rp);
                }
            }
        }
    }
    //Insert uncondensed nodes in hap_map
    for (const auto& hap: uncondensed_nodes)
        hap_map.insert({hap_map.size(), hap});
    uncondensed_nodes.clear();

    //MAP of sites covered by reads
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map.find(i)->second;
        for (int j = rp->start; j <= rp->end; j++)
            site_read_map.insert(j);
    }
    
    computeDistance(T, hap_map, vcf_samples, freyja_lineage_abun_map, site_read_map);
    //placeReads(T, ref_seq, read_map, hap_map);
}

//Computes distance between peaks and given samples
void computeDistance(const MAT::Tree &T, const std::unordered_map<size_t, struct read_info*> &hap_map, const std::vector<std::string> &vcf_samples, const std::unordered_map<std::string, double> &freyja_lineage_abun_map, const std::unordered_set<int> &site_read_map) {
    printf("Haplotypes: %d\n\n", (int)hap_map.size());
    printf("\nMUTATION DISTANCE NEW:\n");
    //Closest distance of sample from peak
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = getMutations(T, sample);
        //Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end()) {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        std::string best_node = "";
        std::vector<MAT::Mutation> best_hap_mutations;
        auto hap_itr = hap_map.begin();
        while (hap_itr != hap_map.end()) {
            size_t last_underscore = hap_itr->second->read.find_last_of('_');
            std::string curr_hap_name = hap_itr->second->read.substr(0, last_underscore);
            while (T.get_node(curr_hap_name) == NULL) {
                last_underscore = curr_hap_name.find_last_of('_');
                curr_hap_name = curr_hap_name.substr(0, last_underscore);
            }
            //Getting hap_mutations from the Tree
            auto hap_mutations = getMutations(T, curr_hap_name);
            //Remove mutations from hap_mutations that are not present in site_read_map
            auto mut_itr = hap_mutations.begin();
            while (mut_itr != hap_mutations.end()) {
                if (site_read_map.find(mut_itr->position) == site_read_map.end())
                    mut_itr = hap_mutations.erase(mut_itr);
                else
                    mut_itr++;
            }

            int curr_dist = mutationDistance(sample_mutations, hap_mutations);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = hap_itr->second->read;
                best_hap_mutations = hap_mutations;
            }
            hap_itr++;
        }
        printf("\nNode: %s, Closest_node: %s, mutationDistance: %d", sample.c_str(), best_node.c_str(), min_dist);
    
    }

    //Closest distance of sample from root node of lineage
    printf("\n\nFREYJA AVG MUTATION DISTANCE:\n");
    //Depth first expansion to get all nodes in the tree
    std::vector<MAT::Node*> dfs = T.depth_first_expansion(); 
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        auto sample_mutations = getMutations(T, sample);
        //Remove mutations from sample_mutations that are not present in site_read_map
        auto mut_itr = sample_mutations.begin();
        while (mut_itr != sample_mutations.end()) {
            if (site_read_map.find(mut_itr->position) == site_read_map.end())
                mut_itr = sample_mutations.erase(mut_itr);
            else
                mut_itr++;
        }

        std::string best_node = "";
        for (const auto& lin_abun: freyja_lineage_abun_map) {
            auto curr_lin = lin_abun.first;
            //Find root haplotype of the lineage
            for (const auto& curr_node: dfs) {
                if (curr_node->clade_annotations[1] == curr_lin) {
                    auto hap_mutations = getMutations(T, curr_node->identifier);
                    //Remove mutations from hap_mutations that are not present in site_read_map
                    auto mut_itr = hap_mutations.begin();
                    while (mut_itr != hap_mutations.end()) {
                        if (site_read_map.find(mut_itr->position) == site_read_map.end())
                            mut_itr = hap_mutations.erase(mut_itr);
                        else
                            mut_itr++;
                    }

                    int curr_dist = mutationDistance(hap_mutations, sample_mutations);
                    if (curr_dist < min_dist) {
                        min_dist = curr_dist;
                        best_node = curr_node->identifier;
                    }
                    break;
                }
            }
        }
        printf("\nNode: %s, Closest_node: %s, mutationDistance: %d", sample.c_str(), best_node.c_str(), min_dist);
    }
}
