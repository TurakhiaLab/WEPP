#include "src/matUtils/experiment.hpp"

void compute_distance(const MAT::Tree &T, const std::unordered_map<int, struct read_info*> &read_map, const std::vector<std::string> &vcf_samples) {
    fprintf(stderr, "Haplotypes: %d\n\n", (int)read_map.size());
    //Get Mutations of samples
    for (auto sample: vcf_samples) {
        int min_dist = 100000;
        std::vector<MAT::Mutation> sample_mutations;
        for (auto anc: T.rsearch(sample, true)) { //Checking all ancestors of a node
            for (auto mut: anc->mutations) {
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
        //Remove Back-Mutations
        auto sm_itr = sample_mutations.begin();
        while (sm_itr != sample_mutations.end()) {
            if (sm_itr->ref_nuc == sm_itr->mut_nuc)
                sm_itr = sample_mutations.erase(sm_itr);
            else
                sm_itr++;
        }
        std::string best_node = "";
        auto itr = read_map.begin();
        while (itr != read_map.end()) {
            int curr_dist = mutation_distance(sample_mutations, itr->second->mutations);
            if (curr_dist < min_dist) {
                min_dist = curr_dist;
                best_node = itr->second->read;
            }
            itr++;
        }
        printf("Node: %s, Closest_node: %s, mutation_distance: %d\n", sample.c_str(), best_node.c_str(), min_dist);
    }
}

void simulate_and_place_reads (po::parsed_options parsed) {
    // Load input MAT and uncondense tree
    MAT::Tree T;
    T = MAT::load_mutation_annotated_tree("public-2021-05-31.all.masked.nextclade.pangolin.pb");
    T.uncondense_leaves();

    //Depth first expansion to get all nodes in the tree and 
    // comparison with given lineage to get all nodes of the required lineage 
    std::vector<MAT::Node*> dfs = T.depth_first_expansion(T.root); 
    
    //Checking how close are input samples with peaks
    std::unordered_map<int, struct read_info*> read_map;
    std::vector<std::string> vcf_samples;
    read_sample_vcf(vcf_samples, "my_vcf_samples.vcf");
    read_vcf(T, dfs, read_map, "my_vcf_haplotypes.vcf");
    compute_distance(T, read_map, vcf_samples);
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

void read_vcf(const MAT::Tree &T, const std::vector<MAT::Node*> &dfs, std::unordered_map<int, struct read_info*> &read_map, const std::string vcf_filename_reads) {
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

//Function to calculation distance between two nodes
int mutation_distance(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations) {
    std::sort(node1_mutations.begin(), node1_mutations.end(), compare_mutations);
    std::sort(node2_mutations.begin(), node2_mutations.end(), compare_mutations);
    auto n1_itr = node1_mutations.begin();
    while (n1_itr != node1_mutations.end()) {
        bool found_both = false;
        auto n2_itr = node2_mutations.begin();
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

//Comparing mutations for sorting a vector 
bool compare_mutations(const MAT::Mutation &a, const MAT::Mutation &b) {
    return a.position < b.position;
}