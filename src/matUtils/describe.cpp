#include "describe.hpp"
#include "translate.hpp"

std::vector<std::string> mutation_paths(MAT::Tree* T, std::vector<std::string> samples) {
    std::vector<std::string> mpaths;
    mpaths.push_back("sample_id\tpath_from_root");
    for (auto sample: samples) {
        std::string cpath = sample + "\t";
        auto root_to_sample = T->rsearch(sample, true);
        std::reverse(root_to_sample.begin(), root_to_sample.end());
        for (auto n: root_to_sample) {
            for (size_t i=0; i<n->mutations.size(); i++) {
                std::string mutation_string = n->mutations[i].get_string();
                cpath += mutation_string;
                if (i+1 < n->mutations.size()) {
                    cpath += ",";
                }
            }
            if (n != root_to_sample.back()) {
                //note, for uncondensed samples with parsimony score 0,
                //this will leave a > at the end. That's basically fine.
                //cpath += " (" + n->identifier + ") > ";
                cpath += ";";
            }
        }
        mpaths.push_back(cpath);
    }
    return mpaths;
}

std::vector<std::string> mutation_paths_all(MAT:: Tree*T) {
    std::vector<std::string> mpaths;
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        std::string cpath = n->identifier + "\t";
        auto root_to_sample = T->rsearch(n->identifier, true);
        std::reverse(root_to_sample.begin(), root_to_sample.end());
        bool first_mut_found = false;
        for (size_t j = 0; j < root_to_sample.size(); j++) {
            auto n2 = root_to_sample[j];
            //Handling special characters
            if (!n2->mutations.empty()) {
                if (!first_mut_found)
                    first_mut_found = true;
                else 
                    cpath += ">";
            }

            //Iterating through n2 mutations
            for (size_t i=0; i < n2->mutations.size(); i++) {
                cpath += n2->mutations[i].get_string();
                if (i+1 < n2->mutations.size()) {
                    cpath += ",";
                }
            }
        }
        mpaths.push_back(cpath);
    }
    return mpaths;
}

std::vector<std::string> mutation_paths_lineages(MAT::Tree* T, std::vector<std::string> lineages) {
    std::vector<std::string> mpaths;
    std::vector<MAT::Node*> dfs = T->depth_first_expansion(); 
    for (auto lineage: lineages) {
        for (auto curr_node: dfs) {
            if (curr_node->clade_annotations[1] == lineage) {
                std::vector<MAT::Mutation> mut_list;
                for (auto anc: T->rsearch(curr_node->identifier, true)) {
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
                std::string cpath = lineage + "\t";

                //Reversing mutations to get from root to leaf    
                std::vector<MAT::Node*> ancestor_list = T->rsearch(curr_node->identifier, true);

                MAT::Node* prev_node = NULL;
                int j = ancestor_list.size() - 1;
                for (int i = mut_list.size() - 1; i >= 0; --i) {
                    auto curr_mut = mut_list[i];
                    bool found = false;
                    while (j >= 0) {
                        MAT::Node* curr_node = ancestor_list[j];
                        for (const auto& mut: curr_node->mutations) {
                            if ((mut.position == curr_mut.position) && (mut.mut_nuc == curr_mut.mut_nuc)) {
                                found = true;
                                break;
                            }
                        }
                        if (found) {
                            if (prev_node == curr_node) {
                                cpath += ",";
                                cpath += MAT::get_nuc(curr_mut.ref_nuc) + std::to_string(curr_mut.position) + MAT::get_nuc(curr_mut.mut_nuc);
                               
                            }
                            else {
                                cpath += ">";
                                cpath += MAT::get_nuc(curr_mut.ref_nuc) + std::to_string(curr_mut.position) + MAT::get_nuc(curr_mut.mut_nuc);
                             
                            }
                            prev_node = curr_node;
                            break;
                        }
                        j--;
                    }
                }
                mpaths.push_back(cpath);

                mut_list.clear();
                break;
            }
        }
    }
    return mpaths;
}

std::vector<std::string> all_nodes_paths(MAT::Tree* T) {
    std::vector<std::string> dfs_strings;
    auto dfs = T->depth_first_expansion();
    for (auto n: dfs) {
        std::string node_path = n->identifier + ": ";
        for (size_t i=0; i<n->mutations.size(); i++) {
            node_path += n->mutations[i].get_string();
            if (i+1 < n->mutations.size()) {
                node_path += ",";
            }
        }
        dfs_strings.push_back(node_path);
    }
    return dfs_strings;
}
