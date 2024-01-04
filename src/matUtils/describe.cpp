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
        for (auto n2: root_to_sample) {
            for (size_t i=0; i<n2->mutations.size(); i++) {
                cpath += n2->mutations[i].get_string();
                if (i+1 < n2->mutations.size()) {
                    cpath += ",";
                }
            }
            if (n2->mutations.size() > 0 && n2 != root_to_sample.back()) {
                cpath += ">";
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

std::vector<std::string> clade_paths(MAT::Tree* T) {
    //get the set of clade path strings for printing
    //similar to the above, construct a series of strings to be printed or redirected later on
    std::vector<std::string> clpaths;
    clpaths.push_back("clade\troot_id\tfrom_tree_root\n");
    //do a breadth-first search
    //clades are annotated only at the root, so when we see an annotation, add it to the list.

    auto bfs = T->breadth_first_expansion();
    for (auto n: bfs) {
        for (auto ann: n->clade_annotations) {
            if (ann != "") {
                std::string curpath;
                //record the name of the clade
                curpath += ann + "\t";
                curpath += n->identifier + "\t";
                //get the ancestral mutations back to the root
                std::string root = "";
                auto root_to_sample = T->rsearch(n->identifier, true);
                std::reverse(root_to_sample.begin(), root_to_sample.end());
                for (auto an: root_to_sample) {
                    for (size_t i=0; i<an->mutations.size(); i++) {
                        root += an->mutations[i].get_string();
                        if (i+1 < an->mutations.size()) {
                            root += ",";
                        }
                    }
                    if (an != root_to_sample.back()) {
                        root += " > ";
                    }
                }
                //save values to the string and save the string
                curpath += root;
                clpaths.push_back(curpath + "\n");
            }
        }
    }
    return clpaths;
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
