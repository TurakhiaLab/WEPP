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
    fprintf(stderr, "Nodes: %ld\n", dfs.size());

    ////////////////////////////////////////////////////GIVES full path
    std::ifstream file("aa_2023_12_15_node_mutations.tsv");
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return mpaths;
    }

    std::unordered_map<std::string, std::string> node_mutations;
    std::string line, node_id, aa_mutations, nt_mutations, codon_changes, leaves_sharing_mutations;

    // Read the header line (if needed)
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        
        std::getline(ss, node_id, '\t');
        std::getline(ss, aa_mutations, '\t');
        std::getline(ss, nt_mutations, '\t');
        std::getline(ss, codon_changes, '\t');
        std::getline(ss, leaves_sharing_mutations, '\t');

        //Check for redundant mutation by separating aa_mutations based on ";"
        std::string filtered_mutations;
        std::istringstream aa_stream(aa_mutations);
        std::string mutation;

        // Iterate over each part and split further based on ":"
        while (std::getline(aa_stream, mutation, ';')) {
            // Split the mutation based on `:`
            std::istringstream mutation_stream(mutation);
            std::string part1, part2;

            std::getline(mutation_stream, part1, ':');
            std::getline(mutation_stream, part2);
            
            // Check if the first and last characters of `part2` are the same
            if ((!part2.empty()) && (part2.front() != part2.back())) {
                // If they are not the same, include this mutation
                if (!filtered_mutations.empty()) {
                    filtered_mutations += ";";
                }
                filtered_mutations += mutation;
            }
        }

        if (!filtered_mutations.empty())
            node_mutations[node_id] = filtered_mutations;
    }
    file.close();

    for (const auto& curr_node: dfs) {
        if (curr_node->is_leaf()) {
            bool mut_found = false;
            std::string cpath = curr_node->identifier + "\t";
            auto path_nodes = T->rsearch(curr_node->identifier, true);
            std::reverse(path_nodes.begin(), path_nodes.end());
            for (const auto& curr_node: path_nodes) {
                if (node_mutations.find(curr_node->identifier) != node_mutations.end()) {
                    if (mut_found)
                        cpath += ">";
                    cpath += node_mutations[curr_node->identifier];
                    mut_found = true;
                }
            }
            if (mut_found)
                mpaths.emplace_back(cpath);
        }
    }
    
    
    


    /*
    //////////////////////////////////////////////////REMOVE
    // Get full path nucleotide mutations
    for (const auto& curr_node: dfs) {
        if (curr_node->is_leaf()) {
            bool mut_found = false;
            std::string cpath = curr_node->identifier + "\t";
            auto path_nodes = T->rsearch(curr_node->identifier, true);
            std::reverse(path_nodes.begin(), path_nodes.end());
            for (const auto& curr_node: path_nodes) {
                for (const auto& mut: curr_node->mutations) {
                    if (mut.get_string() != curr_node->mutations.front().get_string())
                        cpath += ",";
                    else if (mut_found)
                        cpath += ">";
                   cpath += mut.get_string();
                   mut_found = true; 
                }
            }
            mpaths.emplace_back(cpath);
        }
    }
    //////////////////////////////////////////////////
    */

    //for (const auto& curr_node: dfs) {
    //    if (curr_node->is_leaf()) {
    //        std::string cpath = curr_node->identifier + "\t";
    //        std::vector<MAT::Mutation> mut_list;
    //        for (auto anc: T->rsearch(curr_node->identifier, true)) {
    //            for (auto c_mut: anc->mutations) {
    //                //Only consider the mutation closest to root at a site
    //                auto m_itr = mut_list.begin();
    //                while (m_itr != mut_list.end()) {
    //                    if (m_itr->position == c_mut.position)
    //                        break;
    //                    m_itr++;
    //                }
    //                if (m_itr == mut_list.end())
    //                    mut_list.emplace_back(c_mut);
    //            }
    //        }
    //        //Remove mutations where ref_nuc is same as mut_nuc -> Back Mutation
    //        auto m_itr = mut_list.begin();
    //        while (m_itr != mut_list.end()) {
    //            if (m_itr->mut_nuc == m_itr->ref_nuc)
    //                m_itr = mut_list.erase(m_itr);
    //            else
    //                m_itr++;
    //        }

    //        //Reversing mutations to get from root to leaf    
    //        std::vector<MAT::Node*> ancestor_list = T->rsearch(curr_node->identifier, true);
    //        MAT::Node* prev_node = NULL;
    //        int j = ancestor_list.size() - 1;
    //        for (int i = mut_list.size() - 1; i >= 0; --i) {
    //            auto curr_mut = mut_list[i];
    //            bool found = false;
    //            while (j >= 0) {
    //                MAT::Node* curr_node = ancestor_list[j];
    //                for (const auto& mut: curr_node->mutations) {
    //                    if ((mut.position == curr_mut.position) && (mut.mut_nuc == curr_mut.mut_nuc)) {
    //                        found = true;
    //                        break;
    //                    }
    //                }
    //                if (found) {
    //                    if (prev_node == curr_node) {
    //                        cpath += ",";
    //                        cpath += MAT::get_nuc(curr_mut.ref_nuc) + std::to_string(curr_mut.position) + MAT::get_nuc(curr_mut.mut_nuc);
    //                       
    //                    }
    //                    else {
    //                        if (prev_node != NULL)
    //                            cpath += ">";
    //                        cpath += MAT::get_nuc(curr_mut.ref_nuc) + std::to_string(curr_mut.position) + MAT::get_nuc(curr_mut.mut_nuc);
    //                     
    //                    }
    //                    prev_node = curr_node;
    //                    break;
    //                }
    //                j--;
    //            }
    //        }
    //        mpaths.emplace_back(cpath);
    //    }
    //}
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
