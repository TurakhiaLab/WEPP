#include <iostream>
#include <boost/filesystem.hpp>
#include <unordered_set>

#include "util.hpp"

//Get the names of samples prsent in samples.vcf
std::vector<std::string> read_sample_vcf(const std::string vcf_filename_samples) {
    std::vector<std::string> vcf_samples;

    // Boost library used to stream the contents of the input VCF file
    boost::filesystem::ifstream fileHandler(vcf_filename_samples);
    std::string s;
    while (getline(fileHandler, s)) {
        std::vector<std::string> words;
        MAT::string_split(s, words);
        if (words.size() > 1) {
            //Checking for header
            if (words[1] == "POS") {
                //Leave certain fields based on our VCF format
                for (int j=9; j < (int)words.size(); j++)
                    vcf_samples.emplace_back(words[j]);
            }
        }
    }

    return vcf_samples;
}

void create_condensed_tree(MAT::Node* ref_root, const std::vector<read> &read_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings, MAT::Tree &T) {
    //MAP of site_read
    std::unordered_set<int> site_read_map;
    for (size_t i = 0; i < read_map.size(); i++) {
        auto rp = read_map[i].second;
        for (int j = rp->start; j <= rp->end; j++)
            site_read_map.insert(j);
    }

    //REMOVE sites not covered by reads
    std::queue<std::pair<MAT::Node*, MAT::Node*>> remaining_nodes;
    auto new_node = T.create_node("DUMMY-CONDENSED", -1.0, ref_root->clade_annotations.size());
    node_mappings[new_node] = std::vector<MAT::Node*>();
    remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(ref_root, new_node));
    
    //Add the new_node to the node_mappings
    while(remaining_nodes.size() > 0) {
        auto r_curr_node = remaining_nodes.front().first;
        auto n_parent_node = remaining_nodes.front().second;
        remaining_nodes.pop();
        std::vector<MAT::Mutation> covered_mutations;
        for (const auto& mut: r_curr_node->mutations) {
            if (site_read_map.find(mut.position) != site_read_map.end())
                covered_mutations.emplace_back(mut);
        }

        //Add a new node to tree if mutation found in range
        if (!covered_mutations.empty()) {
            //Create a new_node
            auto new_node = T.create_node(r_curr_node->identifier, n_parent_node, -1.0);
            //Add mutations to new_node
            new_node->mutations.reserve(covered_mutations.size());
            new_node->mutations.insert(new_node->mutations.end(), covered_mutations.begin(), covered_mutations.end());
            covered_mutations.clear();
            //Add new_node to the node_mappings
            node_mappings[new_node] = {r_curr_node};
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, new_node));
        }
        //Without any mutation, can only be Placed as a child if r_curr_node is NOT leaf node
        else {
            //Add current_node to the n_parent_node's list in node_mappings
            node_mappings[n_parent_node].emplace_back(r_curr_node);
            //Add children to remaining_nodes    
            for (auto child: r_curr_node->children)
                remaining_nodes.push(std::pair<MAT::Node*, MAT::Node*>(child, n_parent_node));
        }
    }
}
