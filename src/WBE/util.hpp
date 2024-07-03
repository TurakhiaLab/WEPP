#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <boost/program_options.hpp>

#include "src/usher_graph.hpp"

#include "read.hpp"

std::vector<std::string> 
read_sample_vcf(const std::string& vcf_filename_samples);

MAT::Tree 
create_condensed_tree(MAT::Node* ref_root, const std::vector<raw_read> &read_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings);

boost::program_options::variables_map 
parseWBEcommand(boost::program_options::parsed_options parsed);

int 
mutation_distance(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations);

std::vector<MAT::Mutation> 
get_mutations(const MAT::Tree& T, const std::string sample);

size_t 
get_num_leaves(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, MAT::Node* condensed_node);