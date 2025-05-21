#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <boost/program_options.hpp>

#include "src/usher_graph.hpp"

#include "read.hpp"

std::vector<std::pair<std::string, std::vector<MAT::Mutation>>>
read_sample_vcf(const std::string& vcf_filename_samples);

MAT::Tree 
create_condensed_tree(MAT::Node* ref_root, const std::unordered_set<int>&site_read_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings);

extern boost::program_options::options_description conv_desc;
void initializeConvDesc();

boost::program_options::variables_map parseWEPPcommand(boost::program_options::parsed_options parsed);

int mutation_distance(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations);

std::vector<MAT::Mutation> mutation_distance_vector(std::vector<MAT::Mutation> node1_mutations, std::vector<MAT::Mutation> node2_mutations, const std::string& reference);

std::vector<MAT::Mutation> get_mutations(const MAT::Tree& T, const std::string sample);

size_t get_num_leaves(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &condensed_node_mappings, MAT::Node* condensed_node);