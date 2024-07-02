#pragma once

#include <vector>
#include <string>
#include <unordered_map>

#include "src/usher_graph.hpp"

#include "read.hpp"

std::vector<std::string> read_sample_vcf(const std::string& vcf_filename_samples);
void create_condensed_tree(MAT::Node* ref_root, const std::vector<read> &read_map, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &node_mappings, MAT::Tree &T);
