#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <boost/program_options.hpp>

#include "panman/panman.hpp"

#include "mutation.hpp"
#include "read.hpp"

std::vector<std::string> 
read_sample_vcf(const std::string& vcf_filename_samples);

panmanUtils::Tree 
create_condensed_tree(panmanUtils::Node* ref_root, const std::vector<raw_read> &read_map, std::unordered_map<panmanUtils::Node*, std::vector<panmanUtils::Node*>> &node_mappings);

boost::program_options::variables_map 
parseWBEcommand(boost::program_options::parsed_options parsed);

int 
mutation_distance(std::vector<mutation> node1_mutations, std::vector<mutation> node2_mutations);

std::vector<mutation> 
get_mutations(const panmanUtils::Tree& T, const std::string sample);