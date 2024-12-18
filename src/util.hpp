#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <boost/program_options.hpp>

#include "panman/panmanUtils.hpp"

#include "panman_bridge.hpp"
#include "mutation.hpp"
#include "read.hpp"

std::vector<std::string> 
read_sample_vcf(const std::string& vcf_filename_samples);

boost::program_options::variables_map 
parseWBEcommand(boost::program_options::parsed_options parsed);

void 
get_single_mutations(std::vector<mutation>& mutations, const std::string& ref, const panmanUtils::Node* node, const coord_converter &coord);

float mutation_distance(std::vector<mutation> const& node1_mutations, std::vector<mutation> const& node2_mutations);

std::vector<mutation> mutation_distance_vector(std::vector<mutation> const& node1_mutations, std::vector<mutation> const& node2_mutations);