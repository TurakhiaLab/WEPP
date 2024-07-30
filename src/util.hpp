#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <boost/program_options.hpp>

#include "panman/panmanUtils.hpp"

#include "mutation.hpp"
#include "read.hpp"

std::vector<std::string> 
read_sample_vcf(const std::string& vcf_filename_samples);

boost::program_options::variables_map 
parseWBEcommand(boost::program_options::parsed_options parsed);
