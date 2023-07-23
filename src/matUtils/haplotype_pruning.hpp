#include "filter.hpp"
#include "select.hpp"
#include "experiment.hpp"

po::variables_map parse_haplotype_pruning_command(po::parsed_options);

void haplotype_pruning(po::parsed_options);

std::vector<std::string> samples_outside_mut_dist(const MAT::Tree&, std::vector<std::string>, int);
