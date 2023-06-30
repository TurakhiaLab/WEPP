#include "experiment.hpp"

po::variables_map parse_post_processing_command(po::parsed_options);

int mutation_distance(std::vector<MAT::Mutation>, std::vector<MAT::Mutation>);

void compute_distance(const MAT::Tree &, const std::unordered_map<int, struct read_info*> &, const std::vector<std::string> &);

void post_processing(po::parsed_options);


