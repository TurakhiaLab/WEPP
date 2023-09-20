#include "experiment.hpp"

po::variables_map parse_post_processing_command(po::parsed_options);

void post_processing(po::parsed_options);

void read_csv(std::unordered_map<std::string, double> &, const std::string);

void read_csv(std::unordered_map<std::string, std::string> &, const std::string);

int mutation_distance(std::vector<MAT::Mutation>, std::vector<MAT::Mutation>);

void compute_distance(const MAT::Tree &, const std::unordered_map<size_t, struct read_info*> &, const std::vector<std::string> &);

void compute_abundance(const std::unordered_map<std::string, double> &, const std::unordered_map<std::string, std::string> &);

std::vector<MAT::Mutation> get_mutations(const MAT::Tree&, const std::string);

void place_reads(const std::unordered_map<int, struct read_info*>&, const std::unordered_map<int, struct read_info*>&, const std::unordered_map<std::string, double>&, tbb::concurrent_hash_map<std::string, std::vector<size_t>>&);