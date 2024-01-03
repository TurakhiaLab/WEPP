#pragma once
#include "experiment.hpp"
#include <string>
#include <random>
#include <regex>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>

po::variables_map parse_post_processing_command(po::parsed_options parsed);

void read_csv(std::unordered_map<std::string, double> &, const std::string);

void read_csv(std::unordered_map<std::string, std::string> &, const std::string);

int mutation_distance(std::vector<MAT::Mutation>, std::vector<MAT::Mutation>);

void compute_distance(const MAT::Tree &, const std::vector<MAT::Node*> &, const std::unordered_map<int, struct read_info*> &, const std::vector<std::string> &);

void post_processing(po::parsed_options parsed);

void compute_abundance(const MAT::Tree &, const std::unordered_map<int, struct read_info*> &, const std::unordered_map<std::string, double> &, const std::unordered_map<std::string, std::string> &);

std::vector<MAT::Mutation> get_mutations(const MAT::Tree&, const std::string);
