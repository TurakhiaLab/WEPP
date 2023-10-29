#pragma once

#include "convert.hpp"
#include "select.hpp"
#include <string>
#include <regex>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>

struct read_info {
   std::string read;
   std::vector<MAT::Mutation> mutations;
   int start;
   int end;
};

struct parsimony {
   std::vector<MAT::Mutation> p_node_par;
   MAT::Node* curr_node; 
};

struct min_parsimony {
   std::vector<int> idx_list;
   std::vector<std::vector<MAT::Mutation>> par_list;
};

struct node_branch {
   MAT::Node* node;
   int branch_length;
};

po::variables_map parse_peak_search_command(po::parsed_options parsed);

void peak_search(po::parsed_options parsed);

void read_sample_vcf(std::vector<std::string> &, const std::string);

void read_vcf(std::unordered_map<size_t, struct read_info*> &, const std::string);

int place_reads(const std::vector<MAT::Node*> &, struct read_info*, const MAT::Node*, tbb::concurrent_hash_map<MAT::Node*, double> &);

void analyze_reads(const MAT::Tree &, const MAT::Tree &, const std::unordered_map<size_t, struct read_info*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, const std::vector<std::string> &, const std::string &, const std::string &);

int mutation_distance(const MAT::Tree &, const MAT::Tree &, const MAT::Node*, const MAT::Node*);

std::string get_clade(const MAT::Tree &, MAT::Node*);

void add_neighbor_peaks(const MAT::Tree &, std::vector<MAT::Node*> &, const int, const int);

void child_nodes_addition(my_mutex_t &, const std::vector<MAT::Node*> &, std::vector<MAT::Node*> &, int &, std::vector<MAT::Node*> &, const MAT::Node*, const MAT::Node*, const int &, const int &, const int &, const int &);

void generate_regression_abundance_data(const MAT::Tree &, const std::vector<MAT::Node*> &, const std::unordered_map<size_t, struct read_info*> &, const std::string &, const std::string &);

bool compare_mutations(const MAT::Mutation &, const MAT::Mutation &);