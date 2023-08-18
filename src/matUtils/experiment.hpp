#pragma once

#include "convert.hpp"
#include "select.hpp"
#include <string>
#include <random>
#include <regex>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>

struct ances_sample_list {
   Mutation_Annotated_Tree::Node * ancestor_node; 
   std::vector<Mutation_Annotated_Tree::Node*> * sample_nodes;
};

struct sample_ances_list {
   Mutation_Annotated_Tree::Node * sample; 
   std::vector<Mutation_Annotated_Tree::Node*> * ances;
};

struct sample_read_pair {
   Mutation_Annotated_Tree::Node * sample;
   std::string read;
};

struct pos_misread {
   int pos;
   std::string read;
   bool used;
};

struct read_mut_pair {
   std::string read_name;
   char mut_nuc;
};

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

po::variables_map parse_extract_command(po::parsed_options parsed);
void simulate_and_place_reads(po::parsed_options parsed);

void read_sample_vcf(std::vector<std::string> &, const std::string);

void read_vcf(std::unordered_map<int, struct read_info*> &, const std::string);

int place_reads(const MAT::Tree &, const std::vector<MAT::Node*> &, struct read_info*, const MAT::Node*, tbb::concurrent_hash_map<MAT::Node*, double> &, const std::vector<MAT::Node*> &, std::vector<int> &, const int);

void analyze_reads(const MAT::Tree &, const MAT::Tree &, const std::vector<MAT::Node*> &, const std::unordered_map<int, struct read_info*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, const std::vector<std::string> &, const std::string &, const std::string &, const std::string &);

bool check_peaks_neighbourhood (const MAT::Tree &, const MAT::Node*, const std::vector<MAT::Node*> &, const int);
int mutation_distance(const MAT::Tree &, const MAT::Tree &, const MAT::Node*, const MAT::Node*);
std::string get_clade(const MAT::Tree &, MAT::Node*);

void add_neighbor_peaks(const MAT::Tree &, std::vector<MAT::Node*> &, const int, const int);

void generate_EM_data(const MAT::Tree &, const std::vector<MAT::Node*> &, const std::unordered_map<int, struct read_info*> &, const std::vector<MAT::Node*> &, const std::string &);
void generate_regression_abundance_data(const MAT::Tree &, const std::vector<MAT::Node*> &, const std::unordered_map<int, struct read_info*> &, const std::string &, const std::string &);
bool compare_mutations(const MAT::Mutation &, const MAT::Mutation &);