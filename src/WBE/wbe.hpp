#include "src/usher_graph.hpp"
#include <regex>
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace po = boost::program_options;

extern Timer timer;

// Define your mutex type
using my_mutex_t = tbb::queuing_mutex;

struct read_info {
   std::string read;
   std::vector<MAT::Mutation> mutations;
   int start;
   int end;
};

struct min_parsimony {
   std::vector<size_t> idx_list;
   std::vector<std::vector<MAT::Mutation>> par_list;
};

struct node_pair_clade {
   MAT::Node* ref_tree_node;
   MAT::Node* new_tree_parent;
   std::string ref_tree_parent_lineage;
};

struct node_branch {
   MAT::Node* node;
   int branch_length;
};

po::variables_map parseWBEcommand(po::parsed_options);

void sam2VCF(po::parsed_options);

void selectHaplotypes (po::parsed_options);

void filterLineages(po::parsed_options);

void detectPeaks (po::parsed_options);

void refinePeaks(po::parsed_options);

void readSampleVCF(std::vector<std::string> &, const std::string);

void readVCF(std::unordered_map<size_t, struct read_info*> &, const std::string &, const size_t &, const bool &);

void readCSV(std::unordered_map<std::string, double>& , const std::string &);

void readCSV(std::unordered_map<std::string, std::vector<std::string>>&, const std::string &);

void placeReadHelper(MAT::Node*, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::unordered_map<size_t, struct read_info*> &, std::vector<size_t>, const std::vector<MAT::Node*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, std::vector<size_t>&, const int &, const int &, const int &);

void placeReads(const MAT::Tree &, const std::vector<size_t> &, const std::unordered_map<size_t, struct read_info*> &, const std::vector<MAT::Node*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, std::vector<bool> &);

void updateParsimony(struct min_parsimony &, const std::vector<MAT::Mutation> &, const int &);

void analyzeReads(const MAT::Tree &, const MAT::Tree &, const std::string &, const std::unordered_map<size_t, struct read_info*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, const std::vector<std::string> &, const std::string &, const std::string &, const std::string &);

int mutationDistance(const MAT::Tree &, const MAT::Tree &, const MAT::Node*, const MAT::Node*);

int mutationDistance(std::vector<MAT::Mutation>, std::vector<MAT::Mutation>);

std::vector<MAT::Mutation> getMutations(const MAT::Tree&, const std::string);

std::string getLineage(const MAT::Tree &, MAT::Node*);

void updateProhibitedNodes(const MAT::Tree &, const std::vector<MAT::Node*> &, std::vector<MAT::Node*> &, const int&);

void getProhibitedNodes(const MAT::Tree &, const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Node* , std::vector<MAT::Node*> &, const int&);

std::vector<MAT::Node*> updateNeighborNodes(const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::vector<MAT::Node*> &, const std::vector<MAT::Node*> &, const tbb::concurrent_hash_map<MAT::Node*, double> &, std::vector<MAT::Node*> &, const int&, const int&);

void addNeighborNodes(const MAT::Tree &, std::vector<MAT::Node*> &, const int&, const int&);

int childNodesAddition(const MAT::Tree &, my_mutex_t &, const std::vector<MAT::Node*> &, std::vector<MAT::Node*> &, int &, std::vector<MAT::Node*> &, MAT::Node *, MAT::Node *, const int &, const int &, const int &, const int &);

void generateFilteringData(const MAT::Tree &, const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::string &, const std::vector<MAT::Node*> &, const std::unordered_map<size_t, struct read_info*> &, const std::string &, const std::string &, const std::string &);

bool compareMutations(const MAT::Mutation &, const MAT::Mutation &);

void sortNodeScore(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const tbb::concurrent_hash_map<MAT::Node*, double> &, tbb::concurrent_vector<std::pair<MAT::Node*, double>> &);

bool compareNodeScore(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::pair<MAT::Node*, double>&, const std::pair<MAT::Node*, double>&);

bool compareReadPos (const std::pair<size_t, std::string>&, const std::pair<size_t, std::string>&);

size_t getNumLeaves(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Node*);

void createRangeTree(MAT::Node*, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const int &, const int &, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Tree &);

void createLineageTree(MAT::Node*, const std::vector<std::string> &, MAT::Tree &);

void createCondensedTree(MAT::Node*, const std::unordered_map<size_t, struct read_info*> &, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Tree &);

void computeDistance(const MAT::Tree &, const std::unordered_map<size_t, struct read_info*> &, const std::vector<std::string> &, const std::unordered_map<std::string, double> &, const std::unordered_set<int> &);

void placeReads(const MAT::Tree &, const std::string &, const std::unordered_map<size_t, struct read_info*> &, const std::unordered_map<size_t, struct read_info*> &);

void readSAM(const std::string &, const std::string &, std::unordered_map<size_t, std::string> &, std::unordered_map<int, std::vector<std::tuple<std::string, std::string, std::vector<size_t>>>> &);

bool compareIdx(const std::pair<int, size_t> &, const std::pair<int, size_t> &);