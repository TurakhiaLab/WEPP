#include "src/usher_graph.hpp"
#include <regex>
#include <array>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>
#include <random>
#include <string>
#include <map>
#include <optional>
#include <utility>
#include <ostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include "sam.pb.h"  

namespace po = boost::program_options;

extern Timer timer;

// Define your mutex type
using my_mutex_t = tbb::queuing_mutex;

struct read_info {
   std::string read;
   std::vector<MAT::Mutation> mutations;
   int start;
   int end;
   int degree;
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

void load_reads_from_proto(std::string const &, std::unordered_map<size_t, struct read_info *> &, std::unordered_map<std::string, std::vector<std::string>> &);

void sam2PB(po::parsed_options);

void selectHaplotypes (po::parsed_options);

void filterLineages(po::parsed_options);

void detectPeaks (po::parsed_options);

void refinePeaks(po::parsed_options);

void readSampleVCF(std::vector<std::string> &, const std::string);

void readVCF(std::unordered_map<size_t, struct read_info*> &, const std::string &, const size_t &);

void readCSV(std::unordered_map<std::string, double>& , const std::string &);

void readCSV(std::unordered_map<std::string, std::vector<std::string>>&, const std::string &);

std::vector<int> getRandomElements(const int &, const int &); 

void placeReadHelper(MAT::Node*, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::unordered_map<size_t, struct read_info*> &, std::vector<size_t>, const std::vector<MAT::Node*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, std::vector<size_t>&, const int &, const int &, const int &, const bool &);

void placeReads(const MAT::Tree &, const std::vector<size_t> &, const std::unordered_map<size_t, struct read_info*> &, const std::vector<MAT::Node*> &, tbb::concurrent_hash_map<MAT::Node*, double> &, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, std::vector<bool> &, const bool&);

void updateParsimony(struct min_parsimony &, const std::vector<MAT::Mutation> &, const int &);

void analyzeReads(
    const MAT::Tree &T, 
    const std::string &ref_seq, 
    std::unordered_map<size_t, struct read_info*> &read_map,
    const std::vector<std::string> &vcf_samples,
    const std::string &condensed_nodes_csv,
    const std::string &barcode_file,
    const std::string &read_mutation_depth_vcf
);
