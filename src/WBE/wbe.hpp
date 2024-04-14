#include "src/usher_graph.hpp"
#include <regex>
#include <array>
#include <fstream>
#include <iostream>
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

struct SAM_read {
    std::string raw_name;
    int start_idx;
    int degree; // how many repeats;
    std::string aligned_string;

    std::string degree_name() const {
        return raw_name + 
            "_READ_" + std::to_string(start_idx + 1) + 
            "_" + std::to_string(start_idx + 1 + aligned_string.size() - 1) + 
            "_" + std::to_string(degree);
    }

    std::string degreeless_name() const {
        return raw_name + 
            "_READ_" + std::to_string(start_idx + 1) + 
            "_" + std::to_string(start_idx + 1 + aligned_string.size() - 1); 
    }

    bool operator <(const SAM_read& rhs) {
        if (start_idx != rhs.start_idx) {
            return start_idx < rhs.start_idx;
        }

        if (aligned_string.size() != rhs.aligned_string.size()) {
            return aligned_string.size() < rhs.aligned_string.size();
        }

        return aligned_string < rhs.aligned_string;
    }

    bool operator ==(const SAM_read& rhs) {
        return start_idx == rhs.start_idx && aligned_string == rhs.aligned_string;
    }
};

const std::string GENOME_STRING{"ACGTN"};
typedef std::vector<std::array<int, 5>> sub_table;

struct SAM {
   private:
       const std::string reference_seq;
   
       /* aligned/padded reads (possibly merged) */
       std::vector<SAM_read> aligned_reads;
   
       /*
         these are the only the things that we need,
         since we can derive all other information from here
        */
   
       /* for any index, highlight all possible SNP */
       /* does NOT take into account read correction */
       sub_table frequency_table;
   
       /* frequency + collapsed indels to reference */
       /* may take into account read correction */
       sub_table collapsed_frequency_table;
   
       /* for any starting index, length possible indels (length in reference, replacement string)  */
       std::vector<std::map<std::pair<size_t, std::string>, int>> indel_frequency_table;
   
       /* merging unmap information, column -> preimage of column */
       std::map<std::string, std::vector<std::string>> reverse_merge;
   
       void read_correction();
       void merge_duplicates();
   
   public:
       SAM(const std::string& ref);
   
       void add_read(const std::string& line);
       void build();
   
       void dump_proto(std::string const& filename);
       void dump_vcf(std::ostream& out);
       void dump_reverse_merge(std::ostream& out);
       void dump_freyja(std::ostream& dout, std::ostream& vout);
};

po::variables_map parseWBEcommand(po::parsed_options);

void load_reads_from_proto(std::string const &, std::unordered_map<size_t, struct read_info *> &, std::unordered_map<std::string, std::vector<std::string>> &);

void sam2VCF(po::parsed_options);

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

int mutationDistance(const MAT::Tree &, const MAT::Tree &, const MAT::Node*, const MAT::Node*);

int mutationDistance(std::vector<MAT::Mutation>, std::vector<MAT::Mutation>);

std::vector<MAT::Mutation> getMutations(const MAT::Tree&, const std::string);

std::string getLineage(const MAT::Tree &, MAT::Node*);

void updateProhibitedNodes(const MAT::Tree &, const std::vector<MAT::Node*> &, std::vector<MAT::Node*> &, const int&);

void getProhibitedNodes(const MAT::Tree &, const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Node* , std::vector<MAT::Node*> &, const int&);

std::vector<MAT::Node*> updateNeighborNodes(const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::vector<MAT::Node*> &, const std::vector<MAT::Node*> &, const tbb::concurrent_hash_map<MAT::Node*, double> &, std::vector<MAT::Node*> &, const int&, const int&);

void addNeighborNodes(const MAT::Tree &, std::vector<MAT::Node*> &, const int&, const int&);

void addNeighborLineages(const MAT::Tree &, const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::vector<MAT::Node*> &, std::vector<std::string> &, const int &);

int childNodesAddition(const MAT::Tree &, my_mutex_t &, const std::vector<MAT::Node*> &, std::vector<MAT::Node*> &, int &, std::vector<MAT::Node*> &, MAT::Node *, MAT::Node *, const int &, const int &, const int &, const int &);

void generateFilteringData(const MAT::Tree &, const MAT::Tree &, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::string &, const std::vector<MAT::Node*> &, const std::unordered_map<size_t, struct read_info*> &, const std::string &, const std::string &, const std::string &);

void sortNodeScore(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const tbb::concurrent_hash_map<MAT::Node*, double> &, tbb::concurrent_vector<std::pair<MAT::Node*, double>> &);

bool compareNodeScore(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const std::pair<MAT::Node*, double>&, const std::pair<MAT::Node*, double>&);

bool compareReadPos (const std::pair<size_t, std::string>&, const std::pair<size_t, std::string>&);

bool compareMutations(const MAT::Mutation &, const MAT::Mutation &);

size_t getNumLeaves(const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Node*);

bool compareIdx(const std::pair<int, size_t> &, const std::pair<int, size_t> &);

void createRangeTree(MAT::Node*, const std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, const int &, const int &, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Tree &);

void createLineageTree(MAT::Node*, const std::vector<std::string> &, MAT::Tree &);

void createCondensedTree(MAT::Node*, const std::unordered_map<size_t, struct read_info*> &, std::unordered_map<MAT::Node*, std::vector<MAT::Node*>> &, MAT::Tree &);

void computeDistance(const MAT::Tree &, const std::unordered_map<size_t, struct read_info*> &, const std::vector<std::string> &, const std::unordered_map<std::string, double> &, const std::unordered_set<int> &);

void placeReads(const MAT::Tree &, const std::string &, const std::unordered_map<size_t, struct read_info*> &, const std::unordered_map<size_t, struct read_info*> &, const std::unordered_map<std::string, double> &, const std::vector<std::string> &, const std::unordered_set<int> &);
