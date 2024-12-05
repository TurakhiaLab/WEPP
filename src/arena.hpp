#pragma once

#include <vector>
#include <unordered_set>
#include <cmath>
#include <random>

#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/concurrent_hash_map.h>

#include "panman/panmanUtils.hpp"

#include "config.hpp"
#include "haplotype.hpp"
#include "dataset.hpp"
#include "read.hpp"
#include "sam2pb.hpp"
#include "panman_bridge.hpp"
#include "util.hpp"


/* used to sort nodes by their score */
struct score_comparator {
    bool operator() (haplotype* left, haplotype* right) const {
        double const effective_left = left->full_score(); 
        double const effective_right = right->full_score(); 

        if (abs(effective_left - effective_right) > SCORE_EPSILON) {
            return effective_left > effective_right;
        }
        else {
            return left->id > right->id;
        }
    }
};

class arena {
    std::vector<haplotype> nodes;
    std::vector<multi_haplotype> ranged_nodes;
    std::vector<raw_read> raw_reads;
    std::vector<int> masked_sites;

    /* stats */
    size_t read_distribution_bin_size;
    size_t num_reads; // (before merge)
    std::map<std::pair<int, int>, int> ranged_root_map;
    std::array<int, NUM_RANGE_BINS> true_read_counts;
    std::array<double, NUM_RANGE_BINS> true_read_distribution;

    const dataset& ds;

    // note: this is the uncondensed tree
    const panmanUtils::Tree &mat;
    coord_converter coord;
    std::unordered_map<haplotype *, std::vector<panmanUtils::Node *>> condensed_node_mappings;

    haplotype* from_pan(haplotype* parent, panmanUtils::Node* node, const std::unordered_set<int> &site_read_map, std::vector<panmanUtils::Node *> &parent_mapping, const std::vector<mutation>& condensed_n_muts);
    int pan_tree_size(panmanUtils::Node *node);
    int build_range_tree(int parent, haplotype* curr, int start, int end);

    void get_single_mutations(std::vector<mutation>& mutations, const panmanUtils::Node* n);


public:
    arena(const dataset& ds);

    void reset_haplotype_state() {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i != range.end(); i++) {
                    nodes[i].reset_state();
                }
            }
        );
    }

    void recover_haplotype_state() {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i != range.end(); i++) {
                    nodes[i].recover_state();
                }
            }
        );
    }

    void build_range_trees();

    const ::dataset& owned_dataset() {
        return this->ds;
    }

    const std::string& reference() {
        return ds.reference();
    }

    size_t genome_size() {
        return this->reference().size();
    }

    size_t read_dist_bin_size() {
        return read_distribution_bin_size;
    }

    std::vector<multi_haplotype>& ranged_haplotypes() {
        return this->ranged_nodes;
    }

    std::vector<haplotype>& haplotypes() {
        return this->nodes;   
    }

    std::vector<haplotype*> haplotype_pointers() {
        std::vector<haplotype*> ret;
        std::transform(this->nodes.begin(), this->nodes.end(), std::back_inserter(ret),
            [](haplotype& hap) {
                return &hap;
            }
        );
        return ret;
    }

    const std::vector<raw_read>& reads() const {
        return this->raw_reads;
    }

    const size_t num_reads_pre_merge() const {
        return this->num_reads;
    }

    const std::unordered_set<int>& site_read_map() const {
        static std::unordered_set<int> ret;
        if (ret.empty()) {
            for (size_t i = 0; i < this->raw_reads.size(); i++)
            {
                const auto &rp = raw_reads[i];
                std::vector<int> ambiguous_sites;
                for (auto& mut: rp.mutations) {
                    if (mut.mut == NUC_N)
                        ambiguous_sites.emplace_back(mut.pos);
                }

                for (int j = rp.start; j <= rp.end; j++)
                    if (std::find(ambiguous_sites.begin(), ambiguous_sites.end(), j) == ambiguous_sites.end() && std::find(this->masked_sites.begin(), this->masked_sites.end(), j) == this->masked_sites.end())
                        ret.insert(j);
            }
        }
        return ret;
    }

    const std::array<int, NUM_RANGE_BINS>& read_counts() const {
        return this->true_read_counts;
    }

    const std::array<double, NUM_RANGE_BINS>& read_distribution() const {
        return this->true_read_distribution;
    }

    const std::vector<panmanUtils::Node*>& source_nodes(haplotype* haplotype) {
        return this->condensed_node_mappings[haplotype];
    }

    multi_haplotype *find_range_tree_for(const raw_read &read);

    // note: this is by TREE edge distance as opposed to mutation distance
    std::set<haplotype*, score_comparator> closest_neighbors(haplotype * target, int max_radius, int num_limit) const; 

    // greatest score first
    std::set<haplotype*, score_comparator> highest_scoring_neighbors(haplotype* target, bool include_mapped, int max_radius, int num_limit) const;

    haplotype* haplotype_with_id(const std::string& id) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (panmanUtils::Node* source: this->condensed_node_mappings[&nodes[i]])
            if (source->identifier == id) {
                return &nodes[i];
            }
        }
        
        std::cerr << " Could not find " << id;
        return nullptr;
    }

    void get_residual_cooccuring_mutations(int window_size);

    // mutations from root to here
    std::vector<mutation> get_mutations(const panmanUtils::Node* n, bool replace_GAP=false, bool replace_N=false);

    // precondition: true haplotypes of current dataset are known
    void print_mutation_distance(const std::vector<haplotype*>& selected);
    
    void print_flipped_mutation_distance(const std::vector<std::pair<haplotype *, double>> &selected);
    
    void print_full_report(const std::vector<std::pair<haplotype*, double>> & abundance);

    void dump_haplotype_proportion(const std::vector<std::pair<haplotype*, double>> & abundance); 
    
    void resolve_unaccounted_mutations(const std::vector<std::pair<haplotype*, double>> & abundance);  

    void dump_read2haplotype_mapping(const std::vector<std::pair<haplotype*, double>> & abundance);
};