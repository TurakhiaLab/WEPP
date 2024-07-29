#pragma once

#include <vector>
#include <cmath>

#include "panman/panman.hpp"

#include "config.hpp"
#include "haplotype.hpp"
#include "dataset.hpp"
#include "read.hpp"
#include "util.hpp"

/* used to sort nodes by their score */
struct score_comparator {
    bool operator() (haplotype* left, haplotype* right) const {
        double const effective_left = left->full_score(); 
        double const effective_right = right->full_score(); 

        if (abs(effective_left - effective_right) > SCORE_EPSILON) {
            return effective_left > effective_right;
        }
        else if (left->leaf_count != right->leaf_count) {
            return left->leaf_count > right->leaf_count;
        }
        else {
            return left->id > right->id;
        }
    }
};

/* used to sort nodes by their mutation list */
struct mutation_comparator {
    bool operator() (haplotype* left, haplotype* right) const {
        if (left->stack_muts.size() != right->stack_muts.size()) {
            return left->stack_muts.size() < right->stack_muts.size();
        }  
        for (size_t i = 0; i < left->stack_muts.size(); ++i) {
            if (left->stack_muts[i].position != right->stack_muts[i].position) {
                return left->stack_muts[i].position < right->stack_muts[i].position;
            }
            else if (left->stack_muts[i].mut_nuc != right->stack_muts[i].mut_nuc) {
                return left->stack_muts[i].mut_nuc < right->stack_muts[i].mut_nuc;
            }
        }

        return false;
    }
};

class arena {
    std::vector<haplotype> nodes;
    std::vector<multi_haplotype> ranged_nodes;
    std::vector<raw_read> raw_reads;

    /* stats */
    size_t read_distribution_bin_size;
    size_t num_reads; // (before merge)
    std::map<std::pair<int, int>, int> ranged_root_map;
    std::array<int, NUM_RANGE_BINS> true_read_counts;
    std::array<double, NUM_RANGE_BINS> true_read_distribution;

    panmanUtils::Tree mat;
    std::unordered_map<panmanUtils::Node *, std::vector<panmanUtils::Node *>> condensed_node_mappings;

    const dataset& ds;

    haplotype* from_pan(haplotype* parent, panmanUtils::Node* node);
    int pan_tree_size(panmanUtils::Node *node); 
    int build_range_tree(int parent, haplotype* curr, int start, int end);
public:
    arena(const dataset& ds) : ds{ds} {
        this->raw_reads = ds.reads();
        this->mat = ds.mat();
        
        panmanUtils::Tree condensed = create_condensed_tree(this->mat.root, this->raw_reads, this->condensed_node_mappings);

        // create vanilla nodes
        this->nodes.reserve(pan_tree_size(condensed.root));
        this->from_pan(nullptr, condensed.root);
    }

    void reset_haplotype_state() {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, nodes.size()),
            [&](tbb::blocked_range<size_t> range) {
                for (size_t i = range.begin(); i != range.end(); i++) {
                    nodes[i].reset_state();
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
        return ds.reference().size();
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

    const std::array<int, NUM_RANGE_BINS>& read_counts() const {
        return this->true_read_counts;
    }

    const std::array<double, NUM_RANGE_BINS>& read_distribution() const {
        return this->true_read_distribution;
    }

    const std::vector<panmanUtils::Node*>& source_nodes(haplotype* haplotype) {
        return this->condensed_node_mappings[haplotype->condensed_source];
    }

    multi_haplotype *find_range_tree_for(const raw_read &read);

    // note: this is by TREE edge distance as opposed to mutation distance
    std::set<haplotype*, mutation_comparator> closest_neighbors(haplotype * target, int max_radius, int num_limit) const; 

    // greatest score first
    std::set<haplotype*, score_comparator> highest_scoring_neighbors(haplotype* target, bool include_mapped, int max_radius, int num_limit) const;

    haplotype* haplotype_with_id(const std::string& id) {
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (panmanUtils::Node* source: this->condensed_node_mappings[nodes[i].condensed_source])
            if (source->identifier == id) {
                return &nodes[i];
            }
        }
        
        std::cerr << " Could not find " << id;
        return nullptr;
    }

    void print_cooccuring_mutations(int window_size);

    // precondition: true haplotypes of current dataset are known
    void print_mutation_distance(const std::vector<haplotype*>& selected);
    void print_flipped_mutation_distance(const std::vector<haplotype *> &selected);
    void print_full_report(const std::vector<std::pair<haplotype*, double>> & abundance);

    void dump_read2node_mapping(const std::vector<std::pair<haplotype*, double>> & abundance); 
};