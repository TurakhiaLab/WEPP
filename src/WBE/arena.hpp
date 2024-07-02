#pragma once

#include <vector>
#include <cmath>

#include "src/usher_graph.hpp"

#include "config.hpp"
#include "haplotype.hpp"
#include "dataset.hpp"
#include "pipeline.hpp"
#include "read.hpp"

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

    const dataset& ds;
    const config& config;

public:
    arena(const dataset& ds, const config& config) : ds{ds}, config{config} {
        MAT::Tree tree = ds.mat();

    }

    void closest_neighbors(int max_radius, int num_limit) {

    }

    void highest_scoring_neighbors(bool include_mapped, int max_radius, int num_limit) {

    }

    // precondition: true haplotypes of current dataset are known
    void print_mutation_distance(const std::vector<haplotype*>& selected) {
        
    }
};