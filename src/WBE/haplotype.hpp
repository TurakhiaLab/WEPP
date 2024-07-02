#pragma once

#include <vector>

#include "src/usher_graph.hpp"

#include "read.hpp"

struct haplotype {
    // based on uncompressed tree
    int leaf_count;
    bool is_leaf;

    /* children and mutations  */
    haplotype* parent;
    /* muts from root to here */
    std::vector<MAT::Mutation> stack_muts;

    // original tree
    std::string id;
    MAT::Node* condensed_source;

    std::vector<haplotype*> children;

    // raw score
    double score;
    // 'f' score, where it rewards
    // nodes who's corresponding reads come from a 
    // wide variety of genome places
    double dist_divergence;
    // whether or not it has been selected 
    // as a peak or neighbor
    bool mapped;

    bool has_mutations_in_range(int start, int end) {
        MAT::Mutation search;
        search.position = start;
        auto it = std::lower_bound(condensed_source->mutations.begin(), condensed_source->mutations.end(), search);
        return it != condensed_source->mutations.end() && it->position <= end;
    }

    // given a range and reference mutation list
    // tells the (genome) positions where the this mutation list differs from reference
    // min_pos and max_pos are genome positions
    std::vector<int> mutations(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) const {
        std::vector<int> muts;

        const int unknown_nuc = 0b1111;

        MAT::Mutation search;
        search.position = min_pos;
        int i = std::lower_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        search.position = max_pos;
        int last_i = std::upper_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        int j = 0;

        while (i < last_i || j < (int) comp.size()) {
            if (i == last_i) {                
                if (comp[j].mut_nuc != unknown_nuc) {
                    muts.push_back(comp[j].position);
                }
                ++j;
            }
            else if (stack_muts[i].position < min_pos) {
                ++i;
            }
            else if (stack_muts[i].position > max_pos) {
                assert(0);
                return muts;
            }
            else if (j == (int) comp.size()) {
                muts.push_back(stack_muts[i].position);
                ++i;
            }
            else if (stack_muts[i].position < comp[j].position) {
                muts.push_back(stack_muts[i].position);
                ++i;
            }
            else if (stack_muts[i].position > comp[j].position) {
                if (comp[j].mut_nuc != unknown_nuc) {
                    muts.push_back(comp[j].position);
                }
                ++j;
            }
            else if (stack_muts[i].position == comp[j].position && stack_muts[i].mut_nuc != comp[j].mut_nuc && stack_muts[i].mut_nuc != unknown_nuc) {
                muts.push_back(comp[j].position);
                ++i; ++j;
            }
            else {
                ++i; ++j;
            }
        }

        return muts;
    }

    int mutation_distance(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) const {
        return mutations(comp, min_pos, max_pos).size();
    }

    int mutation_distance(read* other) const {
        return mutation_distance(other->mutations, other->start, other->end);
    }

    int mutation_distance(haplotype * other) const {
        return mutation_distance(other->stack_muts, 0, INT_MAX);
    }

    double full_score() const {
        return score * sqrt(dist_divergence);
    }
};

// range tree compression of haplotypes
struct multi_haplotype {
    haplotype* root;
    /* the original haplotypes that correspond to this range tree */
    std::vector<haplotype*> sources;

    /* arena indices of ranged node children */
    std::vector<int> children;

    std::vector<int> mutations(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) {
        return root->mutations(comp, min_pos, max_pos);
    }
};
