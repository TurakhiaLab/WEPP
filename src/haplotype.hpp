#pragma once

#include <vector>

#include "panman/panmanUtils.hpp"

#include "config.hpp"
#include "read.hpp"

struct haplotype {
    /* children and mutations  */
    haplotype* parent;
    /* muts from root to here */
    std::vector<mutation> muts;
    std::vector<mutation> stack_muts;

    // original tree
    std::string id;
    panmanUtils::Node* condensed_source;

    std::vector<haplotype*> children;

    // raw score
    double score;
    // 'f' score, where it rewards
    // nodes who's corresponding reads come from a 
    // wide variety of genome places
    double dist_divergence = 1;
    // whether or not it has been selected 
    // as a peak or neighbor
    bool mapped;
    // note: does not handle overlapping ranges
    // instead, the i'th bucket corresponds to the interval
    // i * max_genome_size / buckets to (i + 1) * max_genome_size / buckets
    std::array<int, NUM_RANGE_BINS> mapped_read_counts{};

    void reset_state() {
        this->score = 0;
        this->dist_divergence = 1;
        this->mapped = false;
        std::fill(this->mapped_read_counts.begin(), this->mapped_read_counts.begin() + NUM_RANGE_BINS, 0);
    }

    bool has_mutations_in_range(int start, int end) {
        mutation search;
        search.pos = start;

        auto it = std::lower_bound(muts.begin(), muts.end(), search);
        return it != muts.end() && it->pos <= end;
    }

    // given a range and reference mutation list
    // tells the (genome) positions where the this mutation list differs from reference
    // min_pos and max_pos are genome positions
    std::vector<int> mutations(std::vector<mutation> const& comp, int min_pos, int max_pos) const {
        std::vector<int> muts;

        mutation search;
        search.pos = min_pos;
        int i = std::lower_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        search.pos = max_pos;
        int last_i = std::upper_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        int j = 0;

        while (i < last_i || j < (int) comp.size()) {
            if (i == last_i) {                
                if (comp[j].mut != NUC_N) {
                    muts.push_back(comp[j].pos);
                }
                ++j;
            }
            else if (stack_muts[i].pos < min_pos) {
                ++i;
            }
            else if (stack_muts[i].pos > max_pos) {
                assert(0);
                return muts;
            }
            else if (j == (int) comp.size()) {
                if (stack_muts[i].mut != NUC_N) {
                    muts.push_back(stack_muts[i].pos);
                }
                ++i;
            }
            else if (stack_muts[i].pos < comp[j].pos) {
                if (stack_muts[i].mut != NUC_N) {
                    muts.push_back(stack_muts[i].pos);
                }
                ++i;
            }
            else if (stack_muts[i].pos > comp[j].pos) {
                if (comp[j].mut != NUC_N) {
                    muts.push_back(comp[j].pos);
                }
                ++j;
            }
            else if (stack_muts[i].pos == comp[j].pos && stack_muts[i].mut != comp[j].mut && stack_muts[i].mut != NUC_N && comp[j].mut != NUC_N) {
                muts.push_back(comp[j].pos);
                ++i; ++j;
            }
            else {
                ++i; ++j;
            }
        }

        return muts;
    }

    // called so much it's worth optimizing out the vector
    int mutation_distance(std::vector<mutation> const& comp, int min_pos, int max_pos) const {
        int muts = 0;

        mutation search;
        search.pos = min_pos;
        int i = std::lower_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        search.pos = max_pos;
        int last_i = std::upper_bound(stack_muts.begin(), stack_muts.end(), search) - stack_muts.begin();
        int j = 0;

        while (i < last_i || j < (int) comp.size()) {
            if (i == last_i) {                
                if (comp[j].mut != NUC_N) {
                    ++muts;
                }
                ++j;
            }
            else if (stack_muts[i].pos < min_pos) {
                ++i;
            }
            else if (stack_muts[i].pos > max_pos) {
                assert(0);
                return muts;
            }
            else if (j == (int) comp.size()) {
                if (stack_muts[i].mut != NUC_N) {
                    ++muts;
                }
                ++i;
            }
            else if (stack_muts[i].pos < comp[j].pos) {
                if (stack_muts[i].mut != NUC_N) {
                    ++muts;
                }
                ++i;
            }
            else if (stack_muts[i].pos > comp[j].pos) {
                if (comp[j].mut != NUC_N) {
                    ++muts;
                }
                ++j;
            }
            else if (stack_muts[i].pos == comp[j].pos && stack_muts[i].mut != comp[j].mut && stack_muts[i].mut != NUC_N && comp[j].mut != NUC_N) {
                ++muts;
                ++i; ++j;
            }
            else {
                ++i; ++j;
            }
        }

        return muts;
    }

    int mutation_distance(const raw_read& other) const {
        return mutation_distance(other.mutations, other.start, other.end);
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

    std::vector<int> mutations(std::vector<mutation> const& comp, int min_pos, int max_pos) {
        return root->mutations(comp, min_pos, max_pos);
    }

    int mutation_distance(std::vector<mutation> const& comp, int min_pos, int max_pos) {
        return root->mutation_distance(comp, min_pos, max_pos);
    }
};
