#pragma once

#include <vector>

#include "src/usher_graph.hpp"

#include "config.hpp"
#include "read.hpp"

struct haplotype {
    // based on uncompressed tree
    int leaf_count;
    bool is_leaf;

    /* children and mutations  */
    haplotype* parent;
    /* muts from root to here */
    std::vector<MAT::Mutation> muts;
    std::vector<MAT::Mutation> stack_muts;

    // original tree
    std::string id;
    MAT::Node* condensed_source;

    std::vector<haplotype*> children;

    //orig score
    double orig_score;
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
        this->orig_score = 0;
        this->score = 0;
        this->dist_divergence = 1;
        this->mapped = false;
        std::fill(this->mapped_read_counts.begin(), this->mapped_read_counts.begin() + NUM_RANGE_BINS, 0);
    }

    void recover_state() {
        this->score = this->orig_score;
        this->mapped = false;
        std::fill(this->mapped_read_counts.begin(), this->mapped_read_counts.begin() + NUM_RANGE_BINS, 0);
    }

    bool has_mutations_in_range(int start, int end) {
        MAT::Mutation search;
        search.position = start;
        search.mut_nuc = 0;

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
        search.mut_nuc = 0;
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
            else if (stack_muts[i].position == comp[j].position && stack_muts[i].mut_nuc != comp[j].mut_nuc && comp[j].mut_nuc != unknown_nuc) {
                muts.push_back(comp[j].position);
                ++i; ++j;
            }
            else {
                ++i; ++j;
            }
        }

        return muts;
    }

    // called so much it's worth optimizing out the vector
    int mutation_distance(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) const {
        int muts = 0;

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
                    ++muts;
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
                ++muts;
                ++i;
            }
            else if (stack_muts[i].position < comp[j].position) {
                ++muts;
                ++i;
            }
            else if (stack_muts[i].position > comp[j].position) {
                if (comp[j].mut_nuc != unknown_nuc) {
                    ++muts;
                }
                ++j;
            }
            else if (stack_muts[i].position == comp[j].position && stack_muts[i].mut_nuc != comp[j].mut_nuc && comp[j].mut_nuc != unknown_nuc) {
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

    std::vector<MAT::Mutation> mutation_distance_vector(const raw_read& other) const {
        std::vector<MAT::Mutation> different_mutations, hap_range_mutations, read_mutations = other.mutations;
        for (auto mut: this->stack_muts) {
            if ((mut.position >= other.start) && (mut.position <= other.end))
                hap_range_mutations.emplace_back(mut);
        }
        
        tbb::parallel_sort(read_mutations.begin(), read_mutations.end(), [](const MAT::Mutation &a, const MAT::Mutation &b) {
            if (a.position != b.position)
                return a.position < b.position;
            else
                return a.mut_nuc < b.mut_nuc;
        });

        tbb::parallel_sort(hap_range_mutations.begin(), hap_range_mutations.end(), [](const MAT::Mutation &a, const MAT::Mutation &b) {
            if (a.position != b.position)
                return a.position < b.position;
            else
                return a.mut_nuc < b.mut_nuc;
        });

        auto read_iterator = read_mutations.begin();
        auto hap_iterator = hap_range_mutations.begin();
        while (read_iterator != read_mutations.end() && hap_iterator != hap_range_mutations.end()) {
            if (read_iterator->position == hap_iterator->position) {
                if ((read_iterator->mut_nuc != 0b1111) && (read_iterator->mut_nuc != hap_iterator->mut_nuc)) {
                    different_mutations.emplace_back(*read_iterator);
                }
                read_iterator++;
                hap_iterator++;
            } else if (read_iterator->position < hap_iterator->position) {
                if (read_iterator->mut_nuc != 0b1111)
                    different_mutations.emplace_back(*read_iterator);
                read_iterator++;
            } else if (hap_iterator->position < read_iterator->position) {
                different_mutations.emplace_back(*hap_iterator);
                hap_iterator++;
            }
        }

        while (read_iterator != read_mutations.end())
        {
            if (read_iterator->mut_nuc != 0b1111)
                different_mutations.emplace_back(*read_iterator);
            read_iterator++;
        }

        while (hap_iterator != hap_range_mutations.end())
        {   
            different_mutations.emplace_back(*hap_iterator);
            hap_iterator++;
        }

        return different_mutations;
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

    int mutation_distance(std::vector<MAT::Mutation> const& comp, int min_pos, int max_pos) {
        return root->mutation_distance(comp, min_pos, max_pos);
    }
};