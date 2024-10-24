#pragma once

#include <vector>

#include "panman/panmanUtils.hpp"

#include "config.hpp"
#include "read.hpp"
#include "util.hpp"

struct haplotype {
    int depth;
    /* children and mutations  */
    haplotype* parent;
    std::vector<mutation> muts;
    std::vector<mutation> stack_muts;

    // original tree
    std::string id;
    panmanUtils::Node* condensed_source;

    std::vector<haplotype*> children;

    // orig score
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
        mutation search;
        search.pos = start;

        auto it = std::lower_bound(muts.begin(), muts.end(), search);
        return it != muts.end() && it->pos <= end;
    }

    std::vector<mutation> get_mutations(int start = 0, int end = INT_MAX) const {
        std::vector<mutation> muts_in_range;
        for (const mutation &mut : stack_muts) {
            if (mut.pos >= start && mut.pos <= end)
                muts_in_range.emplace_back(mut);
        }
        return muts_in_range;
    }

    int mutation_distance(const raw_read& other) const {
        return ::mutation_distance(other.mutations, this->get_mutations(other.start, other.end));
    }

    int mutation_distance(haplotype * other) const {
        return ::mutation_distance(other->stack_muts, this->stack_muts);
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

    int mutation_distance(std::vector<mutation> const& comp, int min_pos, int max_pos) {
        return ::mutation_distance(comp, root->get_mutations(min_pos, max_pos));
    }
};
