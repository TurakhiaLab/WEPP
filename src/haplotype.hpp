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

    std::vector<mutation> stack_muts(int start = 0, int end = INT_MAX) const {
        std::vector<mutation> ret;

        if (parent)
        {
            for (const mutation &mut : parent->stack_muts(start, end))
            {
                /* only need to compare to our muts since parent's are unique */
                bool valid = true;
                for (const mutation &comp : muts)
                {
                    if (comp.pos == mut.pos)
                    {
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    ret.push_back(mut);
                }
            }
        }
        /* maybe rewinded mutation back to original */
        for (const mutation &mut : muts)
        {
            if (mut.ref != mut.mut && mut.pos >= start && mut.pos <= end)
            {
                ret.push_back(mut);
            }
        }
        std::sort(ret.begin(), ret.end());
        return ret;
    }
    
    int mutation_distance(const raw_read& other) const {
        return ::mutation_distance(other.mutations, this->stack_muts(other.start, other.end));
    }

    int mutation_distance(haplotype * other) const {
        std::vector<mutation> us, theirs;
        haplotype const* a = this, *b = other;

        auto add_to = [](haplotype const* src, std::vector<mutation>&dump) {
            for (const auto& mut : src->muts) {
                bool valid = true;
                for (const auto& existing : dump) {
                    if (existing.pos == mut.pos) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    dump.push_back(mut);
                }
            }
        };

        while (a->depth < b->depth) {
            add_to(a, us);
            a = a->parent;
        }

        while (b->depth < a->depth) {
            add_to(b, theirs);
            b = b->parent;
        }

        while (a != b) {
            add_to(a, us);
            add_to(b, theirs);
            a = a->parent;
            b = b->parent;
        }

        auto del_backs = [](std::vector<mutation>& muts) {
            auto it = muts.begin();
            while (it < muts.end()) {
                if (it->ref == it->mut) {
                    it = muts.erase(it);
                }
            }
        };

        del_backs(us);
        del_backs(theirs);

        std::sort(us.begin(), us.end());
        std::sort(theirs.begin(), theirs.end());

        return ::mutation_distance(us, theirs);
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
        return ::mutation_distance(comp, root->stack_muts(min_pos, max_pos));
    }
};
