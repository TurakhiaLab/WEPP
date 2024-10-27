#pragma once

#include <vector>
#include <tbb/parallel_sort.h>

#include "panman/panmanUtils.hpp"

#include "config.hpp"
#include "read.hpp"
#include "util.hpp"

struct haplotype {
    int depth;
    /* children and mutations  */
    haplotype* parent;
    std::vector<mutation> muts;
    std::vector<std::pair<int, int>> n_muts;
    std::vector<mutation> stack_muts;
    std::vector<std::pair<int, int>> stack_n_muts;

    // original tree
    std::string id;
    panmanUtils::Node* condensed_source;

    std::vector<haplotype*> children;

    const std::string* reference;

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

        // Check muts
        auto it = std::lower_bound(muts.begin(), muts.end(), search);
        if (it != muts.end() && it->pos <= end)
            return true;
        // Check n_muts
        else
        {
            auto it = std::lower_bound(n_muts.begin(), n_muts.end(), std::make_pair(start, INT_MIN),
                [](const std::pair<int, int>& range, const std::pair<int, int>& value) {
                    return range.second < value.first; 
            });
            if (it != n_muts.end() && it->first <= end)
                return true;
            else 
                return false;
        }
    }

    std::vector<mutation> get_current_mutations(int start = 0, int end = INT_MAX) const {
        std::vector<mutation> muts_in_range;
        for (const mutation &mut: muts) {
            if (mut.pos >= start && mut.pos <= end)
                muts_in_range.emplace_back(mut);
            else if (mut.pos > end)
                break;
        }
        
        // Store Ns in muts_in_range
        for (const auto& coords: n_muts) {
            int start_pos, end_pos;
            if ((end < coords.first) || (start > coords.second))
                continue;
            else if (start <= coords.first) {
                start_pos = coords.first;
                if (end >= coords.second) {
                    end_pos = coords.second;
                }
                else {
                    end_pos = end;
                }
            }
            else {
                start_pos = start;
                if (end >= coords.second) {
                    end_pos = coords.second;
                }
                else {
                    end_pos = end;
                }
            } 
            
            // Add N muts in muts_in_range
            for (int i = start_pos; i <= end_pos; i++)
            {
                mutation m;
                m.pos = i;
                m.ref = nuc_from_char(reference->at(i - 1));
                m.mut = NUC_N;
                muts_in_range.push_back(m);
            }
        }
        tbb::parallel_sort(muts_in_range.begin(), muts_in_range.end());

        return muts_in_range;
    }
    
    std::vector<mutation> get_mutations(int start = 0, int end = INT_MAX) const {
        std::vector<mutation> muts_in_range;
        for (const mutation &mut: stack_muts) {
            if (mut.pos >= start && mut.pos <= end)
                muts_in_range.emplace_back(mut);
            else if (mut.pos > end)
                break;
        }
        
        // Store Ns in muts_in_range
        for (const auto& coords: stack_n_muts) {
            int start_pos, end_pos;
            if ((end < coords.first) || (start > coords.second))
                continue;
            else if (start <= coords.first) {
                start_pos = coords.first;
                if (end >= coords.second) {
                    end_pos = coords.second;
                }
                else {
                    end_pos = end;
                }
            }
            else {
                start_pos = start;
                if (end >= coords.second) {
                    end_pos = coords.second;
                }
                else {
                    end_pos = end;
                }
            } 
            
            // Add N muts in muts_in_range
            for (int i = start_pos; i <= end_pos; i++)
            {
                mutation m;
                m.pos = i;
                m.ref = nuc_from_char(reference->at(i - 1));
                m.mut = NUC_N;
                muts_in_range.push_back(m);
            }
        }
        tbb::parallel_sort(muts_in_range.begin(), muts_in_range.end());
        
        return muts_in_range;
    }

    int mutation_distance(const raw_read& other) const {
        return ::mutation_distance(other.mutations, this->get_mutations(other.start, other.end));
    }

    int mutation_distance(haplotype * other) const {
        return ::mutation_distance(other->get_mutations(), this->get_mutations());
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
