#pragma once

#include <vector>

#include "arena.hpp"
#include "haplotype.hpp"

class initial_filter {
public:
    virtual std::vector<haplotype*> filter(arena& arena) = 0;
    virtual ~initial_filter() { };
};

class wepp_filter: public initial_filter {
    // optimization settings (no effect on final outpiut)
    size_t max_cached_epp_size = 2048;
    bool high_memory_cartesian_map = true;
    // a mutex corresponds to how many haplotypes?
    int mutex_bin_size = 4096;
    // a thread will process about xx total chunks in its lifespan
    int grain_size_factor = 4;

    // affects final output
    double read_dist_factor_threshold = (double) 0.5 / 100;
    int max_peak_peak_mutation = 4;
    int max_peak_nonpeak_mutation = 2;
    int top_n = 25;
    int max_peaks = 300;
    // for a given peak
    int max_neighbors = 50;

    // given a read index
    // what is its maximum parismony score?
    // i.e max_parismony[i] = best parsimony of ith read
    std::vector<int> max_parismony;
    // given a read index
    // how many epp positions does it map to
    // parsimony_multiplicity[i] = epps of ith read
    std::vector<int> parsimony_multiplicity;
    
    /* if epp_positions_cache[i].size() != multiplicity[i], that means not cached */
    /* since there's too many. */
    std::vector<std::vector<haplotype*>> epp_positions_cache;

    std::set<int> remaining_reads;

    void reset(size_t num_reads) {
        max_parismony.clear();
        max_parismony.resize(num_reads);
        parsimony_multiplicity.clear();
        parsimony_multiplicity.resize(num_reads);
        epp_positions_cache.clear();
        epp_positions_cache.resize(num_reads);

        remaining_reads.clear();
        for (int i = 0; i < (int) num_reads; ++i) remaining_reads.insert(i);
    }

    void cartesian_map(arena& arena, std::vector<haplotype*>& haps, const std::vector<raw_read>& read);
    std::vector<int> find_correspondents(arena& arena, haplotype* hap);
    void remove_read(arena &arena, int read_index, std::vector<tbb::queuing_mutex> &mutex);
    void singular_step(arena &arena, haplotype *hap);

    bool valid_two_tops(haplotype *a, haplotype *b) {
         return a->mutation_distance(b) > max_peak_peak_mutation;
    }
    void clear_neighbors(arena& arena, const std::vector<haplotype*>& consideration, std::set<haplotype*>& peaks, std::set<haplotype*>& nbrs);

    bool step(arena& arena, std::vector<haplotype*>& current, std::set<haplotype*> &peaks, std::set<haplotype*> &nbrs);

    double node_score(int parsimony, int epps, int degree)
    {
        return (double) degree / ((1 + parsimony) * epps);
    }

public:
    std::vector<haplotype*> filter(arena& arena);
};

class lineage_root_filter: public initial_filter {
public:
    std::vector<haplotype*> filter(arena& arena);
};

class lineage_root_except_filter: public initial_filter {
public:
    std::string removed_id;
    int clearance = 4;
    std::vector<haplotype*> filter(arena& arena);
};

class random_nodes_filter: public initial_filter {
public:
    int n = 2000;
    std::vector<haplotype*> filter(arena& arena);
};

class uniform_nodes_filter: public initial_filter {
    void dfs(haplotype* curr, haplotype* last_set, std::vector<haplotype*>& dump);
public:
    int dist = 30;
    std::vector<haplotype*> filter(arena& arena);
};