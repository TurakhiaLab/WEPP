#include <algorithm>
#include <chrono>
#include <random>

#include "initial_filter.hpp"

static constexpr uint8_t N = 0b1111;
/* For a fixed read, we are given the positions where the parent node */
/* differs from the read, and the positions where the current node */
/* differs from the read. This function counts how many difference the parent has */
/* that the child does not. This is NOT the same as the absolute difference in vector sizes */
/* since in this function, positions the child has that the parent does not */
/* are of no relevance */
/* (this is to match reference behavior: can only map if one or mutation positions was removed */
/* even if total size increases) */
static int
mutation_reductions(std::vector<int> const &parent, std::vector<int> const &child)
{
    int reductions = 0;
    for (size_t i = 0, j = 0; i < parent.size(); ++i)
    {
        while (j < child.size() && child[j] < parent[i])
        {
            ++j;
        }
        if (j == child.size() || child[j] != parent[i])
        {
            ++reductions;
        }
        else
        {
            ++j;
        }
    }

    return reductions;
}

/* maps a single read to entire tree, caching aliveness */
/* parent_locations is the positions where the parent has different mutations than the read */
static void
single_read_tree(arena& arena, const std::vector<int>& parent_locations, multi_haplotype *curr, const raw_read& read, std::vector<multi_haplotype *> &max_nodes, int &max_val)
{
    // positions where this node differs from the read */
    std::vector<int> my_locations;
    my_locations.reserve(parent_locations.size());

    // go through all of this haplotypes mutations
    // and see if it increases or decreases overall mutations
    std::vector<MAT::Mutation>& curr_muts = curr->root->muts;

    size_t i = 0; // index of parent location
    MAT::Mutation search;
    search.position = read.start;
    auto j = std::lower_bound(curr_muts.begin(), curr_muts.end(), search); // index of curr_muts

    // keeping track of read mutation position may actually be the dominating factor with so many N
    // so we bin search whenever necessary
    while (i < parent_locations.size() || (j != curr_muts.end() && j->position <= read.end)) {
        bool parent_first = i < parent_locations.size() && (j == curr_muts.end() || j->position > read.end || parent_locations[i] < j->position);
        bool us_first = (j != curr_muts.end() && j->position <= read.end) && (i == parent_locations.size() || j->position < parent_locations[i]);

        if (us_first) {
            auto it = std::lower_bound(read.mutations.begin(), read.mutations.end(), *j);
            uint8_t const read_nuc = it == read.mutations.end() || it->position != j->position ? j->ref_nuc : it->mut_nuc;
            if (read_nuc != N && read_nuc != j->mut_nuc) {
                my_locations.push_back(j->position);
            }
            
            ++j;
        }
        else if (parent_first) {
            my_locations.push_back(parent_locations[i]);

            ++i;
        }
        else {
            auto it = std::lower_bound(read.mutations.begin(), read.mutations.end(), *j);
            uint8_t const read_nuc = it == read.mutations.end() || it->position != j->position ? j->ref_nuc : it->mut_nuc;
            if (read_nuc != N && read_nuc != j->mut_nuc) {
                my_locations.push_back(j->position);
            }
            
            ++j;
            ++i;
        }
    }

    /* map self */
    int parsimony = my_locations.size();
    if (parsimony < max_val)
    {
        max_val = parsimony;
        max_nodes = {curr};
    }
    else if (parsimony == max_val)
    {
        max_nodes.emplace_back(curr);
    }

    for (int child : curr->children)
    {
        single_read_tree(arena, my_locations, &arena.ranged_haplotypes()[child], read, max_nodes, max_val);
    }
}

/* maps a single read to entire tree, caching aliveness */
/* and finding the epp positions of a given read */
/* parent_locations is the positions where the parent has different mutations than the read */
/* max_val corresponds to max parismony */
/* max indices correspond to the indices of the epp nodes */
static void 
single_read_tree(arena& arena, const raw_read& read, std::vector<haplotype*> &max_indices, int &max_val)
{
    std::vector<multi_haplotype *> max_nodes;
    multi_haplotype* root = arena.find_range_tree_for(read);
    
    std::vector<int> root_mutations;
    for (const MAT::Mutation& mut: read.mutations) {
        if (mut.mut_nuc != N) {
            root_mutations.push_back(mut.position);
        }
    }
    single_read_tree(arena, root_mutations, root, read, max_nodes, max_val);

    for (multi_haplotype *rnode : max_nodes)
    {
        for (haplotype *src : rnode->sources)
        {
            if (!src->mapped) {
                max_indices.push_back(src);
            }
        }
    }
}

// calculate initial scores as well as divergences
// haps assumed to be arena indices
void 
wepp_filter::cartesian_map(arena& arena, std::vector<haplotype*>& haps, const std::vector<raw_read>& reads)
{
    tbb::queuing_mutex my_mutex;

    Timer timer; timer.Start();

    int bin_size = arena.genome_size() / NUM_RANGE_BINS;
    int num_threads = arena.owned_dataset().num_threads();
    // avoid many flushes of score and reallocations of the score vector
    // this does come at the cost of some threads missing out on work at the end
    int grain_size = std::max((int)(reads.size() / num_threads / GRAIN_SIZE_FACTOR), 1);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size(), grain_size),
                      [&](tbb::blocked_range<size_t> k)
                      {
                          std::vector<std::pair<double, std::array<int, NUM_RANGE_BINS>>> score;
                          if (HIGH_MEMORY_CARTESIAN_MAP) {
                              score.resize(haps.size()); // otherwise, not used and buffered directly
                          }

                          for (size_t r = k.begin(); r != k.end(); ++r)
                          {
                              std::vector<haplotype *> max_indices;
                              int max_val = INT32_MAX; /* really want minimum value */

                              single_read_tree(arena, reads[r], max_indices, max_val);

                              double const delta = node_score(max_val, max_indices.size(), reads[r].degree);
                              {
                                  int bucket = reads[r].start / bin_size;
                                  if (HIGH_MEMORY_CARTESIAN_MAP) {
                                      haplotype *arena_base = &arena.haplotypes()[0];
                                      for (haplotype *hap : max_indices)
                                      {
                                          int index = (hap - arena_base);
                                          score[index].first += delta;
                                          score[index].second[bucket] += reads[r].degree;
                                      }
                                  }
                                  else {
                                      tbb::queuing_mutex::scoped_lock my_lock{my_mutex};
                                      for (haplotype *hap : max_indices)
                                      {
                                          hap->score += delta;
                                          hap->mapped_read_counts[bucket] += reads[r].degree;
                                      }
                                  }
                              }

                              this->max_parismony[r] = max_val;
                              this->parsimony_multiplicity[r] = max_indices.size();
                              if (max_indices.size() <= MAX_CACHED_EPP_SIZE)
                              {
                                  // ensure it is sorted
                                  std::sort(max_indices.begin(), max_indices.end());
                                  this->epp_positions_cache[r] = std::move(max_indices);
                              }
                          }

                          if (HIGH_MEMORY_CARTESIAN_MAP)
                          {
                              // only use a global mutex in this case
                              tbb::queuing_mutex::scoped_lock my_lock{my_mutex};
                              for (size_t i = 0; i < score.size(); ++i)
                              {
                                  haps[i]->score += score[i].first;
                                  for (size_t j = 0; j < NUM_RANGE_BINS; ++j)
                                  {
                                      haps[i]->mapped_read_counts[j] += score[i].second[j];
                                  }
                              }
                          }
                      });

    for (size_t i = 0; i < haps.size(); ++i)
    {
        /* calculate divergence */
        int divergence = 0, bins_active = 0;
        for (size_t j = 0; j < NUM_RANGE_BINS; ++j)
        {
            if (arena.read_counts()[j])
            {
                bins_active += 1;
            }
            double const proportion = (double) haps[i]->mapped_read_counts[j] / arena.read_counts()[j];
            if (proportion > READ_DIST_FACTOR_THRESHOLD)
            {
                divergence += 1;
            }
        }

        haps[i]->dist_divergence = (double) divergence / bins_active;
        haps[i]->orig_score = haps[i]->score;
    }

    /* initial sort into scores */
    tbb::parallel_sort(haps.begin(), haps.end(), score_comparator());

    std::cout << "--- cartesian mapping took " << (timer.Stop() / 1000) << " seconds " << std::endl;
}

std::vector<int>
wepp_filter::find_correspondents(arena& arena, haplotype* hap)
{
    std::vector<int> correspondents;

    tbb::queuing_mutex my_mutex;

    std::vector<int> remaining(remaining_reads.begin(), remaining_reads.end());
    tbb::parallel_for(tbb::blocked_range<size_t>(0, remaining.size()),
                      [&](tbb::blocked_range<size_t> r)
                      {
                          for (size_t i = r.begin(); i != r.end(); ++i)
                          {
                              size_t const read = remaining[i];
                              /* if cached correspond, use that */
                              auto it = std::lower_bound(epp_positions_cache[read].begin(), epp_positions_cache[read].end(), hap);
                              if (it != epp_positions_cache[read].end() && *it == hap)
                              {
                                  tbb::queuing_mutex::scoped_lock lock{my_mutex};
                                  correspondents.push_back(read);
                              }
                              else if (epp_positions_cache[read].size() == (size_t)parsimony_multiplicity[read])
                              {
                                  /* cached and not found */
                                  continue;
                              }
                              else
                              {
                                  /* recompute to see if correspondent */
                                  const raw_read &r = arena.reads()[read];
                                  int parsimony = hap->mutation_distance(r.mutations, r.start, r.end);
                                  if (parsimony == max_parismony[read])
                                  {
                                      tbb::queuing_mutex::scoped_lock lock{my_mutex};
                                      correspondents.push_back(read);
                                  }
                              }
                          }
                      });
    tbb::parallel_sort(correspondents.begin(), correspondents.end());
    return correspondents;
}

void
wepp_filter::remove_read(arena& arena, int read_index, std::vector<tbb::queuing_mutex>& mutex)
{
    using mutex_t = tbb::queuing_mutex;

    std::vector<haplotype*> recalcuated;
    std::vector<haplotype *> *epps;
    /* if cached, just take that */
    if (this->epp_positions_cache[read_index].size() == (size_t)parsimony_multiplicity[read_index])
    {
        epps = &this->epp_positions_cache[read_index];
    }
    else
    {
        /* recalculate on (almost) entire tree, generally the most optimal */
        int max_val = INT32_MAX;
        single_read_tree(arena, arena.reads()[read_index], recalcuated, max_val);
        std::sort(recalcuated.begin(), recalcuated.end());
        epps = &recalcuated;
    }

    if (epps->empty()) {
        return;
    }

    /* use ORIGINAL size */
    double const delta = node_score(max_parismony[read_index], parsimony_multiplicity[read_index], arena.reads()[read_index].degree);

    // randomly break epps into separate chunks; one mutex for any given chunk
    std::vector<std::pair<size_t, size_t>> chunks;
    size_t i = 0; 
    size_t c = 0;
    size_t first = (*epps)[0] - &arena.haplotypes()[0];
    for (size_t q = 0; q < epps->size(); ++q) {
        size_t const j = (*epps)[q] - &arena.haplotypes()[0];
        if (first / MUTEX_BIN_SIZE != j / MUTEX_BIN_SIZE) {
            chunks.emplace_back(i, c); 
            i = q;
            first = j;
            c = 0;
        }
        ++c;
    }
    chunks.emplace_back(i, c);

    // rotate "randomly" so threads all dont go exactly small to large
    int const random_val = ((*epps)[epps->size() / 2] - &arena.haplotypes()[0]) % chunks.size();
    std::rotate(chunks.begin(), chunks.begin() + random_val, chunks.end());

    for (const auto& [i, c] : chunks) {
        size_t j = ((*epps)[i] - &arena.haplotypes()[0]) / MUTEX_BIN_SIZE;
        mutex_t::scoped_lock lock{mutex[j]};
        
        for (size_t k = i; k < i + c; ++k) {
            haplotype* hap = (*epps)[k];
            hap->score -= delta;
        }
    }
}

void
wepp_filter::singular_step(arena& arena, haplotype* hap)
{
    std::vector<int> correspondents = find_correspondents(arena, hap);

    /* 2. map all said reads onto entire remaining set, adjusting deltas */
    std::vector<tbb::queuing_mutex> mutexes((arena.haplotypes().size() + MUTEX_BIN_SIZE - 1) / MUTEX_BIN_SIZE);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, correspondents.size()),
                      [&](tbb::blocked_range<size_t> range)
                      {
                          for (size_t i = range.begin(); i < range.end(); ++i)
                          {
                              remove_read(arena, correspondents[i], mutexes);
                          }
                      });

    /* 3. mark correspondence */
    for (int correspondent : correspondents)
    {
        this->remaining_reads.erase(correspondent);
    }
}

void
wepp_filter::clear_neighbors(arena& arena, const std::vector<haplotype*>& consideration, std::set<haplotype*>& peaks)
{
    /* 1. add to peaks */    
    peaks.insert(consideration.begin(), consideration.end());

    /* 2. for all nodes within small radius, forcefuly remove all nodes*/
    for (haplotype *pivot : consideration)
    {
        std::set<haplotype*, score_comparator> 
            multisource_radius = arena.highest_scoring_neighbors(pivot, false, MAX_PEAK_PEAK_MUTATION, INT_MAX);

        for (haplotype *node : multisource_radius)
        {
            node->mapped = true;
        }
    }
}

bool
wepp_filter::step(arena& arena, std::vector<haplotype*>& current, std::set<haplotype*> &peaks)
{
    assert(!current.empty() && !remaining_reads.empty());

    /* to get exactly same behavior as before, have a single considering set */
    std::vector<haplotype *> consideration;

    /* get TOP_N (under some special constraints) */
    auto it = current.begin();
    double const min_score = (*it)->full_score();

    /* no available peaks */
    if (min_score < SCORE_EPSILON)
    {
        return true;
    }

    /* while: */
    /* 1. we have current_nodes left */
    /* 2. the score is all the samae */
    /* 3. we have not selected 25 yet */
    /* 4. we have not exceeded max peaks */
    for (int i = 0;
         it != current.end() &&
         abs((*it)->full_score() - min_score) < SCORE_EPSILON &&
         consideration.size() < (size_t) TOP_N &&
         consideration.size() + peaks.size() < (size_t) MAX_PEAKS;
         ++i, ++it)
    {
        bool valid = true;
        for (haplotype *old : consideration)
        {
            if (!valid_two_tops(old, *it))
            {
                valid = false;
            }
        }

        if (valid)
        {
            consideration.push_back(*it);
            (*it)->mapped = true;
            //printf("* raw %.9f divergence %.9f id %s\n", (*it)->score, (*it)->dist_divergence, (*it)->id.c_str());
        }
    }

    /* clear neighbors of selected peaks (marking them as mapped as well) */
    this->clear_neighbors(arena, consideration, peaks);

    for (haplotype *&node : consideration)
    {
        this->singular_step(arena, node);
    }

    /* resort current_nodes (removing mapped nodes as well) */
    current.erase(
        std::remove_if(current.begin(), current.end(),
                       [](haplotype* curr)
                       {
                           return curr->mapped || curr->score <= SCORE_EPSILON;
                       }),
        current.end());
    tbb::parallel_sort(current.begin(), current.end(), score_comparator());

    return peaks.size() >= (size_t) MAX_PEAKS || remaining_reads.empty() || current.empty();
}

std::vector<haplotype*> 
wepp_filter::filter(arena& arena)
{
    arena.reset_haplotype_state();
    arena.build_range_trees();
    reset(arena.reads().size());

    std::vector<haplotype*> initial = arena.haplotype_pointers();

    // cartesian map
    cartesian_map(arena, initial, arena.reads());

    // iterative removal 
    std::set<haplotype*> peaks, nbrs;
    while (!step(arena, initial, peaks)) { }
    
    std::vector<haplotype*> res(peaks.begin(), peaks.end());

    // adding neighbors of peaks
    for (int k = 0; k < 5; k++) {
        std::set<haplotype*> curr_neighbors;
        arena.recover_haplotype_state();

        for (haplotype *pivot : peaks)
        {
            std::set<haplotype*, score_comparator> 
                multisource_radius = arena.highest_scoring_neighbors(pivot, false, (MAX_PEAK_PEAK_MUTATION + k), INT_MAX);

            int i = 0;
            for (haplotype *node : multisource_radius)
            {
                bool const found = peaks.find(node) != peaks.end() || curr_neighbors.find(node) != curr_neighbors.end();
                if (!found)
                {
                    node->mapped = true;
                    curr_neighbors.insert(node);
                    if (++i == MAX_NEIGHBORS_WEPP)
                    {
                        break;
                    }
                }
            }
        }

        if (abs(FREYJA_PEAKS_LIMIT - ((int)curr_neighbors.size() + (int)peaks.size())) < abs(FREYJA_PEAKS_LIMIT - ((int)nbrs.size() + (int)peaks.size()))) {
            nbrs = curr_neighbors;
        }
    }

    res.insert(res.end(), nbrs.begin(), nbrs.end());
    return res;
}

std::vector<haplotype*> 
lineage_root_filter::filter(arena& arena)
{
    auto clade_idx = arena.clade_idx();
    arena.reset_haplotype_state();
    std::vector<haplotype*> initial = arena.haplotype_pointers();
    initial.erase(
        std::remove_if(
            initial.begin(),
            initial.end(),
            [&](haplotype* hap) {
                bool is_root = false;
                for (MAT::Node* node: arena.source_nodes(hap)) {
                    if (node->clade_annotations[clade_idx] != "" && 
                        node->clade_annotations[clade_idx].rfind("misc", 0) == std::string::npos &&
                        node->clade_annotations[clade_idx].rfind("proposed", 0) == std::string::npos) {
                        is_root = true;
                        break;
                    }
                }
                return !is_root;
            }
        ),
        initial.end()
    );

    return initial;
}

std::vector<haplotype*> 
lineage_root_except_filter::filter(arena& arena) {
    arena.reset_haplotype_state();
    lineage_root_filter filter;
    std::vector<haplotype*> all_lineages = filter.filter(arena);

    haplotype* ref = arena.haplotype_with_id(this->removed_id);

    all_lineages.erase(
        std::remove_if(
            all_lineages.begin(), all_lineages.end(),
            [&](haplotype* hap) {
                return ref->mutation_distance(hap) <= this->clearance;
            }
        ),
        all_lineages.end()
    );

    return all_lineages;
}

std::vector<haplotype*> 
random_nodes_filter::filter(arena& arena)
{
    arena.reset_haplotype_state();
    std::vector<haplotype*> initial = arena.haplotype_pointers();
    // ensure root is always selected
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);
    std::shuffle(initial.begin() + 1, initial.end(), rng);
    initial.erase(initial.begin() + n, initial.end());

    std::vector<haplotype*> other = arena.haplotype_pointers();
    std::shuffle(other.begin(), other.end(), rng);
    other.erase(other.begin() + 10, other.end());

    double average_dist = 0;
    for (haplotype* hap: other) {
        int min_dist = INT_MAX;
        for (haplotype* hap2: initial) {
            min_dist = std::min(min_dist, hap->mutation_distance(hap2));
        }

        average_dist += (double) min_dist / other.size();
    }

    std::cout << "distance " << average_dist << std::endl;

    return initial;
}


void
uniform_nodes_filter::dfs(haplotype* curr, haplotype* last_set, std::vector<haplotype*>& dump)
{
    if (curr->mutation_distance(last_set) >= this->dist) {
        dump.push_back(curr);
        last_set = curr;
    }
    for (haplotype* child: curr->children) {
        dfs(child, last_set, dump);
    }
}
std::vector<haplotype*> 
uniform_nodes_filter::filter(arena& arena)
{
    arena.reset_haplotype_state();
    std::vector<haplotype*> dump;
    dump.push_back(&arena.haplotypes()[0]);
    dfs(&arena.haplotypes()[0], &arena.haplotypes()[0], dump);

    return dump;
}