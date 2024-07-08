#include <algorithm>

#include "initial_filter.hpp"

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
single_read_tree(arena& arena, multi_haplotype *curr, const raw_read& read, std::vector<multi_haplotype *> &max_nodes, int &max_val)
{
    // positions where this node differs from the read */
    // std::vector<int> const my_locations = curr->mutations(read.mutations, read.start, read.end);

    /* map self */
    // this basically ensures semantics are the exact same
    // as the previous algorithm
    // int parsimony, reductions = mutation_reductions(parent_locations, my_locations);
    // if (reductions)
    // {
    //     parsimony = parent_locations.size() - reductions;
    // }
    // else
    // {
    //     parsimony = my_locations.size();
    // }

    int const parsimony = curr->mutation_distance(read.mutations, read.start, read.end);
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
        single_read_tree(arena, &arena.ranged_haplotypes()[child], read, max_nodes, max_val);
    }
}

/* maps a single read to entire tree, caching aliveness */
/* and finding the epp positions of a given read */
/* parent_locations is the positions where the parent has different mutations than the read */
/* max_val corresponds to max parismony */
/* max indices correspond to the indices of the epp nodes */
static void 
single_read_tree(arena& arena, const raw_read& read, std::set<haplotype*> &max_indices, int &max_val)
{
    std::vector<multi_haplotype *> max_nodes;
    multi_haplotype* root = arena.find_range_tree_for(read);
    single_read_tree(arena, root, read, max_nodes, max_val);

    for (multi_haplotype *rnode : max_nodes)
    {
        for (haplotype *src : rnode->sources)
        {
            max_indices.emplace(src);
        }
    }
}

// calculate initial scores as well as divergences
// haps assumed to be arena indices
void 
wepp_filter::cartesian_map(arena& arena, std::vector<haplotype*>& haps, const std::vector<raw_read>& reads)
{
    std::vector<tbb::queuing_mutex> my_mutex(this->num_mutexes);

    int bin_size = arena.genome_size() / NUM_RANGE_BINS;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                      [&](tbb::blocked_range<size_t> k)
                      {
                          std::vector<std::pair<double, std::array<int, NUM_RANGE_BINS>>> score;
                          if (this->high_memory_cartesian_map) {
                              score.resize(haps.size()); // otherwise, not used and buffered directly
                          }

                          for (size_t r = k.begin(); r != k.end(); ++r)
                          {
                              std::set<haplotype *> max_indices;
                              int max_val = INT32_MAX; /* really want minimum value */

                              single_read_tree(arena, reads[r], max_indices, max_val);

                              double const delta = node_score(max_val, max_indices.size(), reads[r].degree);
                              {
                                  int bucket = reads[r].start / bin_size;
                                  if (this->high_memory_cartesian_map) {
                                      haplotype *arena_base = &arena.haplotypes()[0];
                                      for (haplotype *hap : max_indices)
                                      {
                                          int index = (hap - arena_base);
                                          score[index].first += delta;
                                          score[index].second[bucket] += reads[r].degree;
                                      }
                                  }
                                  else {
                                      for (haplotype *hap : max_indices)
                                      {
                                          int index = (hap - arena_base);
                                          int mutex_bucket = index % this->num_mutexes;
                                          tbb::queuing_mutex::scoped_lock my_lock{my_mutex[mutex_bucket]};
                                          hap->score += delta;
                                          hap->mapped_read_counts[bucket] += reads[r].degree;
                                      }
                                  }
                              }

                              this->max_parismony[r] = max_val;
                              this->parsimony_multiplicity[r] = max_indices.size();
                              if (max_indices.size() <= max_cached_epp_size)
                              {
                                  this->epp_positions_cache[r] = std::move(max_indices);
                              }
                          }

                          if (this->high_memory_cartesian_map)
                          {
                              // only use a global mutex in this case
                              tbb::queuing_mutex::scoped_lock my_lock{my_mutex[0]};
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
        double divergence = 0;
        for (size_t j = 0; j < NUM_RANGE_BINS; ++j)
        {
            double const proportion = (double) haps[i]->mapped_read_counts[j] / arena.read_distribution()[j];
            if (proportion > read_dist_factor_threshold)
            {
                divergence += 1.0 / NUM_RANGE_BINS;
            }
        }

        haps[i]->dist_divergence = divergence;
    }

    /* initial sort into scores */
    std::sort(haps.begin(), haps.end(), score_comparator());
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
                              if (epp_positions_cache[read].find(hap) != epp_positions_cache[read].end())
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
                                  const raw_read& r = arena.reads()[read];
                                  //   std::vector<int> parent = hap->parent ? hap->parent->mutations(r.mutations, r.start, r.end) : std::vector<int>();
                                  //   std::vector<int> ours = hap->mutations(r.mutations, r.start, r.end);

                                  //   int parsimony, reductions = mutation_reductions(parent, ours);
                                  //   if (reductions)
                                  //   {
                                  //       parsimony = parent.size() - reductions;
                                  //   }
                                  //   else
                                  //   {
                                  //       parsimony = ours.size();
                                  //   }

                                  int parsimony = hap->mutation_distance(r);
                                  if (parsimony == max_parismony[read])
                                  {
                                      tbb::queuing_mutex::scoped_lock lock{my_mutex};
                                      correspondents.push_back(read);
                                  }
                              }
                          }
                      });
    std::sort(correspondents.begin(), correspondents.end());

    return correspondents;
}

void
wepp_filter::remove_read(arena& arena, int read_index, tbb::queuing_mutex* mutex)
{
    using mutex_t = tbb::queuing_mutex;

    std::set<haplotype*> recalcuated;
    std::set<haplotype *> *epps;
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
        epps = &recalcuated;
    }

    /* use ORIGINAL size */
    mutex_t::scoped_lock lock{*mutex};
    double const delta = node_score(max_parismony[read_index], parsimony_multiplicity[read_index], arena.reads()[read_index].degree);
    for (haplotype* hap: (*epps))
    {
        if (hap->mapped)
        {
            continue;
        }

        hap->score -= delta;
    }
}

void
wepp_filter::singular_step(arena& arena, haplotype* hap)
{
    std::vector<int> correspondents = find_correspondents(arena, hap);

    /* 2. map all said reads onto entire remaining set, adjusting deltas */
    tbb::queuing_mutex my_mutex;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, correspondents.size()),
                      [&](tbb::blocked_range<size_t> range)
                      {
                          for (size_t i = range.begin(); i < range.end(); ++i)
                          {
                              remove_read(arena, correspondents[i], &my_mutex);
                          }
                      });

    /* 3. mark correspondence */
    for (int correspondent : correspondents)
    {
        this->remaining_reads.erase(correspondent);
    }
}

void
wepp_filter::clear_neighbors(arena& arena, const std::vector<haplotype*>& consideration, std::set<haplotype*, mutation_comparator>& peaks, std::set<haplotype*, mutation_comparator>& nbrs)
{
    /* 1. add to peaks */    
    peaks.insert(consideration.begin(), consideration.end());

    std::vector<haplotype*> added;
    /* 2. find candidates within radius (taking only top n for each pivot) */
    /* buffer the selected neighbors into added */
    for (haplotype *pivot : consideration)
    {
        std::set<haplotype*, score_comparator> 
            multisource_radius = arena.highest_scoring_neighbors(pivot, false, max_peak_nonpeak_mutation, INT_MAX);

        int i = 0;
        for (haplotype *node : multisource_radius)
        {
            bool const found = peaks.find(node) != peaks.end() || nbrs.find(node) != nbrs.end();
            if (!found)
            {
                node->mapped = true;
                added.emplace_back(node);
                if (++i == max_neighbors)
                {
                    break;
                }
            }
        }
    }

    /* 3. for all nodes within small radius, forcefuly remove even if more than 100 */
    for (haplotype *pivot : consideration)
    {
        std::set<haplotype*, score_comparator> 
            multisource_radius = arena.highest_scoring_neighbors(pivot, false, max_peak_peak_mutation, INT_MAX);

        for (haplotype *node : multisource_radius)
        {
            node->mapped = true;
        }
    }

    nbrs.insert(added.begin(), added.end());
}

bool
wepp_filter::step(arena& arena, std::vector<haplotype*>& current, std::set<haplotype*, mutation_comparator> &peaks, std::set<haplotype*, mutation_comparator> &nbrs)
{
    assert(!current.empty() && !remaining_reads.empty());

    /* to get exactly same behavior as before, have a single considering set */
    std::vector<haplotype *> consideration;

    /* get top_n (under some special constraints) */
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
         consideration.size() < (size_t) top_n &&
         consideration.size() + peaks.size() < (size_t) max_peaks;
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
            // printf("%.9f raw %.9f divergence id: %s\n", (*it)->score, (*it)->dist_divergence, (*it)->id.c_str());
        }
    }

    /* clear neighbors of selected peaks (marking them as mapped as well) */
    this->clear_neighbors(arena, consideration, peaks, nbrs);

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
    std::sort(current.begin(), current.end(), score_comparator());

    return peaks.size() >= (size_t) max_peaks || remaining_reads.empty() || current.empty();
}

std::vector<haplotype*> 
wepp_filter::filter(arena& arena)
{
    arena.reset_haplotype_state();
    reset(arena.reads().size());

    std::vector<haplotype*> initial = arena.haplotype_pointers();

    // cartesian map
    cartesian_map(arena, initial, arena.reads());

    // iterative removal 
    std::set<haplotype*, mutation_comparator> peaks, nbrs;
    while (!step(arena, initial, peaks, nbrs)) { }

    std::vector<haplotype*> res(peaks.begin(), peaks.end());
    res.insert(res.end(), nbrs.begin(), nbrs.end());

    return res;
}

std::vector<haplotype*> 
lineage_root_filter::filter(arena& arena)
{
    std::vector<haplotype*> initial = arena.haplotype_pointers();
    initial.erase(
        std::remove_if(
            initial.begin(),
            initial.end(),
            [&](haplotype* hap) {
                bool is_root = false;
                for (MAT::Node* node: arena.source_nodes(hap)) {
                    if (node->clade_annotations[1] != "" && 
                        node->clade_annotations[1].rfind("misc", 0) == std::string::npos &&
                        node->clade_annotations[1].rfind("proposed", 0) == std::string::npos) {
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