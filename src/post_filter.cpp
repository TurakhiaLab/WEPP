#include "post_filter.hpp"

#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>
#include <tbb/parallel_for.h>

#include <numeric>
#include <random>

std::vector<std::pair<haplotype*, double>>
likelihood_post_filter::filter(arena& arena, std::vector<haplotype*> input)
{
    arena.reset_haplotype_state();

    const std::vector<raw_read>& reads = arena.reads();

    std::vector<std::vector<double>> parsimony(reads.size(), std::vector<double>(input.size()));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                      [&](tbb::blocked_range<size_t> range)
                      {
                          for (size_t i = range.begin(); i < range.end(); ++i)
                          {
                              for (size_t j = 0; j < reads.size(); ++j)
                              {
                                  double dist = input[i]->mutation_distance(reads[j]) / sigma;
                                  parsimony[j][i] = reads[j].degree * (coeff - dist * dist);
                                  input[i]->score += parsimony[j][i];
                              }
                          }
                      });

    std::vector<haplotype *> consideration, current_nodes = input;
    for (size_t q = 0; q < input.size() && q < max_peaks; ++q)
    {
        std::sort(current_nodes.begin(), current_nodes.end(), score_comparator());

        current_nodes[0]->mapped = true;
        consideration.emplace_back(current_nodes[0]);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t i = range.begin(); i < range.end(); ++i)
                              {
                                  input[i]->score = 0;
                              }
                          });

        std::vector<double> comp(reads.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t i = range.begin(); i < range.end(); ++i)
                              {
                                  double dist = current_nodes[0]->mutation_distance(reads[i]) / sigma;
                                  comp[i] = reads[i].degree * (coeff - dist * dist);
                              }
                          });

        tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t i = range.begin(); i < range.end(); ++i)
                              {
                                  for (size_t j = 0; j < reads.size(); ++j)
                                  {
                                      parsimony[j][i] = std::max(parsimony[j][i], comp[j]);
                                      input[i]->score += parsimony[j][i];
                                  }
                              }
                          });

        current_nodes.erase(
            std::remove_if(
                current_nodes.begin(),
                current_nodes.end(),
                [](haplotype*node) -> bool
                {
                    return node->mapped;
                }),
            current_nodes.end());
    }

    std::vector<std::pair<haplotype *, double>> abundance;
    std::transform(consideration.begin(), consideration.end(), std::back_inserter(abundance), 
        [&](haplotype* hap) {
            return std::make_pair(hap, 1.0 / consideration.size());
        }
    );

    return abundance;
}

void 
freyja_post_filter::dump_barcode(arena& a, const std::vector<haplotype*>& haplotypes) 
{
    std::set<std::string> mutations;
    std::vector<std::set<std::string>> node_muts;

    std::map<size_t, size_t> ungapped_pos;
    size_t non_gap_reference_size = 0;
    for (size_t i = 0; i < a.reference().size(); ++i) {
        if (a.reference()[i] != '_')
        {
            ungapped_pos[i + 1] = ++non_gap_reference_size;
        }
    }

    auto ungapped_string = [&](mutation const& m) {
        char ref = char_from_nuc(m.ref);
        int pos = ungapped_pos.at(m.pos);
        char mut = char_from_nuc(m.mut);
        return ref + std::to_string(pos) + mut;
    };

    for (haplotype *n : haplotypes)
    {
        std::set<std::string> my_muts;
        for (const mutation& m : n->stack_muts())
        {
            if (m.mut != NUC_N && !m.is_indel()) {
                std::string build = ungapped_string(m);
                mutations.insert(build);
                my_muts.insert(std::move(build));
            }
        }
        node_muts.push_back(std::move(my_muts));
    }

    for (const raw_read& read: a.reads()) {
        for (mutation mut : read.mutations)
        {
            if (mut.mut != NUC_N && !mut.is_indel()) {
                mutations.insert(ungapped_string(mut));
            }
        }
    }

    std::vector<std::string> mutation_vec(mutations.begin(), mutations.end());

    std::ofstream outfile("Freyja/data/usher_barcodes.csv");
    for (const std::string &mut : mutations)
    {
        outfile << "," << mut;
    }
    outfile << std::endl;
    for (size_t i = 0; i < haplotypes.size(); ++i)
    {
        haplotype *n = haplotypes[i];
        outfile << 'N' << (n - &a.haplotypes()[0]);
        for (const std::string &mut : mutations)
        {
            bool contains = node_muts[i].find(mut) != node_muts[i].end();
            outfile << "," << (contains ? '1' : '0');
        }
        outfile << std::endl;
    }
}

std::vector<std::pair<haplotype*, double>>
freyja_post_filter::filter(arena& arena, std::vector<haplotype*> input)
{
    dump_barcode(arena, input);

    if (std::system(
            "bash -c \""
                "source ~/miniconda3/etc/profile.d/conda.sh && "
                "conda activate freyja-env && "
                "freyja demix Freyja/cwap_variants.tsv Freyja/cwap_depth.tsv --barcodes Freyja/data/usher_barcodes.csv --output Freyja/my_output_latest.txt --eps 0.001"
            "\""
        ) != 0)
    {
        std::cerr << "Failed to run freyja" << std::endl;
        return {};
    }

    std::vector<std::pair<haplotype *, double>> freyja_nodes;

    std::ifstream fin("Freyja/my_output_latest.txt");
    std::string tmp;
    std::getline(fin, tmp);
    std::getline(fin, tmp);
    fin >> tmp;
    // all of the selected ids
    std::getline(fin, tmp);
    std::stringstream ss{tmp};
    std::string index;
    while (ss >> index)
    {
        // skip past the 'N'
        int ind = atoi(index.c_str() + 1);
        freyja_nodes.emplace_back(&arena.haplotypes()[ind], 0);
    }

    fin >> tmp;
    std::getline(fin, tmp);
    ss = std::stringstream{tmp};
    double sum = 0;
    for (size_t i = 0; i < freyja_nodes.size(); ++i) {
        double abundance; ss >> abundance;
        freyja_nodes[i].second = abundance;
        sum += abundance;
    }

    for (size_t i = 0; i < freyja_nodes.size(); ++i) {
        freyja_nodes[i].second /= sum;
    }

    return freyja_nodes;
}

std::vector<std::pair<haplotype*, double>>
em_post_filter::filter(arena& arena, std::vector<haplotype*> input)
{
    arena.reset_haplotype_state();

    const std::vector<raw_read>& reads = arena.reads();
    std::vector<haplotype*> subset = std::move(input);

    std::vector<std::vector<double>> q(subset.size(), std::vector<double>(reads.size()));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()),
                      [&](tbb::blocked_range<size_t> range)
                      {
                          for (size_t i = range.begin(); i < range.end(); ++i)
                          {
                              for (size_t j = 0; j < reads.size(); ++j)
                              {
                                  double dist = subset[i]->mutation_distance(reads[j]);
                                  q[i][j] = std::pow(this->alpha, dist) * std::pow(1 - this->alpha, reads[j].end - reads[j].start + 1 - dist);
                              }
                          }
                      });

    double prev = 0, curr = 0;
    std::vector<double> p(subset.size(), 1.0 / subset.size());
    do
    {
        prev = curr;

        std::vector<double> denom(reads.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t j = range.begin(); j < range.end(); ++j)
                              {
                                  double sum = 0;
                                  for (size_t l = 0; l < subset.size(); ++l)
                                  {
                                      sum += p[l] * q[l][j];
                                  }
                                  denom[j] = sum;
                              }
                          });

        tbb::parallel_for(tbb::blocked_range<size_t>(0, subset.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t i = range.begin(); i < range.end(); ++i)
                              {
                                  double sum = 0;
                                  size_t count = 0;
                                  for (size_t j = 0; j < reads.size(); ++j)
                                  {
                                      sum += reads[j].degree * p[i] * q[i][j] / denom[j];
                                      count += reads[j].degree;
                                  }
                                  // don't even need multiple p
                                  p[i] = sum / count;
                              }
                          });

        std::vector<double> logl(reads.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t j = range.begin(); j < range.end(); ++j)
                              {
                                  double sum = 0;
                                  for (size_t l = 0; l < subset.size(); ++l)
                                  {
                                      sum += p[l] * q[l][j];
                                  }
                                  logl[j] = std::log(sum);
                              }
                          });

        curr = std::accumulate(logl.begin(), logl.end(), 0.0);
    } while (abs(curr - prev) > epsilon && --max_it);

    for (size_t i = 0; i < p.size(); ++i)
    {
        subset[i]->score = p[i];
    }
    std::sort(subset.begin(), subset.end(), score_comparator());
    auto it = subset.begin();
    double total = 0;
    while (it != subset.end() && (*it)->score / ((*it)->score + total) >= min_proportion)
    {
        total += (*it)->score;
        ++it;
    }
    subset.erase(it, subset.end());

    std::vector<std::pair<haplotype *, double>> abundance;
    std::transform(subset.begin(), subset.end(), std::back_inserter(abundance), 
        [&](haplotype* hap) {
            return std::make_pair(hap, hap->score / total);
        }
    );

    return abundance;
}

std::vector<std::pair<haplotype*, double>>
kmeans_post_filter::filter(arena& arena, std::vector<haplotype*> input)
{
    using my_mutex_t = tbb::queuing_mutex;
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    // const size_t n = (size_t)(1 + 1 / min_abundance);

    arena.reset_haplotype_state();
    const std::vector<raw_read>& reads = arena.reads();

    size_t const total_reads_degree = arena.num_reads_pre_merge();

    for (int k = 0; k < max_it; ++k)
    {
        my_mutex_t mutex;
        std::vector<std::vector<int>> corresponding_reads(input.size());
        std::vector<std::pair<double, int>> index_set(input.size());
        for (size_t i = 0; i < input.size(); ++i)
        {
            index_set[i].second = i;
        }

        // map all reads to current peaks
        // and create epp sets
        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t i = range.begin(); i < range.end(); ++i)
                              {
                                  int min_dist = INT_MAX;
                                  std::vector<size_t> best;
                                  for (size_t j = 0; j < input.size(); ++j)
                                  {
                                      int dist = input[j]->mutation_distance(reads[i]);
                                      if (dist < min_dist)
                                      {
                                          min_dist = dist;
                                          best = {j};
                                      }
                                      else if (dist == min_dist)
                                      {
                                          best.emplace_back(j);
                                      }
                                  }

                                  // avoiding data race on both rng and corresponding_reads
                                  my_mutex_t::scoped_lock _lock{mutex};
                                  int ind = best[rng() % best.size()];
                                  corresponding_reads[ind].emplace_back(i);
                                  index_set[ind].first += (double)reads[ind].degree / total_reads_degree;
                              }
                          });

        std::sort(index_set.begin(), index_set.end());
        // index of previous selection
        std::vector<size_t> seeds;
        // size_t j = 0;
        // double running_abundance = 0;
        // we delete nodes with too little abundance
        // and split nodes with too high abudance
        // idea is basically that
        // < abundance / 2 = delete
        // > abundance * 2 = split (theoretically to same node)
        // while (j < index_set.size()) {
        //     if (index_set[j].first + running_abundance + 1e-6 >= min_abundance) {
        //         seeds.emplace_back(index_set[j].second);
        //         index_set[j].first -= min_abundance - running_abundance;
        //         running_abundance = 0;
        //     }
        //     else {
        //         running_abundance += index_set[j].first;
        //         index_set[j].first = 0;
        //         ++j;
        //     }
        // }
        std::vector<haplotype *> new_selection;

        // debug (ignore running abundance)
        // seeds = std::vector<size_t>(selected.size());
        // std::iota(seeds.begin(), seeds.end(), 0);
        seeds = {};
        for (size_t j = 0; j < index_set.size(); ++j)
        {
            if (index_set[j].first < min_abundance * 2)
            {
                seeds.push_back(index_set[j].second);
            }
            else
            {
                seeds.push_back(index_set[j].second);
            }
        }

        // map n -> n'
        for (size_t i = 0; i < seeds.size(); ++i)
        {
            std::set<haplotype*> 
                all_neighbors = arena.closest_neighbors(input[seeds[i]], explore_rad, INT_MAX);

            // randomly select weighted neighbor based on scores
            my_mutex_t mutex;

            std::vector<haplotype*> flat(all_neighbors.begin(), all_neighbors.end());
            std::vector<size_t> distances;
            size_t best_distance = SIZE_MAX;
            for (haplotype *node : flat)
            {
                size_t total_sum = 0;
                tbb::parallel_for(tbb::blocked_range<size_t>(0, corresponding_reads[seeds[i]].size()), [&](tbb::blocked_range<size_t> range)
                                  {
                        size_t sum = 0;
                        for (size_t j = range.begin(); j < range.end(); ++j)
                        {
                            raw_read read = reads[corresponding_reads[seeds[i]][j]];
                            sum += (size_t) read.degree * node->mutation_distance(read);
                        }

                        my_mutex_t::scoped_lock _lock{mutex};
                        total_sum += sum; 
                        });

                distances.push_back(total_sum);
                best_distance = std::min(best_distance, total_sum);
            }

            // lower average -> higher chance of selection
            std::vector<double> prob;
            double normalization = 0;
            for (size_t j = 0; j < flat.size(); ++j)
            {
                normalization += exp(-stay_factor * (double)(distances[j] - best_distance) / flat.size());
                prob.push_back(normalization);
            }
            std::uniform_real_distribution<double> unif(0, normalization);
            double rand = unif(rng);
            double prev = 0;
            haplotype *best_node = flat.back();
            for (size_t j = 0; j < prob.size(); ++j)
            {
                if (prev <= rand && rand < prob[j])
                {
                    best_node = flat[j];
                    break;
                }
                prev = prob[j];
            }
            new_selection.emplace_back(best_node);
        }

        input = new_selection;
    }

    std::vector<std::pair<haplotype *, double>> abundance;
    std::transform(input.begin(), input.end(), std::back_inserter(abundance), 
        [&](haplotype* hap) {
            return std::make_pair(hap, 1.0 / input.size());
        }
    );

    return abundance;
}