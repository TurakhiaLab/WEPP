#pragma once

#include <vector>
#include <cmath>

#include "haplotype.hpp"
#include "arena.hpp"
#include "config.hpp"

class post_filter {
public:
    int num_filter_rounds = 1;
    int freeze_round = 1;
    int max_nbrs = MAX_NEIGHBORS_FREYJA;
    int max_rad = MAX_NEIGHBOR_MUTATION;

    virtual std::vector<std::pair<haplotype*, double>> filter(arena& arena, std::vector<haplotype*> input) = 0;

    std::vector<std::pair<haplotype*, double>> iterative_filter(arena& arena, std::vector<haplotype*> input) {
        assert(num_filter_rounds >= 1);

        std::set<haplotype*> frozen, last_round;
        for (int i = 0; i < num_filter_rounds; ++i) {
            std::vector<haplotype *> full_input = input;
            for (auto hap: frozen) {
                if (std::find(full_input.begin(), full_input.end(), hap) == full_input.end())
                    full_input.emplace_back(hap);
            }

            std::vector<std::pair<haplotype*, double>> filtered = this->filter(arena, full_input); 

            std::vector<haplotype *> this_round;
            std::transform(filtered.begin(), filtered.end(), std::back_inserter(this_round), 
                [](const auto& p) {
                    return p.first;
                }
            );

            //arena.print_full_report(filtered);
            arena.print_mutation_distance(this_round);

            std::sort(this_round.begin(), this_round.end());
            if (i == num_filter_rounds - 1 || 
                std::includes(last_round.begin(), last_round.end(), this_round.begin(), this_round.end())) {
                return filtered;
            }


            // common should go into frozen
            if (i >= freeze_round) {
                std::set_intersection(this_round.begin(), this_round.end(), last_round.begin(), last_round.end(),
                                      std::inserter(frozen, frozen.end()));
            }
            last_round = std::set<haplotype *>(this_round.begin(), this_round.end());

            // add neighbors
            {
                std::set<haplotype *, score_comparator> build;
                for (haplotype* hap: this_round) {
                    std::set<haplotype *, score_comparator> nbrs = arena.closest_neighbors(hap, max_rad, max_nbrs);
                    build.insert(nbrs.begin(), nbrs.end());
                }
                input = std::vector<haplotype*>(build.begin(), build.end());
            }
        }
        
        return {};
    }

    virtual ~post_filter() { };
};

class identity_post_filter: public post_filter {
public:
    std::vector<std::pair<haplotype*, double>> filter(arena& arena, std::vector<haplotype*> input) override {
        std::vector<std::pair<haplotype*, double>> ret;
        std::transform(input.begin(), input.end(), std::back_inserter(ret), 
            [&](haplotype* in)  {
                return std::make_pair(in, 1.0 / input.size());
            }
        );

        return ret;
    }
};


class freyja_post_filter: public post_filter {
    void dump_barcode(arena& a, const std::vector<haplotype*>& inputs);
public:
    std::vector<std::pair<haplotype*, double>> filter(arena& arena, std::vector<haplotype*> input) override;
};


class likelihood_post_filter: public post_filter {
public:
    double sigma = 4.0;
    double coeff = log(1 / (sigma * sqrt(2 * M_PI)));
    size_t max_peaks = 150;

    std::vector<std::pair<haplotype*, double>> filter(arena& arena, std::vector<haplotype*> input) override;
};

class em_post_filter: public post_filter {
public:
    double alpha = 0.005;
    double epsilon = 1e-3;
    double min_proportion = 1 / 200.0;
    int max_it = 50;

    std::vector<std::pair<haplotype*, double>> filter(arena& arena, std::vector<haplotype*> input) override;
};

class kmeans_post_filter: public post_filter {
public:
    int max_it = 15, explore_rad = 5;
    double stay_factor = 0.75;
    double min_abundance = 1 / 200.0;

    std::vector<std::pair<haplotype *, double>> filter(arena &arena, std::vector<haplotype *> input) override;
};