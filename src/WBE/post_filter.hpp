#pragma once

#include <vector>
#include <cmath>

#include "haplotype.hpp"
#include "arena.hpp"

class post_filter {
public:
    virtual std::vector<haplotype*> filter(arena& arena, std::vector<haplotype*> input) = 0;

    std::vector<haplotype*> iterative_filter(arena& arena, std::vector<haplotype*> input, int rounds) {
        for (int i = 0; i < rounds; ++i) {
            input = this->filter(arena, input);
            // TODO add neighbors at some point
        }
        return input;
    }
};

class freyja_post_filter: public post_filter {
    void dump_barcode(arena& a, const std::vector<haplotype*>& inputs);
public:
    std::vector<haplotype*> filter(arena& arena, std::vector<haplotype*> input) override;
};


class likelihood_post_filter: public post_filter {
    double sigma = 4.0;
    double coeff = log(1 / (sigma * sqrt(2 * M_PI)));
    size_t max_peaks = 150;
public:
    std::vector<haplotype*> filter(arena& arena, std::vector<haplotype*> input) override;
};

class em_post_filter: public post_filter {
    double alpha = 0.005;
    double epsilon = 1e-3;
    double min_proportion = 1 / 200.0;
    int max_it = 50;
public:
    std::vector<haplotype*> filter(arena& arena, std::vector<haplotype*> input) override;
};

class kmeans_post_filter: public post_filter {
    int max_it = 15, explore_rad = 5;
    double stay_factor = 0.75;
    double min_abundance = 1 / 200.0;
public:
    std::vector<haplotype *> filter(arena &arena, std::vector<haplotype *> input) override;
};