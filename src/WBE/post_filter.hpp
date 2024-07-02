#pragma once
#include <vector>

#include "haplotype.hpp"
#include "arena.hpp"

class post_filter {
public:
    virtual std::vector<haplotype*> filter(const arena& arena, std::vector<haplotype*> input) const = 0;

    std::vector<haplotype*> iterative_filter(const arena& arena, std::vector<haplotype*> input, int rounds) {
        for (int i = 0; i < rounds; ++i) {
            input = this->filter(arena, input);
            // TODO add neighbors at some point
        }
        return input;
    }
};

class freyja_post_filter: public post_filter {

};

class likelihood_post_filter: public post_filter {

};

class em_post_filter: public post_filter {

};

class kmeans_post_filter: public post_filter {

};