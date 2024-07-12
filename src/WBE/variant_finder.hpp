#pragma once

#include <vector>
#include <utility>

#include "arena.hpp"
#include "haplotype.hpp"

class variant_finder {
public:
    // frequency of a subset for us to include it
    double freq_cutoff = 0.1;
    size_t k = 10;

    void find(arena& arena, const std::vector<std::pair<haplotype*, double>>& selected);
};

