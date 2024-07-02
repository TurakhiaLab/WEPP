#pragma once

#include <vector>

#include "arena.hpp"
#include "haplotype.hpp"

class initial_filter {

    virtual std::vector<haplotype> filter(const arena& arena, std::vector<haplotype> initial) = 0;
};

class wepp_filter: public initial_filter {
    std::vector<haplotype> filter(const arena& arena, std::vector<haplotype> initial);
};