#pragma once

#include <vector>
#include <memory>

#include "dataset.hpp"
#include "arena.hpp"
#include "initial_filter.hpp"
#include "post_filter.hpp"

class pipeline {
    dataset ds;
    arena arena;
    std::unique_ptr<initial_filter> main;
    std::unique_ptr<post_filter> post;

    pipeline(po::parsed_options options, std::unique_ptr<initial_filter>&& main, std::unique_ptr<post_filter>&& post) 
        ds{directory}
    {
        this->main = std::move(main);
        this->post = std::move(post);
    }

    void run() {
        std::vector<haplotype*> running = main_filter->filter();
        running = post_filter.filter(running);
    }
};

