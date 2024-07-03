#pragma once

#include <vector>
#include <iostream>
#include <memory>

#include "dataset.hpp"
#include "arena.hpp"
#include "initial_filter.hpp"
#include "post_filter.hpp"

class pipeline {
    const dataset& ds;
    arena a;
    std::unique_ptr<initial_filter> main;
    std::unique_ptr<post_filter> post;

public:
    pipeline(const dataset& ds, std::unique_ptr<initial_filter>&& main, std::unique_ptr<post_filter>&& post) 
        : ds{ds}, a{ds}, main{std::move(main)}, post{std::move(post)}
    {
    }

    void run();
};

void detect_peaks(const dataset& d);
