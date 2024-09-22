#pragma once

#include <vector>
#include <iostream>
#include <memory>

#include "dataset.hpp"
#include "arena.hpp"
#include "initial_filter.hpp"
#include "post_filter.hpp"

class pipeline {
public:
    const dataset& ds;
    arena a;
    std::unique_ptr<initial_filter> main;
    std::unique_ptr<post_filter> post;

    pipeline(const dataset& ds, std::unique_ptr<initial_filter>&& main, std::unique_ptr<post_filter>&& post) 
        : ds{ds}, a{ds}, main{std::move(main)}, post{std::move(post)}
    {
    }

    void save(const std::vector<haplotype*>& running, const std::string& path) {
        std::ofstream fout(path);
        for (haplotype* hap : running) {
            fout << (hap - &a.haplotypes()[0]) << std::endl;
        }
    }
    
    std::vector<haplotype*> recover(const std::string& path)
    {
        std::ifstream fin(path);
        std::vector<haplotype*> running;
        int index;
        while (fin >> index) {
            running.push_back(&a.haplotypes()[index]);
        }

        return running;
    }

    void run_initial();
    void run_from_last_initial(bool is_full_run);
};

void detect_peaks(const dataset& d, bool is_initial);