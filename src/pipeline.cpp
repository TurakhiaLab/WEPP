#include <memory> 

#include <tbb/blocked_range.h>
#include <tbb/queuing_mutex.h>
#include <tbb/parallel_for.h>

#include "timer.hpp"

#include "pipeline.hpp"

void detect_peaks(const dataset& d, bool is_initial) {
    auto main = std::make_unique<wepp_filter>();
    //auto main = std::make_unique<lineage_root_filter>();
    
    auto post = std::make_unique<freyja_post_filter>();
    //post->num_filter_rounds = 1;
    post->num_filter_rounds = 10;

    pipeline p{d, std::move(main), std::move(post)};
    if (is_initial) {
        p.run_initial();
    }
    else {
        p.run_from_last_initial(false);
    }
}

void pipeline::run_initial() {
    std::vector<haplotype *> running = a.haplotype_pointers();

    {
        std::cout << "----- [running initial filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        timer t;
        running = main->filter(a);
        // a.print_mutation_distance(running);
        std::cout << "--- initial filter took " << t.seconds() << " seconds " << std::endl << std::endl;
    }

    // save peaks
    save(running, ds.first_checkpoint_path());
}

void pipeline::run_from_last_initial(bool is_full_run) {
    std::vector<haplotype *> running = this->recover(ds.first_checkpoint_path()); 

    {
        std::cout << "----- [running post filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        timer t; 
        std::vector<std::pair<haplotype*, double>> full = post->iterative_filter(a, running);
        a.print_full_report(full);

        std::cout << "--- post filter took " << t.seconds() << " seconds " << std::endl;

        a.dump_haplotype_proportion(full);
        a.dump_read2node_mapping(full);
    }
}