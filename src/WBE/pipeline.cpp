#include <memory> 

#include "pipeline.hpp"

void detect_peaks(const dataset& d) {
    auto main = std::make_unique<wepp_filter>();
    auto post = std::make_unique<likelihood_post_filter>();

    pipeline p{d, std::move(main), std::move(post)};
    // p.run();
    p.run_from_last_initial();
}

void pipeline::run() {
    tbb::task_scheduler_init init(this->ds.num_threads());

    std::vector<haplotype *> running = a.haplotype_pointers();

    {
        std::cout << "----- [running main filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        Timer timer; timer.Start();
        running = main->filter(a);
        a.print_mutation_distance(running);
        std::cout << "--- main filter took " << timer.Stop() / 1000 << " seconds " << std::endl << std::endl;
    }

    // save peaks
    save(running, ds.first_checkpoint_path());

    {
        std::cout << "----- [running post filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        timer.Start();
        running = post->filter(a, running);
        a.print_mutation_distance(running);
        std::cout << "--- post filter took " << timer.Stop() / 1000 << " seconds " << std::endl;
    }

    save(running, ds.last_checkpoint_path());
}

void pipeline::run_from_last_initial() {
    tbb::task_scheduler_init init(this->ds.num_threads());

    std::vector<haplotype *> running = this->recover(ds.checkpoint_path()); 

    {
        std::cout << "----- [running post filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        Timer timer; timer.Start();
        running = post->filter(a, running);
        a.print_mutation_distance(running);
        std::cout << "--- post filter took " << timer.Stop() / 1000 << " seconds " << std::endl;
    }

    save(running, ds.last_checkpoint_path());
}
