#include <memory> 

#include "pipeline.hpp"

void detect_peaks(const dataset& d) {
    auto main = std::make_unique<wepp_filter>();
    // main->removed_id = "hCoV-19/Hong";

    auto post = std::make_unique<freyja_post_filter>();
    post->num_filter_rounds = 10;

    pipeline p{d, std::move(main), std::move(post)};
    // p.run();
    p.run_from_last_initial();
}

void pipeline::run() {
    std::vector<haplotype *> running = a.haplotype_pointers();

    Timer timer;
    {
        std::cout << "----- [running initial filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        timer.Start();
        running = main->filter(a);
        a.print_mutation_distance(running);
        std::cout << "--- initial filter took " << timer.Stop() / 1000 << " seconds " << std::endl << std::endl;
    }

    // save peaks
    save(running, ds.first_checkpoint_path());

    this->run_from_last_initial();
}

void pipeline::run_from_last_initial() {
    std::vector<haplotype *> running = this->recover(ds.first_checkpoint_path()); 

    {
        std::cout << "----- [running post filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        Timer timer; timer.Start();
        std::vector<std::pair<haplotype*, double>> full = post->iterative_filter(a, running);
        
        a.print_full_report(full);
        std::cout << "--- post filter took " << timer.Stop() / 1000 << " seconds " << std::endl;

        vf.find(a, full);
    }
}
