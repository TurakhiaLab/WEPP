#include <memory> 

#include "pipeline.hpp"

void detect_peaks(const dataset& d) {
    std::unique_ptr<initial_filter> main;
    main = std::make_unique<em_filter>();
    auto post = std::make_unique<em_post_filter>();
    post->num_filter_rounds = MAX_NEIGHBOR_ITERATIONS;

    pipeline p{d, std::move(main), std::move(post)};
    
    if (FULL_RUN)
        p.run();
    else
        p.run_from_last_initial(false);
}

void pipeline::run() {
    std::vector<haplotype *> running = a.haplotype_pointers();

    Timer timer;
    {
        std::cout << "----- [running initial filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        timer.Start();
        running = main->filter(a);
        std::cout << "--- initial filter took " << timer.Stop() / 1000 << " seconds " << std::endl << std::endl;
    }

    // save peaks
    save(running, ds.checkpoint_path());

    this->run_from_last_initial(true);
}

void pipeline::run_from_last_initial(bool is_full_run) {
    std::vector<haplotype *> running = this->recover(ds.checkpoint_path());
    {
        std::cout << "----- [running post filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        Timer timer; timer.Start();
        
        if (is_full_run)
            a.recover_haplotype_state();
        else 
            a.reset_haplotype_state();

        std::vector<std::pair<haplotype*, double>> full = post->iterative_filter(a, running);
        
        a.print_full_report(full);

        if (SIMULATED_DATA)
        {
            std::vector<haplotype*> current;
            current.reserve(full.size());
            std::transform(full.begin(), full.end(), std::back_inserter(current),
                       [](const std::pair<haplotype*, double>& pair) {
                           return pair.first;
                       });
            a.print_mutation_distance(current);
        }

        std::cout << "--- post filter took " << timer.Stop() / 1000 << " seconds " << std::endl;

        a.dump_haplotype_proportion(full);
        a.dump_haplotype_uncertainty(full);
        a.dump_lineage_proportion(full);
        //a.resolve_unaccounted_mutations(full);
        //a.dump_haplotypes(full);
        //a.dump_read2haplotype_mapping(full);

        std::cout << "--- RUN COMPLETED" << std::endl;
    }
}