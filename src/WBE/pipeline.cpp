#include <memory> 

#include "pipeline.hpp"

void detect_peaks(const dataset& d) {
    auto main = std::make_unique<wepp_filter>();
    //auto main = std::make_unique<lineage_root_filter>();
    
    auto post = std::make_unique<freyja_post_filter>();
    //post->num_filter_rounds = 1;
    post->num_filter_rounds = 10;

    pipeline p{d, std::move(main), std::move(post)};
    // p.a.print_cooccuring_mutations(600);
    p.run();
    //p.run_from_last_initial(false);
}

void pipeline::run() {
    std::vector<haplotype *> running = a.haplotype_pointers();

    Timer timer;
    {
        std::cout << "----- [running initial filter] -----" << std::endl;
        std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;

        timer.Start();
        running = main->filter(a);
        //a.print_mutation_distance(running);
        std::cout << "--- initial filter took " << timer.Stop() / 1000 << " seconds " << std::endl << std::endl;
    }

    // save peaks
    save(running, ds.first_checkpoint_path());

    this->run_from_last_initial(true);
}

void pipeline::run_from_last_initial(bool is_full_run) {
    std::vector<haplotype *> running = this->recover(ds.first_checkpoint_path()); 

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

        ////////////////////////////////////REMOVE when working with real data
        //std::vector<haplotype*> current;
        //current.reserve(full.size());
        //std::transform(full.begin(), full.end(), std::back_inserter(current),
        //           [](const std::pair<haplotype*, double>& pair) {
        //               return pair.first;
        //           });
        //a.print_mutation_distance(current);
        /////////////////////////////////////////

        std::cout << "--- post filter took " << timer.Stop() / 1000 << " seconds " << std::endl;

        a.dump_haplotype_proportion(full);
        a.dump_read2node_mapping(full);
    }
}