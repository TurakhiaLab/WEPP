#include <memory> 

#include "pipeline.hpp"

void detect_peaks(const dataset& d) {
    auto main = std::make_unique<wepp_filter>();
    auto post = std::make_unique<em_post_filter>();

    pipeline p{d, std::move(main), std::move(post)};
    p.run();
}

void pipeline::run() {
    tbb::task_scheduler_init init(this->ds.num_threads());

    std::vector<haplotype *> running = a.haplotype_pointers();

    std::cout << "--- [running main filter] ---" << std::endl;
    std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;
    running = main->filter(a);
    a.print_mutation_distance(running);

    std::cout << "--- [running post filter] ---" << std::endl;
    std::cout << "--- in: " << running.size() << " haplotypes" << std::endl;
    running = post->filter(a, running);
    a.print_mutation_distance(running);
}
