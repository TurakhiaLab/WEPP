#include "post_filter.hpp"

#include <numeric>
#include <random>

void 
freyja_post_filter::dump_barcode(arena& a, const std::vector<haplotype*>& haplotypes) 
{
    std::set<std::string> mutations;
    std::vector<std::set<std::string>> node_muts;

    for (haplotype *n : haplotypes)
    {
        std::set<std::string> my_muts;
        for (MAT::Mutation mut : n->stack_muts)
        {
            std::string build = MAT::get_nuc(mut.ref_nuc) + std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc);
            mutations.insert(build);
            my_muts.insert(std::move(build));
        }
        node_muts.push_back(std::move(my_muts));
    }

    for (const raw_read& read: a.reads()) {
        for (MAT::Mutation mut : read.mutations)
        {
            uint8_t const N = 0b1111;
            if (mut.mut_nuc != N) {
                std::string build = MAT::get_nuc(mut.ref_nuc) + std::to_string(mut.position) + MAT::get_nuc(mut.mut_nuc);
                mutations.insert(build);
            }
        }
    }

    std::vector<std::string> mutation_vec(mutations.begin(), mutations.end());

    std::ofstream outfile(a.owned_dataset().intermediate_directory() + a.owned_dataset().file_prefix() + "_barcodes.csv");
    for (const std::string &mut : mutations)
    {
        outfile << "," << mut;
    }
    outfile << std::endl;
    for (size_t i = 0; i < haplotypes.size(); ++i)
    {
        haplotype *n = haplotypes[i];
        outfile << 'N' << (n - &a.haplotypes()[0]);
        for (const std::string &mut : mutations)
        {
            bool contains = node_muts[i].find(mut) != node_muts[i].end();
            outfile << "," << (contains ? '1' : '0');
        }
        outfile << std::endl;
    }
}

std::vector<std::pair<haplotype*, double>>
freyja_post_filter::filter(arena& arena, std::vector<haplotype*> input)
{
    fprintf(stderr, "%ld peaks selected for Freyja!\n\n", input.size());
    dump_barcode(arena, input);

    auto const& dataset = arena.owned_dataset();

    std::string command = "bash -c \""
                "cd ./src/Freyja/ && "
                "freyja demix "
                    "'../../" + dataset.intermediate_directory() + dataset.file_prefix() + "_corrected_variants.tsv' "
                    "'../../" + dataset.intermediate_directory() + dataset.file_prefix() + "_depth.tsv' "
                    "--barcodes '../../" + dataset.intermediate_directory() + dataset.file_prefix() + "_barcodes.csv' "
                    "--output '../../" + dataset.intermediate_directory() + "freyja_output_latest.txt' --eps 0.005"
                "\"";
    if (std::system(command.c_str()) != 0)
    {
        std::cerr << "Failed to run freyja" << std::endl;
        return {};
    }

    std::vector<std::pair<haplotype *, double>> freyja_nodes;

    std::ifstream fin(dataset.intermediate_directory() + "freyja_output_latest.txt");
    std::string tmp;
    std::getline(fin, tmp);
    std::getline(fin, tmp);
    fin >> tmp;
    // all of the selected ids
    std::getline(fin, tmp);
    std::stringstream ss{tmp};
    std::string index;
    while (ss >> index)
    {
        // skip past the 'N'
        int ind = atoi(index.c_str() + 1);
        freyja_nodes.emplace_back(&arena.haplotypes()[ind], 0);
    }

    fin >> tmp;
    std::getline(fin, tmp);
    ss = std::stringstream{tmp};
    double sum = 0;
    for (size_t i = 0; i < freyja_nodes.size(); ++i) {
        double abundance; ss >> abundance;
        freyja_nodes[i].second = abundance;
        sum += abundance;
    }

    for (size_t i = 0; i < freyja_nodes.size(); ++i) {
        freyja_nodes[i].second /= sum;
    }

    return freyja_nodes;
}

std::vector<std::pair<haplotype*, double>>
em_post_filter::filter(arena& arena, std::vector<haplotype*> input)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    const std::vector<raw_read>& reads = arena.reads();

    std::vector<std::vector<double>> q(input.size(), std::vector<double>(reads.size()));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                      [&](tbb::blocked_range<size_t> range)
                      {
                          for (size_t i = range.begin(); i < range.end(); ++i)
                          {
                              for (size_t j = 0; j < reads.size(); ++j)
                              {
                                  double dist = input[i]->mutation_distance(reads[j]);
                                  q[i][j] = std::pow(this->alpha, dist) * std::pow(1 - this->alpha, reads[j].end - reads[j].start + 1 - dist);
                              }
                          }
                      });

    double prev = 0, curr = 0;
    // Initialize to normal distribtuion
    std::vector<double> p(input.size());
    for (auto& val : p) {
        val = dis(gen);
    }
    // Re-scale to 1
    double sum = std::accumulate(p.begin(), p.end(), 0.0);
    for (auto& val : p) {
        val /= sum;
    }

    int max_iteration = max_it; 
    do
    {
        prev = curr;

        std::vector<double> denom(reads.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t j = range.begin(); j < range.end(); ++j)
                              {
                                  double sum = 0;
                                  for (size_t l = 0; l < input.size(); ++l)
                                  {
                                      sum += p[l] * q[l][j];
                                  }
                                  denom[j] = sum;
                              }
                          });

        tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t i = range.begin(); i < range.end(); ++i)
                              {
                                  double sum = 0;
                                  size_t count = 0;
                                  for (size_t j = 0; j < reads.size(); ++j)
                                  {
                                      sum += reads[j].degree * p[i] * q[i][j] / denom[j];
                                      count += reads[j].degree;
                                  }
                                  // don't even need multiple p
                                  p[i] = sum / count;
                              }
                          });

        std::vector<double> logl(reads.size());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, reads.size()),
                          [&](tbb::blocked_range<size_t> range)
                          {
                              for (size_t j = range.begin(); j < range.end(); ++j)
                              {
                                  double sum = 0;
                                  for (size_t l = 0; l < input.size(); ++l)
                                  {
                                      sum += p[l] * q[l][j];
                                  }
                                  logl[j] = reads[j].degree * std::log(sum);
                              }
                          });

        curr = std::accumulate(logl.begin(), logl.end(), 0.0);
        if (!prev)
            prev = 1e-12;
    } while ((abs(curr - prev) / abs(prev)) >= epsilon && --max_iteration);

    // Select haplotypes with abundance > min_proportion
    std::vector<std::pair<haplotype *, double>> em_nodes;
    for (size_t i = 0; i < p.size(); ++i)
    {
        em_nodes.emplace_back(std::make_pair(input[i], p[i]));
    }

    std::sort(em_nodes.begin(), em_nodes.end(),
    [](const std::pair<haplotype*, double>& a, const std::pair<haplotype*, double>& b) {
        return a.second > b.second;
    });
    auto it = em_nodes.begin();
    double total = 0;
    while ((it != em_nodes.end()) && (it->second / (it->second + total) >= min_proportion))
    {
        total += it->second;
        ++it;
    }
    em_nodes.erase(it, em_nodes.end());

    for (auto& node_score: em_nodes)
    {
        node_score.second /= total;
    }

    return em_nodes;
}