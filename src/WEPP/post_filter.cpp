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

    double af_thresh = 0.0;
    if (arena.min_af() > 0.01)
        af_thresh = arena.min_af();
    
    std::string command = "bash -c \""
                "cd ./src/Freyja/ && "
                "freyja demix "
                    "'../../" + dataset.intermediate_directory() + dataset.file_prefix() + "_corrected_variants.tsv' "
                    "'../../" + dataset.intermediate_directory() + dataset.file_prefix() + "_depth.tsv' "
                    "--barcodes '../../" + dataset.intermediate_directory() + dataset.file_prefix() + "_barcodes.csv' "
                    "--output '../../" + dataset.intermediate_directory() + "freyja_output_latest.txt' --eps " + std::to_string(arena.min_prop()) + " --af " + std::to_string(af_thresh) +
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