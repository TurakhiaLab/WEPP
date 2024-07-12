#include <vector>
#include <unordered_map>

#include "variant_finder.hpp"

void
variant_finder::find(arena& arena, const std::vector<std::pair<haplotype*, double>>& selected) 
{
    std::vector<size_t> occ(arena.genome_size());
    std::vector<std::unordered_map<std::string, size_t>> kmer(arena.genome_size());

    const std::string reference = arena.reference();

    for (const raw_read& read: arena.reads()) {
        for (size_t i = read.start; i <= read.end - this->k; ++i) {
            occ[i] += read.degree;

            std::string build;
            for (size_t j = i; j < i + (size_t) this->k; ++j) {
                MAT::Mutation search;
                search.position = j;
                auto it = std::lower_bound(read.mutations.begin(), read.mutations.end(), search);
                if (it != read.mutations.end() && (size_t) it->position == j) {
                    build += MAT::get_nuc(it->mut_nuc);
                }
                else {
                    build += reference[j];
                }
            }

            kmer[i][build] += read.degree;
        }
    }

    size_t gs = arena.genome_size();
    size_t num_bad = 0;
    for (size_t i = 0; i < gs; ++i) {
        if (occ[i] == 0) continue;

        for (const auto& [k, v] : kmer[i]) {
            if ((double) v / occ[i] > this->freq_cutoff) {
                // ensure it has a match
                bool found = false;
                for (const auto& [hap, _] : selected) {
                    bool matches = true;
                    for (size_t j = i; j < i + (size_t) this->k; ++j) {
                        MAT::Mutation search;
                        search.position = j;
                        auto it = std::lower_bound(hap->stack_muts.begin(), hap->stack_muts.end(), search);
                        char c;
                        if (it != hap->stack_muts.end() && (size_t) it->position == j)
                        {
                            c = MAT::get_nuc(it->mut_nuc);
                        }
                        else
                        {
                            c = reference[j];
                        }

                        if (c != k[j - i]) {
                            matches = false;
                            break;
                        }
                    }
                    if (matches) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    num_bad++;
                }
            }
        }
    }

    std::cout << "NUM BAD" << num_bad << std::endl;
}