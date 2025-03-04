#include "analyze_peaks.hpp"

void read_haplotype_proportion(std::vector<std::tuple<std::string, double, std::vector<MAT::Mutation>>> &abundance, std::string haplotype_proportion_path)
{
    std::ifstream csv(haplotype_proportion_path);
    std::string line;
    while (std::getline(csv, line))
    {
        std::stringstream ss(line);
        std::string id;
        std::string proportion_str;

        // Assuming the format is: haplotype_id,proportion
        std::getline(ss, id, ',');
        std::getline(ss, proportion_str, ',');

        double proportion = std::stod(proportion_str);

        abundance.emplace_back(std::make_tuple(id, proportion, std::vector<MAT::Mutation>()));
    }
}


void update_cluster(const std::vector<std::tuple<std::string, double, std::vector<MAT::Mutation>>>& curr_hap_abun_muts, std::vector<bool>& check_peaks, std::vector<int>& cluster_peaks, double& cluster_abundance)
{
    bool updated = false;
    for (size_t i = 0; i < curr_hap_abun_muts.size(); i++)
    {
        if (check_peaks[i]) {
            auto curr_mutations = std::get<2>(curr_hap_abun_muts[i]);
            for (auto idx: cluster_peaks) {
                auto cmp_mutations = std::get<2>(curr_hap_abun_muts[idx]);
                if (mutation_distance(curr_mutations, cmp_mutations) <= CLUSTER_DIST)
                {
                    updated = true;
                    check_peaks[i] = false;
                    cluster_peaks.emplace_back(i);
                    cluster_abundance += std::get<1>(curr_hap_abun_muts[i]);
                    break;
                }
            }
        }
    }
    if (updated)
        update_cluster(curr_hap_abun_muts, check_peaks, cluster_peaks, cluster_abundance);
}