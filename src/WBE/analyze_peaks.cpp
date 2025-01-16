#include "analyze_peaks.hpp"

void analyze_peaks(const dataset& d) {
    // // Read files
    // std::vector<std::tuple<std::string, double, std::vector<MAT::Mutation>>> curr_hap_abun_muts, cmp_hap_abun_muts; 
    // read_haplotype_proportion(curr_hap_abun_muts, d.haplotype_proportion_path()); 
    // read_haplotype_proportion(cmp_hap_abun_muts, d.comparison_haplotype_proportion_path());
    // auto T_curr = d.mat(); 
    // auto T_cmp = d.cmp_mat();
    
    // std::unordered_map<std::string, std::vector<std::string>> reverse;
    // auto curr_reads = load_reads_from_proto(d.reference(), d.pb_path(), reverse);
    // auto cmp_reads = load_reads_from_proto(d.reference(), d.comparison_pb_path(), reverse);

    // // Create site read maps
    // std::unordered_set<int> site_read_map_curr, site_read_map_cmp;
    // for (size_t i = 0; i < curr_reads.size(); i++)
    // {
    //     const auto &rp = curr_reads[i];
    //     std::unordered_set<int> ambiguous_sites;
    //     for (auto mut: rp.mutations) {
    //         if (mut.mut_nuc == 0b1111)
    //             ambiguous_sites.insert(mut.position);
    //     }
    //     for (int j = rp.start; j <= rp.end; j++) {
    //         if (ambiguous_sites.find(j) == ambiguous_sites.end())
    //             site_read_map_curr.insert(j);
    //     }
    //     ambiguous_sites.clear();
    // }

    // for (size_t i = 0; i < cmp_reads.size(); i++)
    // {
    //     const auto &rp = cmp_reads[i];
    //     std::unordered_set<int> ambiguous_sites;
    //     for (auto mut: rp.mutations) {
    //         if (mut.mut_nuc == 0b1111)
    //             ambiguous_sites.insert(mut.position);
    //     }
    //     for (int j = rp.start; j <= rp.end; j++) {
    //         if (ambiguous_sites.find(j) == ambiguous_sites.end())
    //             site_read_map_cmp.insert(j);
    //     }
    //     ambiguous_sites.clear();
    // }
    // fprintf(stderr, "Sites covered by ref_reads: %ld, cmp_reads: %ld\n", site_read_map_curr.size(), site_read_map_cmp.size());

    // // Store mutation vectors of peaks
    // for (auto& h_a_m: curr_hap_abun_muts) 
    // {
    //     auto sample_mutations = get_mutations(T_curr, std::get<0>(h_a_m));
    //     auto mut_itr = sample_mutations.begin();
    //     while (mut_itr != sample_mutations.end())
    //     {
    //         if (site_read_map_curr.find(mut_itr->position) == site_read_map_curr.end())
    //             mut_itr = sample_mutations.erase(mut_itr);
    //         else
    //             mut_itr++;
    //     }
    //     std::get<2>(h_a_m) = sample_mutations;
    // }

    // for (auto& h_a_m: cmp_hap_abun_muts) 
    // {
    //     auto sample_mutations = get_mutations(T_cmp, std::get<0>(h_a_m));
    //     auto mut_itr = sample_mutations.begin();
    //     while (mut_itr != sample_mutations.end())
    //     {
    //         if (site_read_map_cmp.find(mut_itr->position) == site_read_map_cmp.end())
    //             mut_itr = sample_mutations.erase(mut_itr);
    //         else
    //             mut_itr++;
    //     }
    //     std::get<2>(h_a_m) = sample_mutations;
    // }

    // // Create Clusters from curr_hap_mutations
    // std::vector<std::pair<std::vector<int>, double>> curr_cluster_proportions;
    // std::vector<bool> check_peaks(curr_hap_abun_muts.size(), true);
    // for (int i = 0; i < (int)curr_hap_abun_muts.size(); i++)
    // {
    //     if (check_peaks[i]) {
    //         check_peaks[i] = false;
    //         curr_cluster_proportions.emplace_back(std::make_pair(std::vector<int>{i}, std::get<1>(curr_hap_abun_muts[i])));
    //         update_cluster(curr_hap_abun_muts, check_peaks, curr_cluster_proportions.back().first, curr_cluster_proportions.back().second);
    //     }   
    // }
    // check_peaks.clear();

    // std::sort(curr_cluster_proportions.begin(), curr_cluster_proportions.end(),
    //     [](const std::pair<std::vector<int>, double>& a, const std::pair<std::vector<int>, double>& b) {
    //         return a.second > b.second; 
    //     });


    // // Find neighboring cmp_peaks wrt curr_cluster
    // std::vector<std::tuple<int, std::vector<int>, double>> cmp_cluster_proportions;
    // check_peaks.resize(cmp_hap_abun_muts.size(), true); 
    // for (int i = 0; i < (int)cmp_hap_abun_muts.size(); i++) 
    // {
    //     auto curr_mutations = std::get<2>(cmp_hap_abun_muts[i]);
    //     for (size_t j = 0; j < curr_cluster_proportions.size(); j++)
    //     {
    //         for (const auto& idx: curr_cluster_proportions[j].first)
    //         {
    //             auto ref_mutations = std::get<2>(curr_hap_abun_muts[idx]);
    //             if (mutation_distance(curr_mutations, ref_mutations) <= CLUSTER_DIST)
    //             {
    //                 check_peaks[i] = false;
    //                 bool found = false;
    //                 for (auto& i_c_p: cmp_cluster_proportions) {
    //                     if (std::get<0>(i_c_p) == (int)j) {
    //                         std::get<1>(i_c_p).emplace_back(i);
    //                         std::get<2>(i_c_p) += std::get<1>(cmp_hap_abun_muts[i]);
    //                         found = true;
    //                         break;
    //                     }
    //                 }
    //                 if (!found)
    //                     cmp_cluster_proportions.emplace_back(std::make_tuple(j, std::vector<int>{i}, std::get<1>(cmp_hap_abun_muts[i])));
    //                 break;
    //             }
    //         }
    //         if (!check_peaks[i])
    //             break;
    //     }
    // } 


    // // Update local clusters
    // for (int i = 0; i < (int)cmp_hap_abun_muts.size(); i++) 
    // {
    //     if (check_peaks[i]) {
    //         check_peaks[i] = false;
    //         cmp_cluster_proportions.emplace_back(std::make_tuple((1000 + cmp_cluster_proportions.size()), std::vector<int>{i}, std::get<1>(cmp_hap_abun_muts[i])));
    //         update_cluster(cmp_hap_abun_muts, check_peaks, std::get<1>(cmp_cluster_proportions.back()), std::get<2>(cmp_cluster_proportions.back()));
    //     }
    // }


    // // Dump the clusters found
    // std::ofstream txt(d.haplotype_growth_path());
    // std::vector<double> proportion_difference;
    // for (size_t i = 0; i < curr_cluster_proportions.size(); i ++) {
    //     double prop_diff = 0.0;
    //     std::string curr_cluster_nodes = "";
    //     for (auto idx: curr_cluster_proportions[i].first) {
    //         if (curr_cluster_nodes.size())
    //             curr_cluster_nodes += "__";
    //         curr_cluster_nodes += std::get<0>(curr_hap_abun_muts[idx]);
    //     }
    //     prop_diff = curr_cluster_proportions[i].second;
    //     txt << std::to_string(curr_cluster_proportions[i].second) + "," + curr_cluster_nodes + ",";
                    
    //     bool found = false;
    //     for (auto& i_c_p: cmp_cluster_proportions) {
    //         if (std::get<0>(i_c_p) == (int)i) {
    //             std::string cmp_cluster_nodes = "";
    //             for (auto idx: std::get<1>(i_c_p)) {
    //                 if (cmp_cluster_nodes.size())
    //                     cmp_cluster_nodes += "__";
    //                 cmp_cluster_nodes += std::get<0>(cmp_hap_abun_muts[idx]);
    //             }
    //             prop_diff -= std::get<2>(i_c_p);
    //             txt << std::to_string(std::get<2>(i_c_p)) + "," + cmp_cluster_nodes + "\n";
    //             found = true;
    //             break;
    //         }
    //     }
    //     if (!found)
    //         txt << "0,\n";

    //     proportion_difference.emplace_back(prop_diff);
    // }

    // for (auto& i_c_p: cmp_cluster_proportions) {
    //     if (std::get<0>(i_c_p) >= 1000) {
    //         std::string cmp_cluster_nodes = "";
    //         for (auto idx: std::get<1>(i_c_p)) {
    //             if (cmp_cluster_nodes.size())
    //                 cmp_cluster_nodes += "__";
    //             cmp_cluster_nodes += std::get<0>(cmp_hap_abun_muts[idx]);
    //         }
    //         proportion_difference.emplace_back(std::get<2>(i_c_p));
    //         txt << "0,," + std::to_string(std::get<2>(i_c_p)) + "," + cmp_cluster_nodes + "\n";
    //     }
    // }

    // double mse_both = 0.0, mse_ref = 0.0;
    // for (size_t i = 0; i < proportion_difference.size(); i++) {
    //     double curr_sd = pow(proportion_difference[i], 2);
    //     if (i < curr_cluster_proportions.size())
    //         mse_ref += curr_sd;
    //     mse_both += curr_sd;
    // }
    // mse_ref /= curr_cluster_proportions.size();
    // mse_both /= proportion_difference.size();
    // fprintf(stderr, "MSE (Ref clusters): %f\n", mse_ref);
    // fprintf(stderr, "MSE (Both clusters): %f\n", mse_both);
}


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