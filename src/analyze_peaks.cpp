#include "analyze_peaks.hpp"

static constexpr int NEIGHBOR_DIST = 2;

void analyze_peaks(const dataset& d) {
    // Read files
    // std::vector<std::pair<std::string, double>> curr_hap_abundance, cmp_hap_abundance; 
    // read_haplotype_proportion(curr_hap_abundance, d.haplotype_proportion_path()); 
    // read_haplotype_proportion(cmp_hap_abundance, d.comparison_haplotype_proportion_path());
    // auto T_curr = d.mat(false); 
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
    // std::unordered_map<std::string, std::vector<MAT::Mutation>> curr_hap_mutations, cmp_hap_mutations;
    // std::vector<std::vector<std::string>> curr_hap_neighbors;
    // for (auto h_a: curr_hap_abundance) 
    // {
    //     auto sample_mutations = get_mutations(T_curr, h_a.first);
    //     auto mut_itr = sample_mutations.begin();
    //     while (mut_itr != sample_mutations.end())
    //     {
    //         if (site_read_map_curr.find(mut_itr->position) == site_read_map_curr.end())
    //             mut_itr = sample_mutations.erase(mut_itr);
    //         else
    //             mut_itr++;
    //     }
    //     curr_hap_mutations.insert({h_a.first, sample_mutations});

    //     // Update curr_hap_neighbors
    //     bool found = false;
    //     for (auto& h_n: curr_hap_neighbors) {
    //         for (auto hap: h_n) {
    //             auto cmp_mutations = curr_hap_mutations[hap];
    //             if (mutation_distance(sample_mutations, cmp_mutations) <= NEIGHBOR_DIST)
    //             {
    //                 found = true;
    //                 break;
    //             }
    //         }
    //         if (found) 
    //         {
    //             h_n.emplace_back(h_a.first);
    //             break;
    //         }
    //     }
    //     if (!found)
    //         curr_hap_neighbors.emplace_back(std::vector<std::string>{h_a.first}); 
    // }

    // for (auto h_a: cmp_hap_abundance) 
    // {
    //     auto sample_mutations = get_mutations(T_cmp, h_a.first);
    //     auto mut_itr = sample_mutations.begin();
    //     while (mut_itr != sample_mutations.end())
    //     {
    //         if (site_read_map_cmp.find(mut_itr->position) == site_read_map_cmp.end())
    //             mut_itr = sample_mutations.erase(mut_itr);
    //         else
    //             mut_itr++;
    //     }
    //     cmp_hap_mutations.insert({h_a.first, sample_mutations});
    // }

    // // Find neighboring cmp_peaks wrt curr_peaks
    // std::ofstream txt(d.haplotype_growth_path());
    // for (auto c_h_n: curr_hap_neighbors) 
    // {
    //     std::set<std::string> curr_lineages, cmp_lineages;
    //     std::vector<std::string> cmp_hap_considered;
    //     double curr_abundance = 0.0, cmp_abundance = 0.0;
    //     for (auto curr_hap: c_h_n)
    //     {
    //         auto curr_mutations = curr_hap_mutations[curr_hap];
    //         auto curr_abun_itr = std::find_if(curr_hap_abundance.begin(), curr_hap_abundance.end(),
    //                        [curr_hap](const std::pair<std::string, double>& element) {
    //                            return element.first == curr_hap;
    //                        });
    //         curr_abundance += curr_abun_itr->second; 

    //         // Get curr_lineages 
    //         for (auto anc : T_curr.rsearch(curr_hap, true))
    //         {
    //             const auto &clade = anc->clade_annotations[1];
    //             if (clade != "")
    //             {
    //                 curr_lineages.insert(clade);
    //                 break;
    //             }
    //         }

    //         // Check neighbors in cmp_hap_abundance
    //         for (auto cmp_h_a: cmp_hap_abundance) 
    //         {
    //             if (std::find(cmp_hap_considered.begin(), cmp_hap_considered.end(), cmp_h_a.first) != cmp_hap_considered.end())
    //                 continue;
    //             auto cmp_mutations = cmp_hap_mutations[cmp_h_a.first];
    //             if (mutation_distance(curr_mutations, cmp_mutations) <= NEIGHBOR_DIST) {
    //                 cmp_hap_considered.emplace_back(cmp_h_a.first);
    //                 cmp_abundance += cmp_h_a.second;

    //                 // Get cmp_lineages 
    //                 for (auto anc : T_cmp.rsearch(cmp_h_a.first, true))
    //                 {
    //                     const auto &clade = anc->clade_annotations[1];
    //                     if (clade != "")
    //                     {
    //                         cmp_lineages.insert(clade);
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     if (cmp_abundance)
    //     {
    //         // Write lineage and abundances
    //         txt << "(" + std::to_string(c_h_n.size()) + " Peaks - ";
    //         for (auto it = curr_lineages.begin(); it != curr_lineages.end(); it++) 
    //         {
    //             if (it != curr_lineages.begin())
    //                 txt << ",";
    //             txt << *it;
    //         }
    //         txt << ")\t" + std::to_string(curr_abundance) + "\t->\t(" + std::to_string(cmp_hap_considered.size()) + " Peaks - ";
    //         for (auto it = cmp_lineages.begin(); it != cmp_lineages.end(); it++) 
    //         {
    //             if (it != cmp_lineages.begin())
    //                 txt << ",";
    //             txt << *it;
    //         }
    //         txt << ")\t" + std::to_string(cmp_abundance) + "\n";

    //         // Write peak names
    //         for (size_t i = 0; i < c_h_n.size(); i++)
    //         {
    //             if (i)
    //                 txt << ",";
    //             txt << c_h_n[i];
    //         }
    //         txt << "\t->\t";
    //         for (size_t i = 0; i < cmp_hap_considered.size(); i++)
    //         {
    //             if (i)
    //                 txt << ",";
    //             txt << cmp_hap_considered[i];
    //         }
    //         txt << "\n";
    //     }
    // }
}

void read_haplotype_proportion(std::vector<std::pair<std::string, double>> &abundance, std::string haplotype_proportion_path)
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

        abundance.emplace_back(id, proportion);
    }
}