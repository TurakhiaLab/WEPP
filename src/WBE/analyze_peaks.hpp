#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>

#include "dataset.hpp"
#include "haplotype.hpp"
#include "sam2pb.hpp"
#include "config.hpp"


void analyze_peaks(const dataset& d);
void read_haplotype_proportion(std::vector<std::tuple<std::string, double, std::vector<MAT::Mutation>>> &abundance, std::string haplotype_proportion_path);
void update_cluster(const std::vector<std::tuple<std::string, double, std::vector<MAT::Mutation>>>& curr_hap_abun_muts, std::vector<bool>& check_peaks, std::vector<int>& cluster_peaks, double& cluster_abundance);