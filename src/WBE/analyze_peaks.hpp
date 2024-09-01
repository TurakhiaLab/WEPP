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


void analyze_peaks(const dataset& d);
void read_haplotype_proportion(std::vector<std::pair<std::string, double>> &abundance, std::string haplotype_proportion_path);