#include "analyze_peaks.hpp"

static constexpr int NEIGHBOR_DIST = 2;

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