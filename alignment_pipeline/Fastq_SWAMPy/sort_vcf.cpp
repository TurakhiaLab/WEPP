#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

int main(const int argc, const char* argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_vcf> <output_vcf>" << std::endl;
        return 1;
    }

    std::ifstream inputFile(argv[1]);  
    std::ofstream outputFile(argv[2]);

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files." << std::endl;
        return 1;
    }

    std::vector<std::pair<int, std::string>> dataLines;
    std::string line;

    // Read lines from the VCF file
    while (getline(inputFile, line)) {
        // If it's a meta-info or header line, write directly to the output
        if (line.substr(0, 1) == "#") {
            outputFile << line << std::endl;
        } else {
            std::stringstream ss(line);
            std::string chrom, pos;
            ss >> chrom >> pos;  // Extract CHROM and POS columns

            int position = std::stoi(pos);
            dataLines.push_back({position, line});
        }
    }

    // Sort data lines based on the POS column
    std::sort(dataLines.begin(), dataLines.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    // Write sorted lines to the output
    for (const auto& pair : dataLines) {
        outputFile << pair.second << std::endl;
    }

    inputFile.close();
    outputFile.close();
    return 0;
}
