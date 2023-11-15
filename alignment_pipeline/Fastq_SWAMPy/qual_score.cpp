#include <iostream>
#include <fstream>
#include <vector>

int main() {
    std::ifstream fastq_file("40.fastq");
    if (!fastq_file.is_open()) {
        std::cerr << "Error opening file\n";
        return 1;
    }

    std::string line;
    int line_count = 0;
    while (getline(fastq_file, line)) {
        ++line_count;
        if (line_count % 4 == 0) {  
            std::vector<int> quality_scores;
            for (char c : line) {
                int score = static_cast<int>(c) - 33; // depending on the encoding, this may need to be changed
                quality_scores.push_back(score);
            }
            
            for (int score : quality_scores) {
                std::cout << score << " ";
            }
            std::cout << std::endl;
        }
    }

    fastq_file.close();
    return 0;
}
