#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>

class SmithWaterman {
private:
    int match_score;
    int mismatch_penalty;
    int gap_penalty;

public:
    SmithWaterman(int match_score, int mismatch_penalty, int gap_penalty) {
        this->match_score = match_score;
        this->mismatch_penalty = mismatch_penalty;
        this->gap_penalty = gap_penalty;
    }

    std::vector<std::pair<std::string, int> > align(std::string reference, std::vector<std::string> reads) {
        std::vector<std::pair<std::string, int> > alignment_scores;
        for (size_t i = 0; i < reads.size(); i++) {
            alignment_scores.push_back(smith_waterman(reference, reads[i]));
        }
        return alignment_scores;
    }

private:
    std::pair<std::string, int> smith_waterman(const std::string& reference, const std::string& read) {
        // Initialize the scoring matrix
        int rows = read.length() + 1;
        int cols = reference.length() + 1;
        std::vector<std::vector<int> > scores(rows, std::vector<int>(cols, 0));

        int max_score = 0;
        std::pair<int, int> max_pos;

        // Fill in the scoring matrix
        for (int i = 1; i < rows; i++) {
            for (int j = 1; j < cols; j++) {
                int score_diagonal = scores[i - 1][j - 1] + (read[i - 1] == reference[j - 1] ? match_score : mismatch_penalty);
                int score_up = scores[i - 1][j] + gap_penalty;
                int score_left = scores[i][j - 1] + gap_penalty;

                scores[i][j] = std::max(0, std::max(score_diagonal, std::max(score_up, score_left)));

                if (scores[i][j] > max_score) {
                    max_score = scores[i][j];
                    max_pos = std::make_pair(i, j);
                }
            }
        }

        // Traceback to find the alignment and position
        std::string alignment;
        int i = max_pos.first;
        int j = max_pos.second;
        int reference_position = j - 1;  // Subtract 1 to get the 0-based index

        while (i > 0 && j > 0 && scores[i][j] != 0) {
            if (scores[i][j] == scores[i - 1][j - 1] + (read[i - 1] == reference[j - 1] ? match_score : mismatch_penalty)) {
                alignment = read[i - 1] + alignment;
                i -= 1;
                j -= 1;
            } else if (scores[i][j] == scores[i - 1][j] + gap_penalty) {
                alignment = read[i - 1] + alignment;
                i -= 1;
            } else {
                alignment = reference[j - 1] + alignment;
                j -= 1;
            }
        }

        return std::make_pair(alignment, reference_position);
    }
};

int main(int argc, char** argv) {
    // Check that the correct number of arguments are provided
    if (argc != 4) {
        std::cout << "Usage: ./local_alignment <reference_fasta> <fasta with reads> <output aligned fasta>" << std::endl;
        return 1;
    }

    // Read in the reference fasta file and store the sequence in a string called reference_sequence
    std::string reference_fasta = argv[1];
    std::ifstream ref_fasta(reference_fasta);

    std::string output_file = argv[3];
    std::ofstream output(output_file, std::ios_base::app);
    
    std::vector<std::string> ref_lines;
    std::string line;
    while (std::getline(ref_fasta, line)) {
        ref_lines.push_back(line);
    }
    ref_fasta.close();
    std::string reference_sequence = "";
    for (size_t i = 1; i < ref_lines.size(); i++) {
        reference_sequence += ref_lines[i];
    }

    // Read in the fasta file with reads and store the sequences in a vector called reads
    // that fasta file should have the following format:
    // >read1
    // GTCG.... etc.
    // >read2
    // GCGTAC.... etc.
    std::string reads_fasta = argv[2];
    std::ifstream reads_file(reads_fasta);
    std::vector<std::string> reads;
    std::vector<std::string> read_names;
    std::string read = "";
    while (std::getline(reads_file, line)) {
        if (line[0] == '>') {
            if (read != "") {
                reads.push_back(read);
                read = "";
                
            }
            read_names.push_back(line);
        } else {
            read += line;
        }
    }
    reads.push_back(read);

    SmithWaterman sw(2, -1, -2);
    std::vector<std::pair<std::string, int> > alignment_results = sw.align(reference_sequence, reads);

    for (size_t i = 0; i < reads.size(); i++) {
        std::string alignment = alignment_results[i].first;
        int ref_position = alignment_results[i].second;

        // FOR DEBUGGING
        // std::cout << "Read " << i + 1 << ": " << reads[i] << std::endl;
        // std::cout << "Alignment Score: " << alignment.length() << std::endl;
        // std::cout << "Alignment: " << alignment << std::endl;
        // std::cout << "Starting Position: " << ref_position - alignment.length() + 2 << std::endl;
        // std::cout << "Ending Position: " << ref_position + 1 << std::endl;
        // std::cout << std::endl;

        // Write the aligned reads to the output file
        output << read_names[i] << "_" << ref_position - alignment.length() + 2 << "_" << ref_position + 1 << std::endl;
        output << alignment << std::endl;
        
    }

    output.close();

    return 0;
}
