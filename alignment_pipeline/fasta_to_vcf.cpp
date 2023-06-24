#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <sstream>

std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(line);
    
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

int main(int argc, char* argv[]) {
    // Check that the correct number of arguments are provided
    if (argc != 4) {
        std::cout << "Usage: fasta_to_vcf <input_fasta_file> <reference_fasta> <output_vcf>" << std::endl;
        return 1;
    }

    std::string reference = argv[2];
    std::string input_file = argv[1];
    std::string output_file = argv[3];

    std::ifstream ref_file(reference);
    std::vector<std::string> ref_lines;

    if (!ref_file.is_open()) {
        std::cout << "Error opening reference file" << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(ref_file, line)) {
        ref_lines.push_back(line);
    }
    ref_file.close();

    int num_samples = 0;
    std::vector<std::string> aligned_reads;
    std::vector<std::string> sample_names;

    std::ifstream sample_names_file("./intermediate_files/sample_names.txt");
    if (!sample_names_file.is_open()) {
        std::cout << "Error opening sample names file" << std::endl;
        return 1;
    }

    std::getline(sample_names_file, line);
    sample_names = split(line, ',');
    sample_names_file.close();

    std::ifstream read_file(input_file);
    std::vector<std::string> read_lines;

    if (!read_file.is_open()) {
        std::cout << "Error opening input file" << std::endl;
        return 1;
    }

    while (std::getline(read_file, line)) {
        read_lines.push_back(line);
        if (line[0] == '>') {
            num_samples++;
        } else {
            aligned_reads.push_back(line);
        }
    }
    read_file.close();

    std::string ref_seq;
    std::string chromosome_name;

    for (const auto& line : ref_lines) {
        if (line[0] == '>') {
            chromosome_name = line.substr(1);
        } else {
            ref_seq += line;
        }
    }

    std::ofstream vcf_file(output_file);
    if (!vcf_file.is_open()) {
        std::cout << "Error opening output file" << std::endl;
        return 1;
    }

    vcf_file << "##fileformat=VCFv4.2\n";
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

    for (int i = 0; i < num_samples; i++) {
        vcf_file << sample_names[i] << "\t";
    }

    vcf_file << "\n";

    for (int i = 0; i < ref_seq.length(); i++) {
        bool mismatch = false;
        std::vector<int> mismatch_indices;

        for (int k = 0; k < num_samples; k++) {
            if (ref_seq[i] != aligned_reads[k][i]) {
                mismatch = true;
                mismatch_indices.push_back(k);
            }
        }

        if (mismatch) {
            std::set<char> mutations;
            for (const auto& index : mismatch_indices) {
                mutations.insert(aligned_reads[index][i]);
            }

            for (const auto& mutation : mutations) {
                vcf_file << chromosome_name << "\t";
                vcf_file << i + 1 << "\t";
                char ref = ref_seq[i];
                char alt = mutation;
                std::string id = std::string(1, ref) + std::to_string(i + 1) + std::string(1, alt);
                vcf_file << id << "\t";
                vcf_file << ref << "\t";
                vcf_file << alt << "\t";
                vcf_file << ".\t";
                vcf_file << ".\t";
                vcf_file << ".\t";
                vcf_file << ".\t";

                for (int j = 0; j < num_samples; j++) {
                    if (aligned_reads[j][i] == mutation) {
                        vcf_file << "1\t";
                    } else {
                        vcf_file << "0\t";
                    }
                }

                vcf_file << "\n";
            }
        }
    }

    vcf_file.close();

    return 0;
}
