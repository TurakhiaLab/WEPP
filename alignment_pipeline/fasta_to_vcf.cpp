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
    std::vector<int> read_start_positions;
    std::vector<int> read_end_positions;

    std::ifstream read_file(input_file);
    std::vector<std::string> read_lines;

    if (!read_file.is_open()) {
        std::cout << "Error opening input FASTA file" << std::endl;
        return 1;
    }

    while (std::getline(read_file, line)) {

        read_lines.push_back(line);
        if (line[0] == '>') {
            num_samples++;
            sample_names.push_back(line.substr(1));
        } else {
            aligned_reads.push_back(line);
        }
    }
    
    read_file.close();

    // Iterate over sample names and get start and end positions
    for (int i = 0; i < num_samples; i++) {
        std::string sample_name = sample_names[i];
        std::vector<std::string> tokens = split(sample_name, '_');
        int end = std::stoi(tokens.back());
        tokens.pop_back();
        int start = std::stoi(tokens.back());
        read_start_positions.push_back(start);
        read_end_positions.push_back(end);
    }

    std::string ref_seq;
    std::string chromosome_name;

    for (std::vector<std::string>::const_iterator it = ref_lines.begin(); it != ref_lines.end(); ++it) {
        const std::string& line = *it;
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
            
            // Check if i is within the start and end positions of the read
            if (i >= read_start_positions[k] && i < read_end_positions[k]) {
                
                // std::cout << i - read_start_positions[k] + 1 << std::endl;
                // std::cout << "ref: " << ref_seq[i] << " read: " << aligned_reads[k][i - read_start_positions[k] + 1] << std::endl;
                if (ref_seq[i] != aligned_reads[k][i - read_start_positions[k] + 1]) {
                mismatch = true;
                mismatch_indices.push_back(k);
            }
            
            }
            
        }
        
        if (mismatch) {
            
            std::set<char> mutations;
            for (std::vector<int>::const_iterator it = mismatch_indices.begin(); it != mismatch_indices.end(); ++it) {
                int index = *it;
                mutations.insert(aligned_reads[index][i - read_start_positions[index] + 1]);
            }

            for (std::set<char>::const_iterator it = mutations.begin(); it != mutations.end(); ++it) {
                char mutation = *it;
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
                    if (i >= read_start_positions[j] && i < read_end_positions[j]){
                    if (aligned_reads[j][i - read_start_positions[j] + 1] == mutation) {
                      
                        vcf_file << "1\t";
                    } }
                    else {
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
