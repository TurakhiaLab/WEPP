#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: ./vcf_to_fasta <vcf_file> <reference_fasta> <output_fasta>" << std::endl;
        return 1;
    }

    std::string reference = argv[2];
    std::string vcf_file = argv[1];
    std::string output_file = argv[3];

    std::ifstream ref_fasta(reference);
    std::string line, ref_seq;
    getline(ref_fasta, line); // Skip first line

    while (getline(ref_fasta, line)) {
        ref_seq += line;
    }
    ref_fasta.close();

    std::ifstream vcf_file_stream(vcf_file);
    std::vector<std::string> sample_names;
    while (getline(vcf_file_stream, line)) {
        if (line.find("#CHROM") == 0) {
            std::istringstream iss(line);
            std::string tmp;
            for (size_t i = 0; i < 9; ++i) iss >> tmp;
            while (iss >> tmp) sample_names.push_back(tmp);
            break;
        }
    }

    std::map<std::string, std::string> mutated_seqs;
    std::vector<int> start_positions(sample_names.size()), end_positions(sample_names.size());
    for (size_t sample_index = 0; sample_index < sample_names.size(); ++sample_index) {
        const std::string& sample_name = sample_names[sample_index];
        int start_position = 0;
        int end_position = ref_seq.size();
        //auto last_underscore = sample_name.find_last_of('_');
        //auto second_to_last_underscore = sample_name.find_last_of('_', last_underscore - 1);
        //if (second_to_last_underscore != std::string::npos) {
	//	std::cout << "Sample name: " << sample_name << std::endl;
        //    start_position = stoi(sample_name.substr(second_to_last_underscore + 1, last_underscore - second_to_last_underscore - 1)) - 1;
        //    end_position = stoi(sample_name.substr(last_underscore + 1));
        // }

        start_positions[sample_index] = start_position;
        end_positions[sample_index] = end_position;
       
        mutated_seqs[sample_name] = std::string(ref_seq.begin() + start_position, ref_seq.begin() + end_position);
    }

    while (getline(vcf_file_stream, line)) {
        std::istringstream iss(line);
        std::string chrom, pos_str, id, ref_base, alt_base;
        iss >> chrom >> pos_str >> id >> ref_base >> alt_base;
	//std::cout << "Pos str: " << pos_str << std::endl;
        int pos = stoi(pos_str);
        if (alt_base.find(',') != std::string::npos) {
            alt_base = alt_base.substr(0, alt_base.find(','));
        }

        std::string sample_info;
        for (size_t sample_index = 0; sample_index < sample_names.size(); ++sample_index) {
            iss >> sample_info;
            if (sample_info == "1" && pos >= start_positions[sample_index] && pos <= end_positions[sample_index]) {
                if (mutated_seqs[sample_names[sample_index]][pos - start_positions[sample_index] - 1] == ref_base[0]) {
                    mutated_seqs[sample_names[sample_index]][pos - start_positions[sample_index] - 1] = alt_base[0];
                } else {
                    std::cout << "Reference base (" << ref_base << ") does not match at position " << pos << " for sample " << sample_names[sample_index] << std::endl;
                }
            }
        }
    }

    std::ofstream mutated_fasta(output_file.c_str());
    for (const auto& it : mutated_seqs) {
        mutated_fasta << ">" << it.first << std::endl;
        mutated_fasta << it.second << std::endl;
    }
    mutated_fasta.close();

    return 0;
}
