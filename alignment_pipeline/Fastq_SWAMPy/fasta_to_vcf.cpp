#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <sstream>

  // Split string into words for a specific delimiter delim
void string_split (std::string const& s, char delim, std::vector<std::string>& words) {
    size_t start_pos = 0, end_pos = 0;
    while ((end_pos = s.find(delim, start_pos)) != std::string::npos) {
        // if ((end_pos == start_pos) || end_pos >= s.length()) {
        if (end_pos >= s.length()) {
            break;
        }
        words.push_back(s.substr(start_pos, end_pos-start_pos));
        start_pos = end_pos+1;
    }
    auto last = s.substr(start_pos, s.size()-start_pos);
    if (last != "") {
        words.push_back(std::move(last));
    }
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
    // for (int i = 0; i < num_samples; i++) {
    //     std::string sample_name = sample_names[i];
    //     std::vector<std::string> tokens = split(sample_name, '_');
    //     int end = std::stoi(tokens.back());
    //     tokens.pop_back();
    //     int start = std::stoi(tokens.back());
    //     read_start_positions.push_back(start);
    //     read_end_positions.push_back(end);
    // }

    for (int i = 0; i < num_samples; i++) {
        std::string sample_name = sample_names[i];
        int start_pos = std::stoi(sample_name.substr(sample_name.find("_READ_") + 6, sample_name.find("_", sample_name.find("_READ_") + 6) - (sample_name.find("_READ_") + 6)));
        int end_pos = std::stoi(sample_name.substr(sample_name.find("_", sample_name.find("_READ_") + 6) + 1, sample_name.size() - (sample_name.find("_", sample_name.find("_READ_") + 6) + 1)));
        read_start_positions.push_back(start_pos);
        read_end_positions.push_back(end_pos);
    }

    // print out the start and end positions
    // for (int i = 0; i < num_samples; i++) {
    //     std::cout << sample_names[i] << " " << read_start_positions[i] << " " << read_end_positions[i] << std::endl;
    // }


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

    // print length of reference sequence
    std::cout << "Length of reference: "  << ref_seq.length() << std::endl;

    std::ofstream vcf_file(output_file);
    if (!vcf_file.is_open()) {
        std::cout << "Error opening output file" << std::endl;
        return 1;
    }

    vcf_file << "##fileformat=VCFv4.2\n";
    vcf_file << "##reference=stdin:hCoV-19/Wuhan/Hu-1/2019|EPI_ISL_402125|2019-12-31\n";
    vcf_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

    for (int i = 0; i < num_samples; i++) {
        vcf_file << sample_names[i] << "\t";
    }

    vcf_file << "\n";

    for (int i = 1; i <= ref_seq.length(); i++) {
        bool mismatch = false;
        std::vector<int> mismatch_indices;

        for (int k = 0; k < num_samples; k++) {
            
            // Check if i is within the start and end positions of the read
            if (i >= read_start_positions[k] && i <= read_end_positions[k]) {
                
                // std::cout << i - read_start_positions[k] + 1 << std::endl;
                // std::cout << "ref: " << ref_seq[i] << " read: " << aligned_reads[k][i - read_start_positions[k] + 1] << std::endl;
                if (ref_seq[i-1] != aligned_reads[k][i - read_start_positions[k]]) {
                mismatch = true;
                mismatch_indices.push_back(k);
            }
            
            }
            
        }
        
        if (mismatch) {
            
            std::set<char> mutations;
            for (std::vector<int>::const_iterator it = mismatch_indices.begin(); it != mismatch_indices.end(); ++it) {
                int index = *it;
                mutations.insert(aligned_reads[index][i - read_start_positions[index]]);
            }

            for (std::set<char>::const_iterator it = mutations.begin(); it != mutations.end(); ++it) {
                char mutation = *it;
                vcf_file << chromosome_name << "\t";
                vcf_file << i  << "\t";
                char ref = ref_seq[i-1];
                char alt = mutation;
                std::string id = std::string(1, ref) + std::to_string(i) + std::string(1, alt);
                vcf_file << id << "\t";
                vcf_file << ref << "\t";
                vcf_file << alt << "\t";
                vcf_file << ".\t";
                vcf_file << ".\t";
                vcf_file << ".\t";
                vcf_file << ".\t";

                for (int j = 0; j < num_samples; j++) {
                    if (i >= read_start_positions[j] && i <= read_end_positions[j]){
                    if (aligned_reads[j][i - read_start_positions[j]] == mutation) {
                      
                        vcf_file << "1\t";
                    }

                    else {
                        vcf_file << "0\t";} 
                    
                    }
                    else {
                        vcf_file << "0\t";
                    }
                   
                }

                vcf_file << "\n";
                
            }

            
        }
    }

    vcf_file.close();


    // Read all lines from vcf_file and store it in a variable called vcf_lines
    std::ifstream vcf_file2(output_file);
    std::vector<std::string> vcf_lines;

    if (!vcf_file2.is_open()) {
        std::cout << "Error opening output file" << std::endl;
        return 1;
    }

    while (std::getline(vcf_file2, line)) {
            vcf_lines.push_back(line);
    }

    vcf_file2.close();

    std::vector<std::string> vcf_lines_merged;

    std::cout << "vcf_lines size: " << vcf_lines.size() << std::endl;

    for (int i = 0; i < (int)vcf_lines.size(); i++) {
        if (vcf_lines[i][0] == '#') {
            vcf_lines_merged.push_back(vcf_lines[i]);
            continue;
        }

        std::vector<std::string> words;
        string_split(vcf_lines[i], '\t' , words);
        // print words
        // for (int i = 0; i < words.size(); i++) {
        //     std::cout << words[i] << std::endl;
        // }

        int pos = std::stoi(words[1]);
        std::string alt = words[4];
        std::string id = words[2];
        std::vector<int> sample_i_counts;
        
        for (int k = 9; k < (int)words.size(); k++) {
            sample_i_counts.push_back(std::stoi(words[k]));
        }

        int j;

        for (j = i + 1; j < (int)vcf_lines.size(); j++) {
            std::vector<std::string> words2;
            string_split(vcf_lines[j], '\t', words2);
            int pos2 = std::stoi(words2[1]);
            std::string alt2 = words2[4];
            std::string id2 = words2[2];
            if (pos == pos2) {
                alt = alt + "," + alt2;
                id = id + "," + id2;
                for (int k = 9; k < (int)words2.size(); k++) {
                    if (words[k] == "1") {
                        sample_i_counts[k - 9] = 1;
                    }
                    if (words2[k] == "1") {
                        sample_i_counts[k - 9] = j - i + 1;
                    }
                
            }
            } 
            
            else {
                break;
            }
        }

        i = j - 1;


        vcf_lines_merged.push_back(words[0] + "\t" + words[1] + "\t" + id + "\t" + words[3] + "\t" + alt + "\t" + words[5] + "\t" + words[6] + "\t" + words[7] + "\t" + words[8]);
        for (int k = 0; k < (int)sample_i_counts.size(); k++) {
            vcf_lines_merged.back() += "\t" + std::to_string(sample_i_counts[k]);
        }
    }

    // Write the merged vcf lines, overwriting the original vcf file
    std::ofstream vcf_file3(output_file);
    if (!vcf_file3.is_open()) {
        std::cout << "Error opening output file" << std::endl;
        return 1;
    }

    for (int i = 0; i < (int)vcf_lines_merged.size(); i++) {
        vcf_file3 << vcf_lines_merged[i] << "\n";
    }

    vcf_file3.close();

    return 0;
}
