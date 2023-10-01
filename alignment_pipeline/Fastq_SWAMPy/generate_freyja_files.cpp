#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <sstream>

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
    if (argc != 6) {
        std::cout << "Usage: generate_freyja_files <input_vcf> <reference_fasta> <output_freyja_vcf> <output_depth_file> <input grouped vcf>" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string reference_file = argv[2];
    std::string output_file = argv[3];
    std::string output_depth_file = argv[4];
    std::string input_grouped_vcf = argv[5];

    // Read in the reference sequence and store it in a string called ref_seq
    std::ifstream reference_fasta(reference_file);
    std::string ref_seq;
    std::string line2;
    while (std::getline(reference_fasta, line2)) {
        if (line2[0] != '>') {
            ref_seq += line2;
        }
    }

    std::ifstream input_vcf(input_file);
    std::ofstream output_vcf(output_file);
    std::ofstream output_depth(output_depth_file);

    output_vcf << "##fileformat=VCFv4.2" << std::endl;
    output_vcf << "##source=alignmentpipeline" << std::endl;
    output_vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

    // read in vcf lines, except headers
    std::vector<std::string> vcf_file_reads_lines_list;
    std::vector<std::string> sample_read_list;
    // iterate through the vcf file, and store the lines in a vector, skipping the headers
    std::string line;
    while (std::getline(input_vcf, line)) {
        if (line[0] != '#') {
            vcf_file_reads_lines_list.push_back(line);
        }
        else {
            // if line starts with #CHROM, then store the sample names
            if (line.substr(0, 6) == "#CHROM") {
                std::vector<std::string> words;
                string_split(line, '\t', words);
                for (int i = 9; i < (int)words.size(); i++) {
                    sample_read_list.push_back(words[i]);
                }
            }
        }
    }

    // Make a vector of strings to store the lines for the VCF file for reads for Freyja
    std::vector<std::string> vcf_lines_freyja_list;

    for (auto line: vcf_file_reads_lines_list) {
        std::vector<std::string> alts;
        std::vector<std::string> ids;
        std::vector<std::string> infos;
        std::string current_vcf_line = "";

        // store chrom, pos, ref, qual, filter for the first line in the group
        std::vector<std::string> words;
        string_split(line, '\t', words);

        int num_ones = 0;
        
        for (int i = 9; i < (int)words.size(); i++) {
            if (words[i] != "0") {
                num_ones++;
            }
        }

        int num_reads = 0;
        for (auto sample_read: sample_read_list) {
            int start_pos = std::stoi(sample_read.substr(sample_read.find("_READ_") + 6, sample_read.find("_", sample_read.find("_READ_") + 6) - (sample_read.find("_READ_") + 6)));
            int end_pos = std::stoi(sample_read.substr(sample_read.find("_", sample_read.find("_READ_") + 6) + 1, sample_read.size() - (sample_read.find("_", sample_read.find("_READ_") + 6) + 1)));
            if ((std::stoi(words[1]) >= start_pos) && (std::stoi(words[1]) <= end_pos)) {
                num_reads++;
            }
        }
        
        double ratio = (double)num_ones / (double)num_reads;
        
        current_vcf_line +=  words[0] + "\t" + words[1] + "\t" + words[2] + "\t" + words[3] + "\t" + words[4] + "\t" + words[5] + "\t" + words[6] + "\t" + "AF=" + std::to_string(ratio) + "\n";
        vcf_lines_freyja_list.push_back(current_vcf_line);
    }

    // combine all the lines for VCF freyja into one string, and write all the lines at once
    std::string vcf_lines_freyja_combined;
    for (auto line: vcf_lines_freyja_list) {
        vcf_lines_freyja_combined += line;
    }

    output_vcf << vcf_lines_freyja_combined;
    output_vcf.close();

      std::string reference = "NC_045512v2";
    std::vector<std::string> depth_lines;
    std::vector<int> read_counts;
    for (int i = 1; i <= (int)ref_seq.size(); i++) {
        int read_count = 0;
        for (auto sample_read: sample_read_list) {
            int start_pos = std::stoi(sample_read.substr(sample_read.find("_READ_") + 6, sample_read.find("_", sample_read.find("_READ_") + 6) - (sample_read.find("_READ_") + 6)));
            int end_pos = std::stoi(sample_read.substr(sample_read.find("_", sample_read.find("_READ_") + 6) + 1, sample_read.size() - (sample_read.find("_", sample_read.find("_READ_") + 6) + 1)));
            if ((i >= start_pos) && (i <= end_pos)) {
                read_count++;
            }
        }
        read_counts.push_back(read_count);
        depth_lines.push_back(reference + "\t" + std::to_string(i) + "\t" + ref_seq[i-1] + "\t" + std::to_string(read_count));
    }

    // combine the lines into one string, and write all the lines at once
    std::string depth_lines_combined;
    for (auto line: depth_lines) {
        depth_lines_combined += line + "\n";
    }

    output_depth << depth_lines_combined;
    output_depth.close();

    // iterate through sample_read_list and for each read name, append _ and the average depth (across all positions) to the end of the read name
    std::vector<std::string> sample_read_list_depth;
    for (auto sample_read: sample_read_list) {
        int start_pos = std::stoi(sample_read.substr(sample_read.find("_READ_") + 6, sample_read.find("_", sample_read.find("_READ_") + 6) - (sample_read.find("_READ_") + 6)));
        int end_pos = std::stoi(sample_read.substr(sample_read.find("_", sample_read.find("_READ_") + 6) + 1, sample_read.size() - (sample_read.find("_", sample_read.find("_READ_") + 6) + 1)));
        int total_depth = 0;
        for (int i = start_pos; i <= end_pos; i++) {
            total_depth += read_counts[i-1];
        }
        int average_depth = total_depth / (end_pos - start_pos + 1);
        sample_read_list_depth.push_back(sample_read + "_" + std::to_string(average_depth));
    }

    // Read the input grouped VCF file
    std::ifstream input_grouped_vcf_file(input_grouped_vcf);

    // Read all lines
    std::vector<std::string> input_grouped_vcf_lines;
    std::string line3;
    while (std::getline(input_grouped_vcf_file, line3)) {
        input_grouped_vcf_lines.push_back(line3);
    }

    input_grouped_vcf_file.close();

    // Iterate through the lines and for the line that starts with #CHROM, edit the read names
    for (int i = 0; i < (int)input_grouped_vcf_lines.size(); i++) {
        if (input_grouped_vcf_lines[i].substr(0, 6) == "#CHROM") {
            std::vector<std::string> words;
            string_split(input_grouped_vcf_lines[i], '\t', words);
            for (int j = 9; j < (int)words.size(); j++) {
                words[j] = sample_read_list_depth[j-9];
            }
            std::string new_line = "";
            for (auto word: words) {
                new_line += word + "\t";
            }

            new_line = new_line.substr(0, new_line.size()-1); // remove the last tab

            input_grouped_vcf_lines[i] = new_line;

            break;
        }
    }

    // Combine all the lines into one string and write all the lines at once to the same grouped VCF file
    std::string input_grouped_vcf_lines_combined;
    for (auto line: input_grouped_vcf_lines) {
        input_grouped_vcf_lines_combined += line + "\n";
    }

    std::ofstream input_grouped_vcf_file2(input_grouped_vcf); // overwrite the file
    input_grouped_vcf_file2 << input_grouped_vcf_lines_combined;
    input_grouped_vcf_file2.close();

    return 0;

}