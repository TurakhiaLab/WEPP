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
    if (argc != 5) {
        std::cout << "Usage: generate_freyja_files <input_vcf> <reference_fasta> <output_freyja_vcf> <output_depth_file>" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string reference_file = argv[2];
    std::string output_file = argv[3];
    std::string output_depth_file = argv[4];

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
    output_vcf << "##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

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

    // group the vcf lines by position, making a vector of vectors
    std::vector<std::vector<std::string>> vcf_lines_grouped;
    std::vector<std::string> vcf_lines_grouped_temp;

    for (int i = 0; i < (int)vcf_file_reads_lines_list.size(); i++) {
        std::vector<std::string> words;
        string_split(vcf_file_reads_lines_list[i], '\t', words);
        int pos = std::stoi(words[1]);
        vcf_lines_grouped_temp.push_back(vcf_file_reads_lines_list[i]);
        for (int j = i + 1; j < (int)vcf_file_reads_lines_list.size(); j++) {
            std::vector<std::string> words2;
            string_split(vcf_file_reads_lines_list[j], '\t', words2);
            int pos2 = std::stoi(words2[1]);
            if (pos == pos2) {
                vcf_lines_grouped_temp.push_back(vcf_file_reads_lines_list[j]);
                i = j + 1;
            } else {
                break;
            }
        }
        vcf_lines_grouped.push_back(vcf_lines_grouped_temp);
        vcf_lines_grouped_temp.clear();
    }

    // Make a vector of strings to store the lines for the VCF file for reads for Freyja
    std::vector<std::string> vcf_lines_freyja_list;

    for (auto group: vcf_lines_grouped) {
        std::vector<std::string> alts;
        std::vector<std::string> ids;
        std::vector<std::string> infos;
        std::string current_vcf_line = "";

        // store chrom, pos, ref, qual, filter for the first line in the group
        std::vector<std::string> words_first;
        string_split(group[0], '\t', words_first);

        for (auto line: group) {
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
            infos.push_back("AF=" + std::to_string(ratio));
            alts.push_back(words[4]);
            ids.push_back(words[2]);
        }

        // vcf_file_reads_freyja << words_first[0] << "\t" << words_first[1] << "\t";
        current_vcf_line += words_first[0] + "\t" + words_first[1] + "\t";
        for (int i = 0; i < (int)ids.size(); i++) {
            if (i == 0) {
                current_vcf_line += ids[i];
            } else {
                current_vcf_line += "," + ids[i];
            }
        }
        current_vcf_line += "\t" + words_first[3] + "\t";
        for (int i = 0; i < (int)alts.size(); i++) {
            if (i == 0) {
                current_vcf_line += alts[i];
            } else {
                current_vcf_line += "," + alts[i];
            }
        }
        current_vcf_line += "\t.\t.\t";
        for (int i = 0; i < (int)infos.size(); i++) {
            if (i == 0) {
                current_vcf_line += infos[i];
            } else {
                current_vcf_line += ";" + infos[i];
            }
        }
        current_vcf_line += "\n";
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
    // Make a variable to store all the lines for the depth file, so we can write them all at once
    std::vector<std::string> depth_lines;
    for (int i = 0; i < (int)ref_seq.size(); i++) {
        int read_count = 0;
        for (auto sample_read: sample_read_list) {
            int start_pos = std::stoi(sample_read.substr(sample_read.find("_READ_") + 6, sample_read.find("_", sample_read.find("_READ_") + 6) - (sample_read.find("_READ_") + 6)));
            int end_pos = std::stoi(sample_read.substr(sample_read.find("_", sample_read.find("_READ_") + 6) + 1, sample_read.size() - (sample_read.find("_", sample_read.find("_READ_") + 6) + 1)));
            if ((i >= start_pos) && (i <= end_pos)) {
                read_count++;
            }
        }
        depth_lines.push_back(reference + "\t" + std::to_string(i+1) + "\t" + ref_seq[i] + "\t" + std::to_string(read_count));
    }

    // combine the lines into one string, and write all the lines at once
    std::string depth_lines_combined;
    for (auto line: depth_lines) {
        depth_lines_combined += line + "\n";
    }

    output_depth << depth_lines_combined;
    output_depth.close();
}