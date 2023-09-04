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
    if (argc != 3) {
        std::cout << "Usage: group_vcf <input_vcf> <output_vcf>" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    // Read in vcf lines, except headers
    std::vector<std::string> vcf_lines;
    std::vector<std::string> vcf_headers;
    std::vector<std::string> sample_names;
    std::ifstream vcf_file(input_file);
    std::string line;
    while (std::getline(vcf_file, line)) {
        if (line[0] != '#') {
            // split the vcf line into words, and check if ref and alt are the same
            std::vector<std::string> words;
            string_split(line, '\t', words);
            
            vcf_lines.push_back(line);

            if (words[3] == words[4]) {
                std::cout << "Warning: ref and alt are the same at position " << words[1] << std::endl;
                // print the sample names for which this is the case
                std::cout << "Samples: " << std::endl;
                for (int i = 9; i < (int)words.size(); i++) {
                    if (words[i] == "1") {
                        std::cout << sample_names[i - 9] << " ";
                    }
                }
                std::cout << std::endl;


            }
            
        }
        else {
            vcf_headers.push_back(line);

            // if line starts with #CHROM, then store the sample names
            if (line.substr(0, 5) == "#CHRO") {
                std::vector<std::string> words;
                string_split(line, '\t', words);
                for (int i = 9; i < (int)words.size(); i++) {
                    sample_names.push_back(words[i]);
                }
            }
        }
    }
    vcf_file.close();

    // Group vcf lines by position
    std::vector<std::string> vcf_lines_merged;

    std::cout << "vcf_lines size: " << vcf_lines.size() << std::endl;

    for (int i = 0; i < (int)vcf_lines.size(); i++) {
       
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

    // Write to output file
    std::ofstream output(output_file);

    for (int i = 0; i < (int)vcf_headers.size(); i++) {
        output << vcf_headers[i] << std::endl;
    }


    for (int i = 0; i < (int)vcf_lines_merged.size(); i++) {
        output << vcf_lines_merged[i] << std::endl;
    }

    output.close();

    return 0;

}