#include <iostream>
#include <fstream>
#include <vector>
#include <map>

int main(int argc, char* argv[]) {
    // Check that correct number of arguments are provided
    if (argc != 4) {
        std::cout << "Usage: ./vcf_to_fasta <vcf_file> <reference_fasta> <output_fasta>" << std::endl;
        return 1;
    }

    std::string reference = argv[2];
    std::string vcf_file = argv[1];
    std::string output_file = argv[3];

    // Read reference fasta file into memory
    std::ifstream ref_fasta(reference);
    std::vector<std::string> ref_lines;
    std::string line;
    while (std::getline(ref_fasta, line)) {
        ref_lines.push_back(line);
    }
    ref_fasta.close();
    std::string ref_seq = "";
    for (size_t i = 1; i < ref_lines.size(); i++) {
        ref_seq += ref_lines[i];
    }

    // Read VCF file into memory
    std::ifstream vcf_file_stream(vcf_file);
    std::vector<std::string> vcf_lines;
    while (std::getline(vcf_file_stream, line)) {
        vcf_lines.push_back(line);
    }
    vcf_file_stream.close();

    std::string header_line = "";
    // Get list of sample names from VCF file header
    for (size_t i = 0; i < vcf_lines.size(); i++) {
        const std::string& line = vcf_lines[i];
        if (line.find("#CHROM") == 0) {
            header_line = line;
            break;
        }
    }

    std::vector<std::string> headers;
    size_t start_pos = 0;
    size_t end_pos = header_line.find('\t');
    while (end_pos != std::string::npos) {
        headers.push_back(header_line.substr(start_pos, end_pos - start_pos));
        start_pos = end_pos + 1;
        end_pos = header_line.find('\t', start_pos);
    }
    headers.push_back(header_line.substr(start_pos));

    std::vector<std::string> sample_names(headers.begin() + 9, headers.end());

    // Create a directory called intermediate_files if it does not exist
    system("mkdir -p ./intermediate_files");

    // Save all the sample names in a file
    std::ofstream sample_names_file("./intermediate_files/sample_names.txt");
    for (size_t j = 0; j < sample_names.size() - 1; j++) {
        sample_names_file << sample_names[j] << ',';
    }
    sample_names_file << sample_names[sample_names.size() - 1];
    sample_names_file.close();

    // Create dictionary to store mutated sequences
    std::map<std::string, std::string> mutated_seqs;

    // Loop through each sample in VCF file
    for (size_t sample_index = 0; sample_index < sample_names.size(); sample_index++) {
        const std::string& sample_name = sample_names[sample_index];
        std::vector<char> mutated_seq(ref_seq.begin(), ref_seq.end());

        // Loop through each variant row in VCF file for the current sample
        for (size_t i = 4; i < vcf_lines.size(); i++) {
            const std::string& vcf_line = vcf_lines[i];
            std::vector<std::string> vcf_fields;
            start_pos = 0;
            end_pos = vcf_line.find('\t');
            while (end_pos != std::string::npos) {
                vcf_fields.push_back(vcf_line.substr(start_pos, end_pos - start_pos));
                start_pos = end_pos + 1;
                end_pos = vcf_line.find('\t', start_pos);
            }
            vcf_fields.push_back(vcf_line.substr(start_pos));

            const std::string& chrom = vcf_fields[0];
            int pos = std::stoi(vcf_fields[1]);
            const std::string& ref_base = vcf_fields[3];
            std::string alt_base = vcf_fields[4];
            if (alt_base.find(',') != std::string::npos) {
                alt_base = alt_base.substr(0, alt_base.find(','));
            }

            const std::string& gt = vcf_fields[9 + sample_index];

            // Check if sample is mutated at current position
            if (gt == "1") {
                // Check that reference base at current position matches expected base
                if (mutated_seq[pos - 1] == ref_base[0]) {
                    mutated_seq[pos - 1] = alt_base[0];
                } else {
                    std::cout << "Reference base (" << ref_base << ") does not match at position " << pos
                              << " for sample " << sample_name << std::endl;
                }
            }
        }

        // Create mutated sequence and store in dictionary
        std::string mutated_seq_str(mutated_seq.begin(), mutated_seq.end());
        mutated_seqs[sample_name] = mutated_seq_str;
    }

    // Write mutated sequences to Fasta file
    std::ofstream mutated_fasta(output_file.c_str());
    for (std::map<std::string, std::string>::const_iterator it = mutated_seqs.begin(); it != mutated_seqs.end(); ++it) {
        mutated_fasta << ">" << it->first << std::endl;
        mutated_fasta << it->second << std::endl;
    }
    mutated_fasta.close();

    return 0;
}
