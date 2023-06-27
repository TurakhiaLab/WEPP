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
    size_t first_variant_line = vcf_lines.size();
    // Get list of sample names from VCF file header
    for (size_t i = 0; i < vcf_lines.size(); i++) {
        const std::string& line = vcf_lines[i];
        if (line.find("#CHROM") == 0) {
            header_line = line;
            first_variant_line = i;
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

    // Create dictionary to store mutated sequences
    std::map<std::string, std::string> mutated_seqs;

    // Loop through each sample in VCF file
    for (size_t sample_index = 0; sample_index < sample_names.size(); sample_index++) {
        const std::string& sample_name = sample_names[sample_index];

        int start_position = 0;
        int end_position = 0;

        // check if there are two underscores in the sample name
        if (sample_name.find('_', sample_name.find('_') + 1) == std::string::npos) {
            // std::cout << "Sample name " << sample_name << " does not have two underscores" << std::endl;
            start_position = 0;
            end_position = ref_seq.end() - ref_seq.begin();
        }
        else {
            // std::cout << "Sample name " << sample_name << " has two underscores" << std::endl;
            // Now, isolate the start and end positions; the sample name is of the form <sample Name>_READ_<start position>_<end position>
            start_position = std::stoi(sample_name.substr(sample_name.find('_', sample_name.find('_') + 1) + 1, sample_name.find('_', sample_name.find('_', sample_name.find('_') + 1) + 1) - sample_name.find('_', sample_name.find('_') + 1) - 1)) - 1;
            end_position = std::stoi(sample_name.substr(sample_name.find('_', sample_name.find('_', sample_name.find('_') + 1) + 1) + 1));
            
        }

        // Now, initialize the mutated read to be the reference sequence, from the start position to the end position
        std::vector<char> mutated_seq(ref_seq.begin() + start_position, ref_seq.begin() + end_position);

        // Loop through each variant row in VCF file for the current sample
        for (size_t i = first_variant_line + 1; i < vcf_lines.size(); i++) {
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

            // Check if sample is mutated at current position and if the position is within the sample's region
            if (gt == "1" && pos >= start_position && pos <= end_position) {
                // Check that reference base at current position matches expected base
                // need to shift pos to account for start position of sample
                
                if (mutated_seq[pos - start_position - 1] == ref_base[0]) {
                    mutated_seq[pos - start_position - 1] = alt_base[0];
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
