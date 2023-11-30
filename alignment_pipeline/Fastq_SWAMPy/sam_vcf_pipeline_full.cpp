#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <regex>
#include <thread>
#include <mutex>
#include <cstdlib>
#include <omp.h>
#include <tuple>
#include <cctype>

struct Variant {
    std::string ref_name;
    int position;
    std::string ref_base;
    std::map<std::string, int> alt_bases;
    std::map<std::string, std::set<std::string>> read_names;
};

std::string read_reference(const std::string& reference_file) {
    std::ifstream file(reference_file);
    std::string line, ref_genome = "";
    
    while (std::getline(file, line)) {
        if (line.rfind('>', 0) == 0) continue; // Skip line if it starts with '>'
        ref_genome += line;
    }
    
    file.close();
    return ref_genome;
}

std::vector<std::tuple<int, std::string, std::string>> parse_mdz(const std::string& mdz, const std::string& cigar, const std::string& seq, int pos, const std::string& ref_genome, bool no_indels) {
    std::vector<std::tuple<int, std::string, std::string>> mismatches;
    std::regex token_regex(R"((\d+|\^[A-Z]+|[A-Z]))");
    std::sregex_iterator tokens(mdz.begin(), mdz.end(), token_regex), end;

    int read_pos = 0, ref_pos = pos;

    std::regex cigar_regex(R"((\d+)([MIDNSHPX=]))");
    std::sregex_iterator cigar_tokens(cigar.begin(), cigar.end(), cigar_regex);

    for (; tokens != end; ++tokens) {
        std::string token = (*tokens)[1];
        if (token.find('^') != std::string::npos) { // deletion
   
            int del_len = token.length() - 1;
       
            ref_pos += del_len;
            mismatches.emplace_back(ref_pos, token.substr(1), "-");
        } else if (std::isdigit(token[0])) { // matching bases
        
            int length = std::stoi(token);
      
            read_pos += length;
            ref_pos += length;
        } else { // mismatch
            char ref_base = ref_genome[ref_pos - 1];
            char alt_base = seq[read_pos];
       
            mismatches.emplace_back(ref_pos, std::string(1, ref_base), std::string(1, alt_base));
            read_pos++;
            ref_pos++;
        }
    }

    for (; cigar_tokens != end; ++cigar_tokens) {
        int count = std::stoi((*cigar_tokens)[1]);
        char operation = (*cigar_tokens)[2].str()[0];
        if (operation == 'I') { // Insertion
            std::string alt_base = seq.substr(read_pos, count);
            mismatches.emplace_back(ref_pos, "-", alt_base);
            read_pos += count;
        } else if (operation == 'D' || operation == 'N') { // Deletion
            ref_pos += count;
        }
    }

    return mismatches;
}

std::tuple<std::string, std::string, int, std::string, std::string> parse_sam_line(const std::string& line) {
    std::istringstream iss(line);
    std::string read_name, flag, ref_name, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tag, tags_str;
    iss >> read_name >> flag >> ref_name >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;

    return {read_name, ref_name, std::stoi(pos), seq,  cigar};
}

std::vector<std::tuple<int, std::string, std::string>> parse_cigar(const std::string& cigar, const std::string& seq, int pos, const std::string& ref_genome) {
    std::vector<std::tuple<int, std::string, std::string>> mutations;
    std::regex cigar_regex(R"((\d+)([MIDNSHPX=]))");
    std::sregex_iterator cigar_tokens(cigar.begin(), cigar.end(), cigar_regex);

    int read_pos = 0, ref_pos = pos - 1; // Adjust for 0-based index

    for (; cigar_tokens != std::sregex_iterator(); ++cigar_tokens) {
        int count = std::stoi((*cigar_tokens)[1]);
        char operation = (*cigar_tokens)[2].str()[0];

        switch (operation) {
            case 'M': { // Match or Mismatch
                for (int i = 0; i < count; ++i) {
                    char ref_base = ref_genome[ref_pos + i];
                    char alt_base = seq[read_pos + i];
                    if (ref_base != alt_base) { // Mismatch
                        mutations.emplace_back(ref_pos + i + 1, std::string(1, ref_base), std::string(1, alt_base));
                    }
                }
                read_pos += count;
                ref_pos += count;
                break;
            }
            case 'I': { // Insertion
                std::string alt_base = seq.substr(read_pos, count);
                mutations.emplace_back(ref_pos + 1, "-", alt_base);
                read_pos += count;
                break;
            }
            case 'D': { // Deletion
                std::string ref_base = ref_genome.substr(ref_pos, count);
                mutations.emplace_back(ref_pos + 1, ref_base, "-");
                ref_pos += count;
                break;
            }
            // Other CIGAR operations can be handled as needed
        }
    }

    return mutations;
}

// Function to extract start and end positions from the read name
std::pair<int, int> parse_read_positions(const std::string& read_name) {
    std::pair<int, int> positions(0, 0);
    
    size_t readPos = read_name.rfind("_READ_");
    
    if (readPos == std::string::npos) {
        std::cerr << "Invalid format: '_READ_' not found" << std::endl;
        return positions;
    }

    std::string numbersPart = read_name.substr(readPos + 6); // 6 is length of "_READ_"
    size_t underscorePos = numbersPart.find('_');
    if (underscorePos == std::string::npos) {
        std::cerr << "Invalid format: second '_' not found" << std::endl;
        return positions;
    }

    positions.first = std::stoi(numbersPart.substr(0, underscorePos));
    positions.second = std::stoi(numbersPart.substr(underscorePos + 1));

    return positions;
}

// Function to append the average depth to the end of the read name
std::string append_depth(const std::string& read_name, int depth) {
    std::string depth_str = std::to_string(depth);
    return read_name + "_" + depth_str;
}

// Function to find the average depth of a read (i.e. the average value of depth for all positions covered by the read)
int find_read_depth(const std::string& read_name, const std::vector<int>& depth) {
    auto positions = parse_read_positions(read_name);
    int read_depth = 0;
    for (int i = positions.first; i <= positions.second && i <= static_cast<int>(depth.size()); ++i) {
        read_depth += depth[i - 1];
    }
    return read_depth / (positions.second - positions.first + 1);
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <sam_file> <reference_file> <vcf_file> <freyja_vcf_file> <freyja_depth_file> <num_threads> [--no-indels]" << std::endl;
        return 1;
    }

    // Optional argument --no-indels
    bool no_indels = false;
    if (argc == 8) {
        std::string arg = argv[7];
        if (arg == "--no-indels") {
            no_indels = true;
        } else {
            std::cerr << "Usage: " << argv[0] << " <sam_file> <reference_file> <vcf_file> <freyja_vcf_file> <freyja_depth_file> <num_threads> [--no-indels]" << std::endl;
            return 1;
        }
    }

    std::string sam_file = argv[1];
    std::string reference_file = argv[2];
    std::string vcf_file = argv[3];
    std::string freyja_vcf_file = argv[4];
    std::string freyja_depth_file = argv[5];
    int num_threads = std::atoi(argv[6]);

    omp_set_num_threads(num_threads); // Set the number of threads for OpenMP

    std::string ref_genome = read_reference(reference_file);

    std::vector<std::string> lines;
    std::string line;
    std::ifstream file(sam_file);
    while (std::getline(file, line)) {
        if (line.rfind('@', 0) == 0) continue; // Skip header
        lines.push_back(line);
    }

    // Print the number of lines in the SAM file, which is the number of reads
    std::cout << "Number of lines in SAM file: " << lines.size() << std::endl;

    std::map<std::pair<std::string, int>, Variant> variants_map;
    std::set<std::string> reads;
    std::vector<int> depth(ref_genome.size(), 0); // Initialize depth vector with 0s

    #pragma omp parallel
    {
        std::map<std::pair<std::string, int>, Variant> local_variants_map;
        std::set<std::string> local_reads;

        #pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < lines.size(); ++i) {
            // auto [read_name, ref_name, pos, seq, mdz, cigar] = parse_sam_line(lines[i]);
            // auto mismatches = parse_mdz(mdz, cigar, seq, pos, ref_genome, no_indels);
            auto [read_name, ref_name, pos, seq, cigar] = parse_sam_line(lines[i]);
            auto mismatches = parse_cigar(cigar, seq, pos, ref_genome);

            local_reads.insert(read_name); // since this is a set, duplicates will be ignored
            for (const auto& [position, ref_base, alt_base] : mismatches) {
                std::pair<std::string, int> key = {ref_name, position};
                auto& var = local_variants_map[key];
                var.ref_name = ref_name;
                var.position = position;
                var.ref_base = ref_base;

                // Encode each variant with a unique number at this position
                // e.g. A -> 1, C -> 2, G -> 3, T -> 4
                int code = var.alt_bases.size() + 1;
                if (var.alt_bases.find(alt_base) == var.alt_bases.end()) {
                    var.alt_bases[alt_base] = code;
                }

                // Insert read names for this specific alternate allele
                var.read_names[alt_base].insert(read_name);
            }
        }

        #pragma omp critical
        {
            reads.insert(local_reads.begin(), local_reads.end());
            for (const auto& kv : local_variants_map) {
                auto& var = variants_map[kv.first];
                if (var.read_names.empty()) {
                    var = kv.second;
                } else {
                    for (const auto& alt_kv : kv.second.alt_bases) {
                        const auto& alt_base = alt_kv.first;
                        int alt_code = alt_kv.second;
                        var.alt_bases[alt_base] = alt_code;
                        var.read_names[alt_base].insert(kv.second.read_names.at(alt_base).begin(),
                                                        kv.second.read_names.at(alt_base).end());
                    }
                }
            }
        }
    }

    // Generate Freyja VCF and depth files
    std::ofstream freyja_vcf_out(freyja_vcf_file);
    std::ofstream freyja_depth_out(freyja_depth_file);

    // After parallel region, calculate depth for each position
    for (const auto& read : reads) {
        auto positions = parse_read_positions(read);
        for (int i = positions.first; i <= positions.second && i <= static_cast<int>(ref_genome.size()); ++i) {
            depth[i - 1]++;
        }
    }

    std::cout << "Writing Freyja VCF file..." << std::endl;

        // Write headers for Freyja VCF
    freyja_vcf_out << "##fileformat=VCFv4.2\n";
    freyja_vcf_out << "##source=alignmentpipeline\n";
    freyja_vcf_out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    // Process each variant for Freyja VCF and depth file entries
    for (const auto& [key, var] : variants_map) {
        for (const auto& [alt_base, code] : var.alt_bases) {
            // if the ref base equals the alt base, print a message and skip
            if (var.ref_base == alt_base) {
                // std::cout << "Ref base equals alt base at " << var.ref_name << ":" << var.position << std::endl;
                // std::cout << "Ref base: " << var.ref_base << std::endl;
                // std::cout << "Alt base: " << alt_base << std::endl;
                continue;
            }

            if (no_indels && (var.ref_base == "-" || alt_base == "-")) continue; // Skip indels if --no-indels flag is set
            freyja_vcf_out << var.ref_name << '\t'
                           << var.position << '\t'
                           << var.ref_base << var.position << alt_base <<'\t'
                           << var.ref_base << '\t'
                           << alt_base << '\t'
                           << '.' << '\t'
                           << '.' << '\t'
                           << "AF=" << std::to_string(static_cast<double>(var.read_names.at(alt_base).size()) / depth[var.position - 1]) << '\n';
        }
    }

    std::cout << "Writing Freyja depth file..." << std::endl;

    // Write the depth information for each position into Freyja depth file
    for (size_t i = 0; i < ref_genome.size(); ++i) {
        freyja_depth_out << "NC_045512v2" << '\t' << (i + 1) << '\t'
                         << ref_genome[i] << '\t' << depth[i] << '\n';
    }

    freyja_vcf_out.close();
    freyja_depth_out.close();

    std::cout << "Writing WEPP VCF file..." << std::endl;

    std::ofstream out(vcf_file);
    std::string vcf_file_contents = "";
    vcf_file_contents += "##fileformat=VCFv4.2\n";
    vcf_file_contents += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    for (const auto& read : reads) {
        vcf_file_contents += append_depth(read, find_read_depth(read, depth)) + '\t';
    }
    vcf_file_contents += '\n';

    for (const auto& [key, var] : variants_map) {
        if (no_indels && (var.ref_base == "-")) continue; // Skip indels if --no-indels flag is set
        // Concatenate all alternate alleles
        std::string alts;
        std::vector<std::string> alt_list;
        for (const auto& [alt_base, _] : var.alt_bases) {
            if (var.ref_base == alt_base) continue; // Skip if ref base equals alt base (i.e. no variant)
            if (no_indels && (alt_base == "-")) continue; // Skip indels if --no-indels flag is set
            if (!alts.empty()) alts += ",";
            alts += alt_base;
            alt_list.push_back(alt_base);
        }

        vcf_file_contents += var.ref_name + '\t' + std::to_string(var.position) + '\t';
        for (const auto& alt_base : alt_list) {
            vcf_file_contents += var.ref_base + std::to_string(var.position) + alt_base + ",";
        }

        vcf_file_contents.pop_back(); // Remove the last comma
        
        vcf_file_contents += '\t' + var.ref_base + '\t';
        vcf_file_contents += alts + "\t.\t.\t.\t.\t";

        for (const auto& read : reads) {
            int alt_code = 0;
            for (const auto& [alt_base, code] : var.alt_bases) {
                if (var.ref_base == alt_base) continue; // Skip if ref base equals alt base (i.e. no variant)
                if (no_indels && (alt_base == "-")) continue; // Skip indels if --no-indels flag is set
                if (var.read_names.at(alt_base).find(read) != var.read_names.at(alt_base).end()) {
                    alt_code = code; // Get the code of the alt allele present for this read
                    break;
                }
            }
            vcf_file_contents += std::to_string(alt_code) + '\t'; // Write the code (0 if not present, otherwise the allele code)
        }
        vcf_file_contents += '\n';
    }

    out << vcf_file_contents;
    out.close();

    return 0;
}