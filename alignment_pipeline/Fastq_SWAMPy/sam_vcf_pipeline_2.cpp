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

struct Variant {
    std::string ref_name;
    int position;
    std::string ref_base;
    std::map<std::string, int> alt_bases; // Changed to map to hold multiple alternate alleles with their encoding
    std::map<std::string, std::set<std::string>> read_names; // Map from alt_bases to read names
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

std::vector<std::tuple<int, std::string, std::string>> parse_mdz(const std::string& mdz, const std::string& cigar, const std::string& seq, int pos, const std::string& ref_genome) {
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

std::tuple<std::string, std::string, int, std::string, std::string, std::string> parse_sam_line(const std::string& line) {
    std::istringstream iss(line);
    std::string read_name, flag, ref_name, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tag, tags_str;
    iss >> read_name >> flag >> ref_name >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;
    std::map<std::string, std::string> tags;
    while (iss >> tag) {
        if (tag.rfind("MD:Z:", 0) == 0) {
            tags.emplace("MD", tag.substr(5));
            break;
        }
    }
    return {read_name, ref_name, std::stoi(pos), seq, tags["MD"], cigar};
}

int main(int argc, char* argv[]) {
    if (argc != 5 && argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <sam_file> <reference_file> <vcf_file> <num_threads>" << std::endl;
        return 1;
    }

    std::string sam_file = argv[1];
    std::string reference_file = argv[2];
    std::string vcf_file = argv[3];
    int num_threads = (argc == 6) ? std::atoi(argv[4]) : 1;

    omp_set_num_threads(num_threads); // Set the number of threads for OpenMP

    std::string ref_genome = read_reference(reference_file);

    std::vector<std::string> lines;
    std::string line;
    std::ifstream file(sam_file);
    while (std::getline(file, line)) {
        if (line.rfind('@', 0) == 0) continue; // Skip header
        lines.push_back(line);
    }

    std::map<std::pair<std::string, int>, Variant> variants_map;
    std::set<std::string> reads;

    #pragma omp parallel
    {
        std::map<std::pair<std::string, int>, Variant> local_variants_map;
        std::set<std::string> local_reads;

        #pragma omp for schedule(dynamic) nowait
        for (size_t i = 0; i < lines.size(); ++i) {
            auto [read_name, ref_name, pos, seq, mdz, cigar] = parse_sam_line(lines[i]);
            auto mismatches = parse_mdz(mdz, cigar, seq, pos, ref_genome);

            local_reads.insert(read_name);
            for (const auto& [position, ref_base, alt_base] : mismatches) {
                std::pair<std::string, int> key = {ref_name, position};
                auto& var = local_variants_map[key];
                var.ref_name = ref_name;
                var.position = position;
                var.ref_base = ref_base;

                // Encode each variant with a unique number at this position
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

    std::ofstream out(vcf_file);
    out << "##fileformat=VCFv4.2\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    for (const auto& read : reads) {
        out << read << '\t';
    }
    out << '\n';

    for (const auto& [key, var] : variants_map) {
        out << var.ref_name << '\t'
            << var.position << '\t'
            << '.' << '\t' // ID is '.' if unknown
            << var.ref_base << '\t';

        // Concatenate all alternate alleles
        std::string alts;
        for (const auto& [alt_base, _] : var.alt_bases) {
            if (!alts.empty()) alts += ",";
            alts += alt_base;
        }
        out << alts << "\t.\t.\t.\tGT\t";

        for (const auto& read : reads) {
            int alt_code = 0;
            for (const auto& [alt_base, code] : var.alt_bases) {
                if (var.read_names.at(alt_base).find(read) != var.read_names.at(alt_base).end()) {
                    alt_code = code; // Get the code of the alt allele present for this read
                    break;
                }
            }
            out << alt_code << '\t'; // Write the code (0 if not present, otherwise the allele code)
        }
        out << '\n';
    }

    out.close();
    return 0;
}