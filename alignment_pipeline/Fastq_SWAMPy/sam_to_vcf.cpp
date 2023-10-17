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
    while (std::getline(file, line)) {  // Load all lines into memory; ensure your file size is manageable for available RAM
        if (line.rfind('@', 0) == 0) continue; // Skip header
        lines.push_back(line);
    }

    std::map<std::tuple<std::string, int, std::string, std::string>, std::set<std::string>> variants_by_pos;
    std::set<std::string> reads;

    #pragma omp declare reduction (merge : std::set<std::string> : omp_out.insert(omp_in.begin(), omp_in.end()))
    #pragma omp declare reduction (map_merge : std::map<std::tuple<std::string, int, std::string, std::string>, std::set<std::string>> : omp_out.insert(omp_in.begin(), omp_in.end()))

    #pragma omp parallel for reduction(merge: reads) reduction(map_merge: variants_by_pos) schedule(dynamic)
    for (size_t i = 0; i < lines.size(); ++i) {
        auto [read_name, ref_name, pos, seq, mdz, cigar] = parse_sam_line(lines[i]);
        auto mismatches = parse_mdz(mdz, cigar, seq, pos, ref_genome);

        reads.insert(read_name);
        for (const auto& [position, ref_base, alt_base] : mismatches) {
            auto key = std::make_tuple(ref_name, position, ref_base, alt_base);
            variants_by_pos[key].insert(read_name);
        }
    }

    //std::ifstream file(sam_file);
    //std::string line;
    //std::mutex mutex;

    // #pragma omp parallel private(line) shared(file, variants_by_pos, reads, mutex)
    // {
    //     #pragma omp single nowait
    //     {
    //         while (std::getline(file, line)) {
    //             if (line.rfind('@', 0) == 0) continue; // Skip header

    //             #pragma omp task firstprivate(line)
    //             {
    //                 auto [read_name, ref_name, pos, seq, mdz, cigar] = parse_sam_line(line);
    //                 auto mismatches = parse_mdz(mdz, cigar, seq, pos, ref_genome);
                    
    //                 std::lock_guard<std::mutex> lock(mutex);
    //                 reads.insert(read_name);
    //                 for (const auto& [position, ref_base, alt_base] : mismatches) {
    //                     auto key = std::make_tuple(ref_name, position, ref_base, alt_base);
    //                     variants_by_pos[key].insert(read_name);
    //                 }
    //             }
    //         }
    //     }
    // }

    std::ofstream out(vcf_file);
    out << "##fileformat=VCFv4.2\n";
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    for (const auto& read : reads) {
        out << read << '\t';
    }
    out << '\n';

    for (const auto& [key, variant_reads] : variants_by_pos) {
        auto [ref_name, position, ref_base, alt_base] = key;
        out << ref_name << '\t' << position << '\t' << ref_base << position << alt_base << '\t' << ref_base << '\t' << alt_base << "\t.\t.\t.\t.\t";
        for (const auto& read : reads) {
            out << (variant_reads.find(read) != variant_reads.end() ? "1" : "0") << '\t';
        }
        out << '\n';
    }

    out.close();
    return 0;
}
