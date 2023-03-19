#include "amplicon.hpp"
#include "text_parser.hpp"
#include "tbb/concurrent_unordered_set.h"
#include <memory>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

// Get starting and ending genomic coordinates of where the amplicon aligned to
// reference
void get_coordinates(
    std::string fasta_file_path,
    std::vector<std::tuple<int, int>> &samples_start_end_coordinates) {

    // Read in entire fasta file
    text_parser fasta(fasta_file_path);

    // Skip over reference sequence
    fasta.next_line();

    char at = '>';
    for (; !fasta.done(); fasta.next_line()) {
        if (fasta.get_value(0)[0] != at) {
            continue;
        }
        // Get start and end coordinates for each sample from fasta file.
        // Each sequence in amplicons.fa starts with tab delimited header line >
        // Column 1 = amplicon sample name
        // Column 2 = amplicon size
        // Column 3 = start coordinate from alignment
        // Column 4 = end coordinate from alignment
        // eg) >sample_name\t151\t659\t809

        // Create tuple <start,end> coordinates, used for masking in placement
        int start = str_view_to_uint64(fasta.get_value(2));
        int end = str_view_to_uint64(fasta.get_value(3));
        std::tuple<int, int> tup = std::make_tuple(start, end);
        samples_start_end_coordinates.push_back(tup);
    }
}

inline uint64_t str_view_to_uint64(std::string_view str) noexcept {
    uint64_t num = 0;
    for (auto &c : str)
        (num *= 10) += c - '0';
    return num;
}
