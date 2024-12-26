#pragma once

#include <string>
#include <vector>
#include <array>
#include <map>

#include "dataset.hpp"

const std::string GENOME_STRING{"ACGTN_"};
typedef std::vector<std::array<int, 6>> sub_table;

struct sam_read {
    std::string raw_name;
    int start_idx;
    int degree; // how many repeats;
    std::string aligned_string;

    std::string degree_name() const {
        return raw_name + 
            "_READ_" + std::to_string(start_idx + 1) + 
            "_" + std::to_string(start_idx + 1 + aligned_string.size() - 1) + 
            "_" + std::to_string(degree);
    }

    std::string degreeless_name() const {
        return raw_name + 
            "_READ_" + std::to_string(start_idx + 1) + 
            "_" + std::to_string(start_idx + 1 + aligned_string.size() - 1); 
    }

    bool operator <(const sam_read& rhs) const {
        if (start_idx != rhs.start_idx) {
            return start_idx < rhs.start_idx;
        }

        if (aligned_string.size() != rhs.aligned_string.size()) {
            return aligned_string.size() < rhs.aligned_string.size();
        }

        return aligned_string < rhs.aligned_string;
    }

    bool operator ==(const sam_read& rhs) const {
        return start_idx == rhs.start_idx && aligned_string == rhs.aligned_string;
    }
};

struct sam {
private:
    const std::string reference_seq;

    /* aligned/padded reads (possibly merged) */
    std::vector<sam_read> aligned_reads;

    /*
      these are the only the things that we need,
      since we can derive all other information from here
     */

    /* for any index, highlight all possible SNP */
    /* does NOT take into account read correction */
    sub_table frequency_table;

    /* frequency + collapsed indels to reference */
    /* may take into account read correction */
    sub_table collapsed_frequency_table;

    /* for any starting index, length possible inserions (length in reference, replacement string)  */
    std::vector<std::map<std::pair<size_t, std::string>, int>> insertion_frequency_table;

    /* merging unmap information, column -> preimage of column */
    std::map<std::string, std::vector<std::string>> reverse_merge;

    void read_correction();
    void merge_duplicates();

public:
    sam(const std::string &ref);
    sam(const std::string &ref, const std::vector<sam_read> &raw_reads);

    void add_read(const std::string &line);
    void build();

    void dump_fake_sam(std::string const &filename);
    void dump_proto(std::string const &filename);
    void dump_reverse_merge(std::ostream &out);
};

void sam2PB(const dataset& d);
std::vector<raw_read> load_reads_from_proto(std::string const& reference, std::string const& filename, std::unordered_map<std::string, std::vector<std::string>> &reverse_merge);