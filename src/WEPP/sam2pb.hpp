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
       const int subsampled_reads;
   
       /* aligned/padded reads (possibly merged) */
       std::vector<sam_read> aligned_reads;
   
       /*
         these are the only the things that we need,
         since we can derive all other information from here
        */
   
       /* for any index, highlight all possible SNP */
       /* does NOT take into account read correction */
       sub_table frequency_table;
   
       /* frequency */
       /* may take into account read correction */
       sub_table collapsed_frequency_table;
   
       /* merging unmap information, column -> preimage of column */
       std::map<std::string, std::vector<std::string>> reverse_merge;
   
       void read_correction();
       void merge_duplicates();

       double divergence(std::vector<double> const& p, std::vector<double> const& q) {
           double divergence = 0;
           for (size_t j = 0; j < p.size(); ++j)
           {
               if (p[j] == 0)
                   continue;
               double effective_q = std::max(q[j], 1e-10);

               divergence += p[j] * (std::log(p[j]) - std::log(effective_q));
           }

           return divergence;
       }

       int get_binning_size();
       void subsample();

   public:
       sam(const std::string &ref, int subsampled_reads);

       void add_reads(const std::vector<std::string> &lines, tbb::blocked_range<size_t> range, tbb::queuing_mutex *mutex);
       void build();

       void dump_proto(std::string const& filename);
       void dump_reverse_merge(std::ostream& out);
};

void sam2PB(const dataset& d);
std::vector<raw_read> load_reads_from_proto(std::string const& reference, std::string const& filename, std::unordered_map<std::string, std::vector<std::string>> &reverse_merge);