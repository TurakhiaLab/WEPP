#ifndef AMPLICON_H
#define AMPLICON_H

#include "src/mutation_annotated_tree.hpp"
#include "text_parser.hpp"
#include "src/usher_graph.hpp"

namespace MAT = Mutation_Annotated_Tree;

void get_coordinates(
    std::string input_fasta_filename,
    std::vector<std::tuple<int, int>> &samples_start_end_coordinates);

inline uint64_t str_view_to_uint64(std::string_view str) noexcept;


#endif

