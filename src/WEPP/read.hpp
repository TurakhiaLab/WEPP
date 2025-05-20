#pragma once

#include <string>
#include <vector>

#include "src/usher_graph.hpp"

struct raw_read {
   std::string read;
   std::vector<MAT::Mutation> mutations;
   int start;
   int end;
   int degree;
};