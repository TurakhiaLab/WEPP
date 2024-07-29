#pragma once

#include <string>
#include <vector>

#include "mutation.hpp"

struct raw_read {
   std::string read;
   std::vector<mutation> mutations;
   size_t start;
   size_t end;
   int degree;
};