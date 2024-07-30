#pragma once

#include <string>
#include <vector>

#include "mutation.hpp"

struct raw_read {
   std::string read;
   std::vector<mutation> mutations;
   int start;
   int end;
   int degree;
};