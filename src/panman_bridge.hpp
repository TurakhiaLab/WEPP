#pragma once

#include <vector>
#include <utility>

#include "panman/panmanUtils.hpp"

// https://github.com/amkram/panmap/blob/main/src/tree.hpp#L37
typedef std::vector<std::pair<std::vector<std::pair<int,std::vector<int>>>, std::vector<std::vector<std::pair<int,std::vector<int>>>>>> globalCoords_t;

class coord_converter {
    globalCoords_t map;
public:
    std::string reference;

    coord_converter(const panmanUtils::Tree& t);

    size_t query(const int blockId, const int nucPosition, const int nucGapPosition) const {
        if(nucGapPosition == -1){
            return map[blockId].first[nucPosition].first;
        }
        return map[blockId].first[nucPosition].second[nucGapPosition];
    }
};
