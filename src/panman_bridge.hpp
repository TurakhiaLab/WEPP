#pragma once

#include <vector>
#include <utility>

#include "panman/panmanUtils.hpp"

// https://github.com/amkram/panmap/blob/main/src/tree.hpp#L37
typedef std::vector<std::pair<std::vector<std::pair<int,std::vector<int>>>, std::vector<std::vector<std::pair<int,std::vector<int>>>>>> globalCoords_t;

class coord_converter {
    globalCoords_t map;
public:
    coord_converter(panmanUtils::Tree& t) {
        const panmanUtils::BlockGapList & blockGaps = t.blockGaps;
        const std::vector<panmanUtils::Block>& blocks = t.blocks;
        const std::vector<panmanUtils::GapList>& gaps = t.gaps;

        map.resize(blocks.size() + 1);
        // Assigning block gaps
        for (size_t i = 0; i < blockGaps.blockPosition.size(); i++)
        {
            map[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
        }
        int32_t maxBlockId = 0;
        for (size_t i = 0; i < blocks.size(); i++)
        {
            int32_t blockId = ((int32_t)blocks[i].primaryBlockId);
            maxBlockId = std::max(maxBlockId, blockId);
            for (size_t j = 0; j < blocks[i].consensusSeq.size(); j++)
            {
                bool endFlag = false;
                for (size_t k = 0; k < 8; k++)
                {
                    const int nucCode = (((blocks[i].consensusSeq[j]) >> (4 * (7 - k))) & 15);
                    if (nucCode == 0)
                    {
                        endFlag = true;
                        break;
                    }
                    map[blockId].first.push_back({0, {}});
                }
                if (endFlag)
                {
                    break;
                }
            }
            map[blockId].first.push_back({0, {}});
        }
        map.resize(maxBlockId + 1);
        // Assigning nucleotide gaps
        for (size_t i = 0; i < gaps.size(); i++)
        {
            int32_t blockId = (gaps[i].primaryBlockId);
            for (size_t j = 0; j < gaps[i].nucPosition.size(); j++)
            {
                int len = gaps[i].nucGapLength[j];
                int pos = gaps[i].nucPosition[j];
                map[blockId].first[pos].second.resize(len, 0);
            }
        }
        // Assigning coordinates
        int ctr = 1;
        for (size_t i = 0; i < map.size(); i++)
        {
            for (size_t j = 0; j < map[i].second.size(); j++)
            {
                for (size_t k = 0; k < map[i].second[j].size(); k++)
                {
                    for (size_t w = 0; w < map[i].second[j][k].second.size(); w++)
                    {
                        map[i].second[j][k].second[w] = ctr;
                        ctr++;
                    }
                    map[i].second[j][k].first = ctr;
                    ctr++;
                }
            }
            for (size_t j = 0; j < map[i].first.size(); j++)
            {
                for (size_t k = 0; k < map[i].first[j].second.size(); k++)
                {
                    map[i].first[j].second[k] = ctr;
                    ctr++;
                }
                map[i].first[j].first = ctr;
                ctr++;
            }
        }
    }

    size_t query(const int blockId, const int nucPosition, const int nucGapPosition) {
        if(nucGapPosition == -1){
            return map[blockId].first[nucPosition].first;
        }
        return map[blockId].first[nucPosition].second[nucGapPosition];
    }
};
