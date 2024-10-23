#include "panman_bridge.hpp"

coord_converter::coord_converter(const panmanUtils::Tree& t) {
    const std::vector<panmanUtils::Block> &blocks = t.blocks;
    const std::vector<panmanUtils::GapList> &gaps = t.gaps;

    std::string consensus = "";
    map.resize(blocks.size() + 1);
    int32_t maxBlockId = 0;
    for (size_t block_id = 0; block_id < blocks.size(); block_id++)
    {
        int32_t blockId = ((int32_t)blocks[block_id].primaryBlockId);
        maxBlockId = std::max(maxBlockId, blockId);
        for (size_t nuc_pos = 0; nuc_pos < blocks[block_id].consensusSeq.size(); nuc_pos++)
        {
            bool endFlag = false;
            for (size_t k = 0; k < 8; k++)
            {
                const int nucCode = (((blocks[block_id].consensusSeq[nuc_pos]) >> (4 * (7 - k))) & 15);
                if (nucCode == 0)
                {
                    endFlag = true;
                    break;
                }
                consensus += panmanUtils::getNucleotideFromCode(nucCode);
                map[blockId].push_back({0, {}});
            }
            if (endFlag)
            {
                break;
            }
        }
        // Added to handle pangraph based panmat format
        //map[blockId].push_back({0, {}});
        //consensus += '_';
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
            map[blockId][pos].second.resize(len, 0);
        }
    }
    
    // Assigning coordinates
    int ctr = 1;
    int consensus_ctr = 0;
    for (size_t block_id = 0; block_id < map.size(); block_id++)
    {
        for (size_t j = 0; j < map[block_id].size(); j++)
        {
            for (size_t k = 0; k < map[block_id][j].second.size(); k++)
            {
                map[block_id][j].second[k] = ctr;
                ctr++;

                reference += '_';
            }
            map[block_id][j].first = ctr;
            ctr++;

            reference += consensus[consensus_ctr++];
        }
    }
}