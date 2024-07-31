#include "panman_bridge.hpp"

coord_converter::coord_converter(const panmanUtils::Tree& t) {
    const panmanUtils::BlockGapList &blockGaps = t.blockGaps;
    const std::vector<panmanUtils::Block> &blocks = t.blocks;
    const std::vector<panmanUtils::GapList> &gaps = t.gaps;

    map.resize(blocks.size() + 1);
    // Assigning block gaps
    for (size_t i = 0; i < blockGaps.blockPosition.size(); i++)
    {
        map[blockGaps.blockPosition[i]].second.resize(blockGaps.blockGapLength[i]);
    }

    std::string consensus = "";

    // map[i].first = nuc_pos
    // map[i].first

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
    int consensus_ctr = 0;
    for (size_t block_id = 0; block_id < map.size(); block_id++)
    {
        for (size_t j = 0; j < map[block_id].second.size(); j++)
        {
            for (size_t k = 0; k < map[block_id].second[j].size(); k++)
            {
                for (size_t w = 0; w < map[block_id].second[j][k].second.size(); w++)
                {
                    map[block_id].second[j][k].second[w] = ctr;
                    ctr++;
                    reference += '_';
                }
                map[block_id].second[j][k].first = ctr;
                ctr++;
                reference += '_';
            }
        }
        for (size_t j = 0; j < map[block_id].first.size(); j++)
        {
            for (size_t k = 0; k < map[block_id].first[j].second.size(); k++)
            {
                map[block_id].first[j].second[k] = ctr;
                ctr++;

                reference += '_';
            }
            map[block_id].first[j].first = ctr;
            ctr++;

            reference += consensus[consensus_ctr++];
        }
    }
}