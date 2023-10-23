#pragma once
#include <utility> // for std::pair
#include <queue>   // for std::priority_queue
#include <vector>  // for std::vectorc

#include "mapping_metadata.h"
#include "sequence_batch.h"
class Align{
    public:
    static inline bool cmp(const std::pair<uint32_t, int>& left, const std::pair<uint32_t, int>& right);
    int CandidateRef(const MappingMetadata &mapping_metadata, const SequenceBatch &ref, SequenceBatch &read, uint32_t &i);
    // private:
    size_t valid_mapping_count_ = 0;
    size_t invalid_mapping_count_ = 0;
    std::vector<std::string> unmapped_reads_;
};