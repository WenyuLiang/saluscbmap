#pragma once
#include <queue>   // for std::priority_queue
#include <utility> // for std::pair
#include <vector>  // for std::vectorc

#include "mapping_metadata.h"
#include "sequence_batch.h"

using TripleType = std::tuple<
    uint32_t, 
    uint32_t,
    int
>;

class Align {
public:
  explicit Align(uint32_t edit_distance) : edit_distance_(edit_distance){};
  static inline bool cmp(const std::pair<uint32_t, int> &left,
                         const std::pair<uint32_t, int> &right);
  int CandidateRef(const MappingMetadata &mapping_metadata,
                   SequenceBatch &ref, SequenceBatch &read, uint32_t &i);
  // private:
  size_t valid_mapping_count_ = 0;
  size_t invalid_mapping_count_ = 0;
  const uint32_t edit_distance_ = 0;
  //std::vector<std::string> unmapped_reads_;
};