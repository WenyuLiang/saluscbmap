#ifndef MAPPING_METADATA_H_
#define MAPPING_METADATA_H_

#include "seed.h"
#include <algorithm>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>
struct Candidate {
  // The high 32 bits save the reference sequence index in the reference
  // sequence batch. The low 32 bits save the reference position on that
  // sequence.
  uint64_t position = 0;

  // The number of minimizers supports the position.
  uint8_t count = 0;

  inline uint32_t GetReferenceSequenceIndex() const { return (position >> 32); }

  inline uint32_t GetReferenceSequencePosition() const { return position; }

  inline bool operator<(const Candidate &c) const {
    if (count > c.count) {
      return true;
    }

    if (count < c.count) {
      return false;
    }

    return position < c.position;
  }
};

class MappingMetadata {
public:
  explicit MappingMetadata(int margin) : margin_(margin){};

  inline void PrepareForMappingNextRead(int reserve_size) {
    Clear();
    ReserveSpace(reserve_size);
  }

  inline size_t GetNumSeeds() const { return seed_.size(); }

  inline size_t GetNumPositiveCandidates() const {
    return positive_candidates_.size();
  }

  inline void MoveCandidiatesToBuffer() {
    positive_candidates_.swap(positive_candidates_buffer_);
    positive_candidates_.clear();
  }

  inline void SortCandidates() {
    std::sort(positive_candidates_.begin(), positive_candidates_.end());
  }

  inline int GetMinNumErrors() const { return min_num_errors_; }
  inline int GetSecondMinNumErrors() const { return second_min_num_errors_; }
  inline int GetNumBestMappings() const { return num_best_mappings_; }
  inline int GetNumSecondBestMappings() const {
    return num_second_best_mappings_;
  }

  inline void SetMinNumErrors(int min_num_errors) {
    min_num_errors_ = min_num_errors;
  }
  inline void SetSecondMinNumErrors(int second_min_num_errors) {
    second_min_num_errors_ = second_min_num_errors;
  }
  inline void SetNumBestMappings(int num_best_mappings) {
    num_best_mappings_ = num_best_mappings;
  }
  inline void SetNumSecondBestMappings(int num_second_best_mappings) {
    num_second_best_mappings_ = num_second_best_mappings;
  }
  inline int GetNumRefHits() const { return reference_hit_count_.size(); }

  // static std::vector<int> reference_hit_count_;

  // protected:
  // inline void InitReferenceHitCount(uint32_t num_sequences) {
  //    reference_hit_count_ = std::vector<int>(num_sequences, 0);
  // }

  //  protected:
  inline void ReserveSpace(int reserve_size) {
    seed_.reserve(reserve_size);
    positive_hits_.reserve(reserve_size);
    positive_candidates_.reserve(reserve_size);
    positive_candidates_buffer_.reserve(reserve_size);
  }

  inline void Clear() {
    seed_.clear();
    positive_hits_.clear();
    positive_candidates_.clear();
    positive_candidates_buffer_.clear();
    reference_hit_count_.clear();
  }

  int min_num_errors_, second_min_num_errors_;

  int num_best_mappings_, num_second_best_mappings_;

  const uint32_t margin_;

  std::vector<Seed> seed_;

  std::vector<uint64_t> positive_hits_;

  std::vector<Candidate> positive_candidates_;

  std::vector<Candidate> positive_candidates_buffer_;

  std::map<uint32_t, int> reference_hit_count_;
};

#endif // MAPPING_METADATA_H_