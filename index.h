
#pragma once
#include <limits>
#include <queue>
#include <string>
#include <vector>
#include "khash.h"
#include "sequence_batch.h"
#include "mapping_metadata.h"
#include "index_utils.h" 
#define reflen 30

struct IndexParameters {
  int kmer_size = 9;
  int offset = 4;
  int num_threads = 1;
  std::string reference_file_path;
  std::string index_output_file_path;
};

class Index {
 public:
  Index() = delete;

  // For read mapping.
  Index(const std::string &index_file_path)
      : index_file_path_(index_file_path) {
    lookup_table_ = kh_init(k64);
  }

  // For index construction.
  Index(const IndexParameters &index_parameters)
      : 
        offset_(index_parameters.offset),
        kmer_size_(index_parameters.kmer_size),
        num_threads_(index_parameters.num_threads),
        index_file_path_(index_parameters.index_output_file_path) {
    lookup_table_ = kh_init(k64);
  }

  ~Index() { Destroy(); }

  void Destroy() {
    if (lookup_table_ != nullptr) {
      kh_destroy(k64, lookup_table_);
      lookup_table_ = nullptr;
    }

    std::vector<uint64_t>().swap(occurrence_table_);
  }

  void Construct(uint32_t num_sequences, const SequenceBatch &reference);
  void CheckIndex(uint32_t, const SequenceBatch&) const;
  // // Output index stats.
  // void Statistics(uint32_t num_sequences, const SequenceBatch &reference) const;

  // // Check the index for some reference genome. Only for debug.
  // void CheckIndex(uint32_t num_sequences, const SequenceBatch &reference) const;

  // Return the number of repetitive seeds.
  int GenerateCandidatePositions(MappingMetadata &mapping_metadata, SequenceBatch &ref, SequenceBatch &read, uint32_t &i) const;

  // Return the number of matched seeds on candidate reference id.
  // unordered_map<uint32_t, uint32_t> CandidateStats(MappingMetadata &mapping_metadata, uint32_t candidate_reference_id);

  int GetKmerSize() const { return kmer_size_; }

  uint32_t GetLookupTableSize() const { return kh_size(lookup_table_); }

 private:
  uint64_t GenerateCandidatePositionFromHits(MappingMetadata &mapping_metadata, uint64_t reference_hit,
                                             uint64_t read_hit) const;

  int kmer_size_ = 9;
  int offset_ = 4;
  int length_ = 30;
  // Number of threads to build the index, which is not used right now.
  int num_threads_ = 1;
  const std::string index_file_path_;
  khash_t(k64) *lookup_table_ = nullptr;
  
  std::vector<uint64_t> occurrence_table_;
};