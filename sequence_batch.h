// sequence_batch.h
#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include "kseq.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <zlib.h>

class SequenceBatch {
public:
  KSEQ_INIT(gzFile, gzread);

  SequenceBatch() {}

  ~SequenceBatch() {
    if (sequence_batch_.size() > 0) {
      for (uint32_t i = 0; i < sequence_batch_.size(); ++i) {
        kseq_destroy(sequence_batch_[i]);
      }
    }
  }

  void InitializeLoading(const std::string &sequence_file_path);

  void FinalizeLoading();

  void LoadAllRefSequences();

  void LoadAllReadSequences();

  inline uint64_t GetNumBases() const { return num_bases_; }

  inline std::vector<kseq_t *> &GetSequenceBatch() { return sequence_batch_; }

  inline uint64_t GetNumSequences() const { return num_loaded_sequences_; }

  inline const char *GetSequenceAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.s;
  }

  inline const uint32_t GetSequenceLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.l;
  }

  inline const char *GetSequenceNameAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.s;
  }

  inline const char *GetSequenceCommentAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->comment.s;
  }

  inline uint32_t GetSequenceNameLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.l;
  }

  inline const char *GetSequenceQualAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->qual.s;
  }
  inline uint32_t GetSequenceIdAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->id;
  }

  inline void ModifySequenceAt(uint32_t sequence_index, const char *sequence,
                               uint32_t sequence_length) {
    sequence_batch_[sequence_index]->seq.s = strdup(sequence);
    sequence_batch_[sequence_index]->seq.l = sequence_length;
  }

  std::vector<bool> ref_sequence_keep_;
  std::vector<bool> ref_sequence_keep_modified_;

private:
  uint32_t total_num_loaded_sequences_ = 0;
  uint32_t num_loaded_sequences_ = 0;
  uint64_t num_bases_ = 0;
  gzFile sequence_file_;
  kseq_t *sequence_kseq_ = nullptr;
  std::vector<kseq_t *> sequence_batch_;
};

#endif // SEQUENCEBATCH_H_
