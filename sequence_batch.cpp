// batch.cpp
#include "sequence_batch.h"
#include "utils.h"
#include <tuple>

void SequenceBatch::InitializeLoading(const std::string &sequence_file_path) {
  sequence_file_ = gzopen(sequence_file_path.c_str(), "r");
  if (sequence_file_ == NULL) {
    std::cerr << "Cannot find sequence file " << sequence_file_path << ".\n";
  }
  sequence_kseq_ = kseq_init(sequence_file_);
}

void SequenceBatch::FinalizeLoading() {
  kseq_destroy(sequence_kseq_);
  gzclose(sequence_file_);
}

void SequenceBatch::LoadAllSequences() {
  double real_start_time = GetRealTime();
  sequence_batch_.reserve(30000000);
  num_loaded_sequences_ = 0;
  num_bases_ = 0;
  int length = kseq_read(sequence_kseq_);
  while (length >= 0) {
    if (length > 0) {
      sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
      kseq_t *sequence = sequence_batch_.back();
      std::swap(sequence_kseq_->seq, sequence->seq);
      std::swap(sequence_kseq_->name, sequence->name);     
      std::swap(sequence_kseq_->comment, sequence->comment);
      if (sequence_kseq_->qual.l != 0) { // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
      }
      sequence->id = total_num_loaded_sequences_;

      ++total_num_loaded_sequences_;
      ++num_loaded_sequences_;
      num_bases_ += length;
    }
    length = kseq_read(sequence_kseq_);
  }

  // Make sure to reach the end of the file rather than meet an error.
  if (length != -1) {
    std::cerr
        << "Didn't reach the end of sequence file, which might be corrupted!\n";
  }

  std::cerr << "Loaded all sequences successfully in "
            << GetRealTime() - real_start_time << "s, ";
  std::cerr << "number of sequences: " << num_loaded_sequences_ << ", ";
  std::cerr << "number of bases: " << num_bases_ << ".\n";
}
