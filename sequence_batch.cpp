// batch.cpp
#include "sequence_batch.h"
#include "utils.h"
#include <algorithm>
#include <cstring>
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

void SequenceBatch::LoadAllRefSequences() {
  double real_start_time = GetRealTime();
  sequence_batch_.reserve(90000000);
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
      if (sequence_kseq_->qual.l != 0) [[unlikely]] { // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
      }

      if (sequence->comment.l == 0) {
        std::string name(sequence->name.s);
        // replace name with id to save memory
        // Here, replace the name with "0"
        free(sequence->name.s);
        sequence->name.l = sequence->name.m = 1;
        sequence->name.s = strdup("0");
        // replace name with id to save memory

        size_t lastUnderscore = std::string::npos;
        size_t secondLastUnderscore = std::string::npos;

        // Find the positions of the last two underscores in a single pass
        for (size_t i = name.length(); i-- > 0;) {
          if (name[i] == '_') {
            if (lastUnderscore == std::string::npos) {
              lastUnderscore = i;
            } else {
              secondLastUnderscore = i;
              break;
            }
          }
        }

        if (secondLastUnderscore != std::string::npos) {
          std::string coordinates =
              name.substr(secondLastUnderscore + 1, name.length());
          // Replace underscores with spaces
          std::replace(coordinates.begin(), coordinates.end(), '_', ' ');
          // Assign the extracted coordinates to the comment
          sequence->comment.s = strdup(coordinates.c_str());
          sequence->comment.l = coordinates.length();
        }
      }
      sequence->id = total_num_loaded_sequences_;
      ++total_num_loaded_sequences_;
      ++num_loaded_sequences_;
      num_bases_ += length;
    }
    length = kseq_read(sequence_kseq_);
  }
  sequence_batch_.shrink_to_fit();
  ref_sequence_keep_.resize(sequence_batch_.size(), false);
  ref_sequence_keep_modified_.resize(sequence_batch_.size(), false);

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

bool SequenceBatch::LoadBatchReadSequences() {
  if (sequence_batch_.size() > 0) {
    for (uint32_t i = 0; i < sequence_batch_.size(); ++i) {
      kseq_destroy(sequence_batch_[i]);
    }
  }
  sequence_batch_.clear();
  num_loaded_sequences_ = 0;
  while (num_loaded_sequences_ < batch_size_ && kseq_read(sequence_kseq_) > 0) {
    sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
    kseq_t *sequence = sequence_batch_.back();
    std::swap(sequence_kseq_->seq, sequence->seq);
    std::swap(sequence_kseq_->name, sequence->name);
    std::swap(sequence_kseq_->comment, sequence->comment);
    if (sequence_kseq_->qual.l != 0) { // fastq file
      std::swap(sequence_kseq_->qual, sequence->qual);
    }
    sequence->id = num_loaded_sequences_;
    ++num_loaded_sequences_;
  }
  if (sequence_batch_.size() < batch_size_) {
    sequence_batch_.shrink_to_fit();
    return true;
  }
  return false; // batch full, meaning that there are more reads to be loaded
}

// void SequenceBatch::LoadAllReadSequences() {
//   double real_start_time = GetRealTime();
//   sequence_batch_.reserve(30000000);
//   num_loaded_sequences_ = 0;
//   num_bases_ = 0;
//   int length = kseq_read(sequence_kseq_);
//   while (length >= 0) {
//     if (length > 0) {
//       sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
//       kseq_t *sequence = sequence_batch_.back();
//       std::swap(sequence_kseq_->seq, sequence->seq);
//       std::swap(sequence_kseq_->name, sequence->name);
//       std::swap(sequence_kseq_->comment, sequence->comment);
//       if (sequence_kseq_->qual.l != 0) { // fastq file
//         std::swap(sequence_kseq_->qual, sequence->qual);
//       }
//       sequence->id = total_num_loaded_sequences_;
//       ++total_num_loaded_sequences_;
//       ++num_loaded_sequences_;
//       num_bases_ += length;
//     }
//     length = kseq_read(sequence_kseq_);
//   }
//   sequence_batch_.shrink_to_fit();
//   // Make sure to reach the end of the file rather than meet an error.
//   if (length != -1) {
//     std::cerr
//         << "Didn't reach the end of sequence file, which might be
//         corrupted!\n";
//   }
//   std::cerr << "Loaded all sequences successfully in "
//             << GetRealTime() - real_start_time << "s, ";
//   std::cerr << "number of sequences: " << num_loaded_sequences_ << ", ";
//   std::cerr << "number of bases: " << num_bases_ << ".\n";
// }