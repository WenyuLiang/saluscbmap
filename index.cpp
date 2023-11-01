#include "index.h"

#include <assert.h>

#include "seed.h"
#include <algorithm>
#include <iostream>
#include <omp.h>

// static inline void printBinary(uint64_t num) {
//     for (uint64_t i = 63; i != static_cast<uint64_t>(-1); --i) {
//         if (i == 31) {
//             std::cout << " ";
//         }
//         std::cout << ((num & (1ULL << i)) ? '1' : '0');
//     }
//     std::cout << "\n";
// }

void Index::Construct(uint32_t num_sequences, const SequenceBatch &reference) {
  const double real_start_time = GetRealTime();

  std::vector<Seed> seeds;
  int num_seeds_per_sequence =
      ((reflen_ - 2 * offset_ - kmer_size_ + winsize_ - 1) / winsize_);
  seeds.reserve(reference.GetNumSequences() * num_seeds_per_sequence);
  std::cerr << "Collecting seeds.\n";

  // SeedGenerator seed_generator(offset_, kmer_size_, winsize_, length_);
  // for (uint32_t sequence_index = 0; sequence_index < num_sequences;
  //      ++sequence_index) {
  //   seed_generator.GenerateSeeds(reference, sequence_index, seeds);
  // }
  std::vector<std::vector<Seed>> seeds_per_thread(num_threads_);
  size_t estimated_size =
      (reference.GetNumSequences() * num_seeds_per_sequence) / num_threads_;

  for (auto &seed_vector : seeds_per_thread) {
    seed_vector.reserve(estimated_size);
  }
  // Thread-local instance of SeedGenerator
  SeedGenerator local_seed_generator(offset_, kmer_size_, winsize_, reflen_);

#pragma omp parallel for shared(seeds_per_thread)
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    int thread_id = omp_get_thread_num();
    local_seed_generator.GenerateSeeds(reference, sequence_index,
                                       seeds_per_thread[thread_id]);
  }

  // Merge seeds from all threads
  for (const auto &thread_seeds : seeds_per_thread) {
    seeds.insert(seeds.end(), thread_seeds.begin(), thread_seeds.end());
  }

  std::cerr << "Collected " << seeds.size() << " seeds.\n";
  std::cerr << "Sorting seeds.\n";
  std::stable_sort(seeds.begin(), seeds.end());
  std::cerr << "Sorted all seeds.\n";
  const size_t num_seeds = seeds.size();
  assert(num_seeds > 0);
  // TODO: check this assert!
  // Here I make sure the # seeds is less than the limit of signed int32,
  // so that I can use int to store position later.
  assert(num_seeds <= static_cast<size_t>(INT_MAX));
  occurrence_table_.reserve(num_seeds);
  kh_resize(k64, lookup_table_,
            num_seeds >> 8); // estimate the size of the hash table
  uint64_t previous_lookup_hash = GenerateHashInLookupTable(seeds[0].GetHash());
  uint32_t num_previous_seed_occurrences = 0;
  uint64_t num_nonsingletons = 0;
  uint32_t num_singletons = 0;
  for (size_t mi = 0; mi <= num_seeds; ++mi) {
    const bool is_last_iteration = mi == num_seeds;
    const uint64_t current_lookup_hash =
        is_last_iteration ? previous_lookup_hash + 1
                          : GenerateHashInLookupTable(seeds[mi].GetHash());

    // std::cout << "current_lookup_hash: " <<
    // (uint32_t)(seeds[mi].GetHit()>>32) << std::endl;
    if (current_lookup_hash != previous_lookup_hash) {
      int khash_return_code = 0;
      khiter_t khash_iterator =
          kh_put(k64, lookup_table_, previous_lookup_hash, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);

      // get hash value right here
      // kh_get_val(lookup_table_, khash_iterator) = seeds[mi].GetHit();

      if (num_previous_seed_occurrences == 1) {
        // We set the lowest bit of the key value to 1 if the minimizer only
        // occurs once. And the occurrence is directly saved in the lookup
        // table.
        kh_key(lookup_table_, khash_iterator) |= 1;
        kh_value(lookup_table_, khash_iterator) = occurrence_table_.back();
        // std::cout << "value1: " << (uint32_t)(kh_value(lookup_table_,
        // khash_iterator)>>32) << std::endl;
        occurrence_table_.pop_back();
        ++num_singletons;
      } else {
        kh_value(lookup_table_, khash_iterator) =
            GenerateEntryValueInLookupTable(num_nonsingletons,
                                            num_previous_seed_occurrences);
        num_nonsingletons += num_previous_seed_occurrences;
        // std::cout << "value2: " << (uint32_t)(kh_value(lookup_table_,
        // khash_iterator)>>32) << std::endl;
      }
      num_previous_seed_occurrences = 1;
    } else {
      num_previous_seed_occurrences++;
    }

    if (is_last_iteration) {
      break;
    }

    occurrence_table_.emplace_back(seeds[mi].GetHit());
    previous_lookup_hash = current_lookup_hash;
  }
  assert(num_nonsingletons + num_singletons == num_seeds);
  std::cerr << "Kmer size: " << kmer_size_ << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_)
            << ", # buckets: " << kh_n_buckets(lookup_table_)
            << ", occurrence table size: " << occurrence_table_.size()
            << ", # singletons: " << num_singletons << ".\n";
  std::cerr << "Built index successfully in " << GetRealTime() - real_start_time
            << "s.\n";

  //   for (auto k = kh_begin(h); k != kh_end(lookup_table_); ++k) {
  //     if (kh_exist(lookup_table_, k)) {
  //         //printf("key = %d, value = %d\n", kh_key(lookup_table_, k),
  //         kh_val(lookup_table_, k)); std::cout << "value3: " <<
  //         (uint32_t)(kh_value(lookup_table_, k)>>32) << std::endl;
  //     }
  // }
}

void Index::CheckIndex(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  std::vector<Seed> seeds;
  seeds.reserve(reference.GetNumSequences() *
                (reflen_ - 2 * offset_ - kmer_size_ + 1));
  SeedGenerator seed_generator(offset_, kmer_size_, winsize_, reflen_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    seed_generator.GenerateSeeds(reference, sequence_index, seeds);
  }
  std::cerr << "\nStart checking\n";
  std::cerr << "Collected " << seeds.size() << " seeds.\n";
  std::stable_sort(seeds.begin(), seeds.end());
  std::cerr << "Sorted seeds.\n";

  uint32_t count = 0;
  for (uint32_t i = 0; i < seeds.size(); ++i) {
    khiter_t khash_iterator = kh_get(
        k64, lookup_table_, GenerateHashInLookupTable(seeds[i].GetHash()));
    assert(khash_iterator != kh_end(lookup_table_));
    uint64_t key = kh_key(lookup_table_, khash_iterator);
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    if (key & 1) { // singleton
      assert(seeds[i].GetHit() == value);
      count = 0;
    } else {
      uint32_t offset = GenerateOffsetInOccurrenceTable(value);
      uint32_t num_occ = GenerateNumOccurrenceInOccurrenceTable(value);
      uint64_t value_in_index = occurrence_table_[offset + count];
      assert(value_in_index == seeds[i].GetHit());
      ++count;
      if (count == num_occ) {
        count = 0;
      }
    }
  }
}

int Index::GenerateCandidatePositions(MappingMetadata &mapping_metadata,
                                      SequenceBatch &ref, SequenceBatch &read,
                                      uint32_t &i) const {
  const size_t num_seeds = mapping_metadata.GetNumSeeds();
  const std::vector<Seed> &seeds = mapping_metadata.seed_; // seeds for a read
  for (size_t mi = 0; mi < num_seeds; ++mi) {
    khiter_t khash_iterator = kh_get(
        k64, lookup_table_, GenerateHashInLookupTable(seeds[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      continue;
    }
    const uint64_t read_hit = seeds[mi].GetHit();
    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    if (lookup_key & 1) {
      const uint64_t candidate_position = GenerateCandidatePositionFromHits(
          mapping_metadata, /*reference_hit=*/lookup_value, read_hit);
      if (candidate_position == UINT64_MAX) {
        continue;
      }
      mapping_metadata.positive_hits_.push_back(candidate_position);
    } else {
      // continue;
      uint32_t offset = GenerateOffsetInOccurrenceTable(lookup_value);
      uint32_t num_occ = GenerateNumOccurrenceInOccurrenceTable(lookup_value);

      for (uint32_t j = 0; j < num_occ; ++j) {
        const uint64_t candidate_position = GenerateCandidatePositionFromHits(
            mapping_metadata, /*reference_hit=*/occurrence_table_[offset + j],
            read_hit);
        if (candidate_position == UINT64_MAX) {
          continue;
        }
        mapping_metadata.positive_hits_.push_back(candidate_position);
      }
    }
    // std::cout << "reference_id: " << uint32_t(lookup_value>>32) << std::endl;
  }

  std::sort(mapping_metadata.positive_hits_.begin(),
            mapping_metadata.positive_hits_.end());

  // std::cout << "ref count: " << mapping_metadata.reference_hit_count_.size()
  // << std::endl; std::cout << "read seq: " << read.GetSequenceAt(i) <<
  // std::endl; for (auto const &pair: mapping_metadata.reference_hit_count_){
  //     // std::cout << "reference_id: " << ref.GetSequenceAt(pair.first) << "
  //     count: " << pair.second << "\n" << std::endl; std::cout <<
  //     "reference_id: " << pair.first  << " reference_seq: " <<
  //     ref.GetSequenceAt(pair.first)<< " count: " << pair.second << "\n" <<
  //     std::endl;
  // }

  // std::cout << "\n";
  return mapping_metadata.positive_hits_.size();
  // return mapping_metadata.reference_hit_count_.size();
}

uint64_t
Index::GenerateCandidatePositionFromHits(MappingMetadata &mapping_metadata,
                                         uint64_t reference_hit,
                                         uint64_t read_hit) const {

  const uint32_t reference_position = HitToSequencePosition(reference_hit);
  const uint32_t read_position = HitToSequencePosition(read_hit);
  // std::cout << "ref hit: " << printBinary(reference_hit) << std::endl;
  // printBinary(reference_hit);
  // For now we can't see the reference here. So let us don't validate this
  // candidate position. Instead, we do it later some time when we check the
  // candidates.
  const int reference_start_position = reference_position - read_position;
  const uint64_t reference_id = HitToSequenceIndex(reference_hit);

  if (reference_start_position > mapping_metadata.margin_ ||
      -reference_start_position > mapping_metadata.margin_) {
    // std::cout << "ref pos: " <<reference_position << " read pos: " <<
    // read_position<< " reference_start_position: " << reference_start_position
    // << " mapping_metadata.margin_ " << mapping_metadata.margin_<< std::endl;
    return UINT64_MAX;
  } // filter out the candidates whose reference start position is too far
  // away from the read start position

  mapping_metadata.reference_hit_count_[reference_id]++;
  return SequenceIndexAndPositionToCandidatePosition(reference_id,
                                                     reference_start_position);
}
