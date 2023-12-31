#include "seed.h"
#include "hit_utils.h"
#include "sequence_batch.h"
#include <algorithm>
#include <cstring>
#include <string>

void SeedGenerator::GenerateSeeds(const SequenceBatch &sequence_batch,
                                  uint32_t sequence_index,
                                  std::vector<Seed> &Seeds) const {
  const uint32_t sequence_length =
      sequence_batch.GetSequenceLengthAt(sequence_index);
  std::pair<uint64_t, uint64_t> seed_pair = {UINT64_MAX, UINT64_MAX};
  // If the sequence is shorter than the k-mer size, return early.
  if (sequence_length < kmer_size_) {
    return;
  }
  const uint64_t mask = (((uint64_t)1) << (2 * kmer_size_)) - 1;
  const char *sequence = sequence_batch.GetSequenceAt(sequence_index);
  // Define the start and end positions for k-mer generation.
  int start = offset_;
  int end = sequence_length - offset_ - kmer_size_ + 1;
  uint64_t seed = 0;

  // -------------------------------------------------check
  // seed------------------------------------------------------------------
  // std::string seed_str = "";
  // seed_str.reserve(kmer_size_);
  // -------------------------------------------------check
  // seed------------------------------------------------------------------
  // Initialize the first k-mer.
  for (uint32_t i = start; i < start + kmer_size_; ++i) {
    const uint8_t current_base = CharToUint8(sequence[i]);
    if (current_base < 4) [[likely]] {
      seed = (seed << 2) | current_base;
    }
    // seed_str+= sequence[i];
  }
  // std::cout << "seed: " << seed_str << "\n";
  seed_pair.first = Hash64(seed, mask); // hash of the first k-mer
  seed_pair.second = (((uint64_t)sequence_index) << 32 | (uint32_t)start);
  // The high 31 bits save the sequence index in the sequence batch. The
  // following 32 bits save the end position on that sequence. And the lowest
  // bit encodes the strand (0 for positive).
  if (seed_pair.first != UINT64_MAX) {
    Seeds.emplace_back(seed_pair);
  } else
    std::cerr << "Fail to add seed\n";

  uint64_t num_shifts = generateNumber(winsize_);
  // Slide the window and compute subsequent k-mers.
  // for (uint32_t position = start + 1; position < end; ++position) {
  for (uint32_t position = start + winsize_; position < end;
       position += winsize_) {
    // Remove the leftmost base of the previous k-mer.
    // seed &= ~(((uint64_t)3) << (2 * (kmer_size_ - 1)));
    seed &= ~(((uint64_t)num_shifts) << (2 * (kmer_size_ - winsize_)));

    // seed = (seed << 2) | current_base;

    for (int i = 0; i < winsize_; i++) {
      // Add the new base to the current k-mer.
      const uint8_t current_base =
          CharToUint8(sequence[position + (kmer_size_ - winsize_) + i]);
      if (current_base < 4) [[likely]] {
        seed = (seed << 2) | current_base;
      } else {
        seed = (seed << 2) | 0; // if the base is N, set it to A
        // for (int j = 0; j < 4; j++) {
        //   seed = (seed << 2) | j;
        //   seed_pair.first = Hash64(seed, mask);
        //   seed_pair.second =
        //       (((uint64_t)sequence_index) << 32 | (uint32_t)position);
        //   if (seed_pair.first != UINT64_MAX) {
        //     Seeds.emplace_back(seed_pair);
        //   } else
        //     std::cerr << "Fail to add seed\n";
        // }
      }
    }
    seed_pair.first = Hash64(seed, mask);
    seed_pair.second = (((uint64_t)sequence_index) << 32 | (uint32_t)position);

    // -------------------------------------------------check
    // seed_str.clear();
    // for(size_t i=0; i<kmer_size_*2; i+=2) {
    //     seed_str += Uint8ToChar((seed >> (i)) & 3);
    // }
    // std::reverse(seed_str.begin(), seed_str.end());
    // std::cout << "seed: " << seed_str << "\n";
    // -------------------------------------------------check

    if (seed_pair.first != UINT64_MAX) {
      Seeds.emplace_back(seed_pair);
    } else
      std::cerr << "Fail to add seed\n";
  }
}
