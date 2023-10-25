#pragma once
#include "hit_utils.h"
#include <cassert>
#include <cstdint>
#include <vector>

#include "sequence_batch.h"

class Seed {
public:
  Seed() = delete;

  Seed(std::pair<uint64_t, uint64_t> seed)
      : hash_(seed.first), hit_(seed.second) {}

  Seed(uint64_t hash, uint64_t hit) : hash_(hash), hit_(hit) {}

  ~Seed() = default;

  inline uint64_t GetHash() const { return hash_; }

  inline uint64_t GetHit() const { return hit_; }

  inline uint32_t GetSequenceIndex() const { return HitToSequenceIndex(hit_); }

  inline uint32_t GetSequencePosition() const {
    return HitToSequencePosition(hit_);
  }

  inline bool operator<(const Seed &m) const {
    if (hash_ < m.hash_) {
      return true;
    }

    if (hash_ == m.hash_ && hit_ < m.hit_) {
      return true;
    }

    return false;
  }

private:
  // The hash of the kmer.
  uint64_t hash_ = 0;

  // The high 32 bits save the sequence index in the sequence batch. The
  // following 32 bits save the end position on that sequence.
  uint64_t hit_ = 0;
};

class SeedGenerator {
public:
  SeedGenerator() = delete;
  // offset: offset from both ends of the sequence
  // kmer_size_: size of the k-mer
  // length: length of the sequence
  SeedGenerator(int offset, int kmer_size_, int winsize, int length)
      : offset_(offset), kmer_size_(kmer_size_), winsize_(winsize), length_(length) {
    // 56 bits for a k-mer. So the max kmer size is 28.
    assert(kmer_size_ > 0 && kmer_size_ <= 28);
  }

  ~SeedGenerator() = default;

  void GenerateSeeds(const SequenceBatch &sequence_batch,
                     uint32_t sequence_index, std::vector<Seed> &Seeds) const;

private:
  const int offset_;
  const int kmer_size_;
  const int winsize_;
  const int length_;
};