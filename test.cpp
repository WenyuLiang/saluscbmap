#include <cstdint>
#include <iostream>
// #define WIN_SIZE 4
// void SeedGenerator::GenerateSeeds(const SequenceBatch &sequence_batch,
//                                   uint32_t sequence_index,
//                                   std::vector<Seed> &Seeds) const {
//     const uint32_t sequence_length =
//     sequence_batch.GetSequenceLengthAt(sequence_index); std::pair<uint64_t,
//     uint64_t> seed_pair = {UINT64_MAX, UINT64_MAX};
//     // If the sequence is shorter than the k-mer size, return early.
//     if (sequence_length < kmer_size_) {
//         return;
//     }
//     const uint64_t mask = (((uint64_t)1) << (2 * kmer_size_)) - 1;
//     const char *sequence = sequence_batch.GetSequenceAt(sequence_index);

//     // Define the start and end positions for k-mer generation.
//     int start = offset_;
//     int end = sequence_length - offset_ - kmer_size_ + 1;

//     uint64_t seed = 0;

//     // -------------------------------------------------check
//     seed------------------------------------------------------------------
//     std::string seed_str = "";
//     seed_str.reserve(kmer_size_);
//     // -------------------------------------------------check
//     seed------------------------------------------------------------------
//     // Initialize the first k-mer.
//     for (uint32_t i = start; i < start + kmer_size_; ++i) {
//         const uint8_t current_base = CharToUint8(sequence[i]);
//         if (current_base < 4) {
//             seed = (seed << 2) | current_base;
//         }
//         seed_str+= sequence[i];
//     }
//     if(seed_str == "GCACTCAAATCACGT")
//     {std::cout <<"seq: " << sequence << "\n";
//     std::cout <<"seed exist1: " << seed_str << std::endl;}
//     if(strcmp(sequence, "CAACTAAGCACTCAAATCACGTGCGCTATA")==0) std::cout <<
//     "seedYES1: " << seed_str << std::endl; seed_pair.first = Hash64(seed,
//     mask); // hash of the first k-mer
//     // std::cout << sequence_index << std::endl;
//     // seed_pair.second = ((((uint64_t)sequence_index) << 32 |
//     (uint32_t)start) << 1); seed_pair.second = (((uint64_t)sequence_index) <<
//     32 | (uint32_t)start);
//     //std::cout << "ref id: " << HitToSequenceIndex(seed_pair.second) <<
//     "\tref pos: " << HitToSequencePosition(seed_pair.second) <<"\n"<<
//     std::endl;

//     // The high 31 bits save the sequence index in the sequence batch. The
//     // following 32 bits save the end position on that sequence. And the
//     lowest
//     // bit encodes the strand (0 for positive).
//     if (seed_pair.first != UINT64_MAX) {
//             Seeds.emplace_back(seed_pair);
//         } else std::cerr << "fail!!!!!!!!!\n";

//     // Slide the window and compute subsequent k-mers.
//     for (uint32_t position = start + 1; position < end; ++position) {
//         seed &= ~(((uint64_t)3) << (2 * (kmer_size_ - 1)));

//         // Add the new base to the current k-mer.
//         const uint8_t current_base = CharToUint8(sequence[position +
//         kmer_size_ - 1]); seed = (seed << 2) | current_base;
//         seed_str.clear();
//         for(size_t i=0; i<kmer_size_*2; i+=2) {
//             seed_str += Uint8ToChar((seed >> (i)) & 3);
//         }
//         std::reverse(seed_str.begin(), seed_str.end());
//         if(seed_str == "GCACTCAAATCACGT")
//         std::cout << "seed exist2: " << seed_str << "\n";
//         if(strcmp(sequence, "CAACTAAGCACTCAAATCACGTGCGCTATA")==0) std::cout
//         << "seedYES2: " << seed_str << std::endl; if (current_base < 4) {
//             seed = (seed << 2) | current_base;
//         }

//         seed_pair.first = Hash64(seed, mask);
//         seed_pair.second = (((uint64_t)sequence_index) << 32 |
//         (uint32_t)position); if (seed_pair.first != UINT64_MAX) {
//             Seeds.emplace_back(seed_pair);
//         } else std::cerr << "fail!!!!!!!!!\n";

//         // std::cout << "ref id: " << HitToSequenceIndex(seed_pair.second) <<
//         "\tref pos: " << HitToSequencePosition(seed_pair.second) <<"\n"<<
//         std::endl;
//     }
// std::cout<<"\n";

// constexpr uint64_t generateNumber(int n) {
//     return (1 << n*2) - 1;
// }

typedef struct {
  std::string input;
  std::string output;
  uint32_t kmer;
  uint32_t winsize;
  uint32_t offset;
  uint32_t editD;
} opt;
int main() {
  std::cout << "size\t" << sizeof(opt) << std::endl;
  std::cout << "align\t" << alignof(opt) << std::endl;
}