#pragma once
#include "khash.h"
#include <cstdint>
#include <sys/time.h>
struct uint128_t {
  uint64_t first;
  uint64_t second;
};
KHASH_MAP_INIT_INT64(k128, uint128_t);
KHASH_MAP_INIT_INT64(k64_seq, uint64_t);
KHASH_SET_INIT_INT(k32_set);
KHASH_MAP_INIT_INT64(kmatrix, uint32_t);

inline static double GetRealTime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}

inline static uint64_t Hash64(uint64_t key, const uint64_t mask) {
  key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

static constexpr uint8_t char_to_uint8_table_[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static constexpr char uint8_to_char_table_[8] = {'A', 'C', 'G', 'T',
                                                 'N', 'N', 'N', 'N'};

inline static uint8_t CharToUint8(const char c) {
  return char_to_uint8_table_[(uint8_t)c];
}

inline static char Uint8ToChar(const uint8_t i) {
  return uint8_to_char_table_[i];
}