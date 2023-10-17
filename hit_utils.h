#pragma once

inline static uint32_t HitToSequenceIndex(uint64_t hit) { return (hit >> 32); }

inline static uint32_t HitToSequencePosition(uint64_t hit) {
  return (uint32_t)hit;
}

