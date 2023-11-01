#include "align.h"
#include "edlib.h"
#include <algorithm>
#include <cstdio>
bool Align::cmp(const std::pair<uint32_t, int> &left,
                const std::pair<uint32_t, int> &right) {
  return left.second < right.second;
}

int Align::CandidateRef(const MappingMetadata &mapping_metadata,
                        SequenceBatch &ref, SequenceBatch &read, uint32_t &i) {
  const int numRefHits = mapping_metadata.GetNumRefHits();
  bool isMapped = false;
  const int candidateLimit = 50;
  int candidateCount = 0;
  uint32_t flag = 0;

  std::vector<TripleType> repeat_hits;
  repeat_hits.reserve(candidateLimit); // random estimate
  // std::cout << "numRefHits: " << numRefHits << std::endl;
  if (numRefHits == 0) {
    invalid_mapping_count_++;
    return 0;
  }
  EdlibAlignResult result;
  if (numRefHits > 0) {
    std::priority_queue<std::pair<int, uint32_t>,
                        std::vector<std::pair<int, uint32_t>>,
                        decltype(&Align::cmp)>
        topHits(Align::cmp);
    for (const auto &pair : mapping_metadata.reference_hit_count_) {
      topHits.push({pair.first, pair.second});
    }
    while (!topHits.empty()) {
      const auto &[ref_id, count_] = topHits.top();
      result =
          edlibAlign(ref.GetSequenceAt(ref_id), ref.GetSequenceLengthAt(ref_id),
                     read.GetSequenceAt(i), read.GetSequenceLengthAt(i),
                     edlibDefaultAlignConfig());

      if (result.status == EDLIB_STATUS_OK) {
        if (result.editDistance == 0) {
          valid_mapping_count_++;
          edlibFreeAlignResult(result);
#pragma omp critical
          {  
            if (!ref.ref_sequence_keep_modified_[ref_id]) {
              ref.ref_sequence_keep_[ref_id] = true;
              ref.ref_sequence_keep_modified_[ref_id] = true;
            }
          }
          return 1;
        } else if (result.editDistance > 0 &&
                   result.editDistance <= edit_distance_) {
          if (flag & (1 << result.editDistance))
            flag |= 1; // last bit means repeat
          flag |= 1 << result.editDistance;
          TripleType triple = {ref_id, i, result.editDistance};
          if (candidateCount >= candidateLimit)
            break;
          repeat_hits.emplace_back(triple);
          candidateCount++;
        }
      }
      edlibFreeAlignResult(result);
      topHits.pop();
    }

    if (repeat_hits.size() == 0) {
      invalid_mapping_count_++;
      return 0;
    }

    // Sort the vector by the edit distance
    std::sort(repeat_hits.begin(), repeat_hits.end(),
              [](const TripleType &a, const TripleType &b) {
                return std::get<2>(a) < std::get<2>(b);
              });

    if (flag & 1) {
      float x, y;
      float x0, y0;
      uint32_t ref_id_0 = std::get<0>(repeat_hits[0]); // first hit
      std::vector<uint32_t> ref_ids;
      ref_ids.reserve(repeat_hits.size());

      if (sscanf(ref.GetSequenceCommentAt(ref_id_0), "%f %f", &x0, &y0) != 2) {
        std::cout << "sscanf failed\t" << ref.GetSequenceCommentAt(ref_id_0)
                  << std::endl;
      }

      for (auto it = repeat_hits.begin() + 1; it != repeat_hits.end(); ++it) {
        auto &triple = *it;
        // std::cout << std::get<0>(triple) << "size: " <<
        // ref.GetSequenceCommentAt(std::get<0>(triple)) << std::endl;
        if (sscanf(ref.GetSequenceCommentAt(std::get<0>(triple)), "%f %f", &x,
                   &y) == 2) {
          // std::cout << "yes"<< std::endl;
          if (x0 - x < 2 && x - x0 < 2 && y0 - y < 2 && y - y0 < 2) {
            x0 = (x0 + x) / 2.0;
            y0 = (y0 + y) / 2.0;
            isMapped = true;
            ref_ids.emplace_back(std::get<0>(triple));
          } else {
            // isMapped = false;
            break;
          }
        } else {
          std::cout << "sscanf failed\t"
                    << ref.GetSequenceCommentAt(std::get<0>(repeat_hits[0]))
                    << std::endl;
        }
      }
      if (ref_ids.size() > 0) {
        ref_ids.emplace_back(ref_id_0);
        std::sort(ref_ids.begin(), ref_ids.end());
#pragma omp critical
        {        
          if (!ref.ref_sequence_keep_modified_[ref_ids[0]]){
            ref.ref_sequence_keep_[ref_ids[0]] = true;
            ref.ref_sequence_keep_modified_[ref_ids[0]] = true;
          }
            
          for (auto it = ref_ids.begin() + 1; it != ref_ids.end(); ++it) {
            ref.ref_sequence_keep_[*it] = false;
            if (!ref.ref_sequence_keep_modified_[*it])
              ref.ref_sequence_keep_modified_[*it] = true;
          }
        }
        read.ModifySequenceAt(i, ref.GetSequenceAt(ref_ids[0]),
                              ref.GetSequenceLengthAt(ref_ids[0]));
      }
    } else { // no repeat
      isMapped = true;
      uint32_t ref_id_0 = std::get<0>(repeat_hits[0]); // first hit
      read.ModifySequenceAt(i, ref.GetSequenceAt(ref_id_0),
                              ref.GetSequenceLengthAt(ref_id_0));
#pragma omp critical
      {        
        if (!ref.ref_sequence_keep_modified_[ref_id_0]){
          ref.ref_sequence_keep_[ref_id_0] = true;
          ref.ref_sequence_keep_modified_[ref_id_0] = true;
        }
      }
    }

    // if (isMapped) {
    //   for (auto it = repeat_hits.begin(); it != repeat_hits.end(); ++it) {
    //     auto &triple = *it;
    //     // print info
    //     std::cout << ref.GetSequenceCommentAt(std::get<0>(triple))
    //               << "\tedit_distance("
    //               << ref.GetSequenceAt(std::get<0>(triple)) << ", "
    //               << read.GetSequenceAt(std::get<1>(triple))
    //               << ") = " << std::get<2>(triple) << std::endl;
    //   }
    //   std::cout << "\n" << std::endl;
    // }

    if (!isMapped)
      invalid_mapping_count_++;
    if (isMapped)
      valid_mapping_count_++;
    return 3;
  }
  return -1;
}