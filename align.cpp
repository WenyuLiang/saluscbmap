#include "align.h"
#include "edlib.h"
#include <cstdio>
#include <algorithm>
bool Align::cmp(const std::pair<uint32_t, int> &left,
                const std::pair<uint32_t, int> &right) {
  return left.second < right.second;
}

int Align::CandidateRef(const MappingMetadata &mapping_metadata,
                        SequenceBatch &ref, SequenceBatch &read,
                        uint32_t &i) {
  const int numRefHits = mapping_metadata.GetNumRefHits();
  bool isMapped = false;
  const int candidateLimit = 50;
  int candidateCount = 0;
  uint32_t flag = 0;
  
  std::vector<TripleType> repeat_hits;
  repeat_hits.reserve(candidateLimit); // random estimate
  //std::cout << "numRefHits: " << numRefHits << std::endl;
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
      result = edlibAlign(ref.GetSequenceAt(ref_id), ref.GetSequenceLengthAt(ref_id), read.GetSequenceAt(i),
                          read.GetSequenceLengthAt(i), edlibDefaultAlignConfig());

      if (result.status == EDLIB_STATUS_OK) {
        if (result.editDistance == 0) {
          valid_mapping_count_++;
          edlibFreeAlignResult(result);
          return 1;
        } else if (result.editDistance > 0 && result.editDistance < edit_distance_) { 
          // printf("edit_distance(%s, %s) = %d\n", ref.GetSequenceAt(ref_id),
          //        read.GetSequenceAt(i), result.editDistance);
                  
          if(flag & (1<<result.editDistance)) flag |= 1; // last bit means repeat
          flag |= 1<<result.editDistance;

          TripleType triple = {ref_id, i, result.editDistance};
          
          if(candidateCount >= candidateLimit) break;
          repeat_hits.emplace_back(triple);   
          candidateCount++;    
        }
      }
      edlibFreeAlignResult(result);
      topHits.pop();
    }
    
    if(repeat_hits.size() == 0){
      invalid_mapping_count_++;
      return 0;
    }

    // Sort the vector by the edit distance
    std::sort(repeat_hits.begin(), repeat_hits.end(), 
              [](const TripleType& a, const TripleType& b) {
                  return std::get<2>(a) < std::get<2>(b);
              });
    
    if(flag & 1){
        float x, y;
        float x0, y0;
        if (sscanf(ref.GetSequenceCommentAt(std::get<0>(repeat_hits[0])), "%f %f", &x0, &y0) != 2){
          std::cout << "sscanf failed" << std::endl;
        }
          for(auto it = repeat_hits.begin() + 1; it != repeat_hits.end(); ++it) {
              auto& triple = *it;
              if (sscanf(ref.GetSequenceCommentAt(std::get<0>(triple)), "%f %f", &x, &y) == 2){
                if(x0 - x < 2 && x - x0 < 2 && y0 - y < 2 && y - y0 < 2){
                  x0 = (x0 + x) / 2.0;
                  y0 = (y0 + y) / 2.0;
                  isMapped = true;
                  // replace ref in vicinity with the average of the two

                  // std::cout << "x0: " << x0 << " y0: " << y0  << "\t" << ref.GetSequenceCommentAt(std::get<0>(repeat_hits[0])) << "\tedit_distance(" << ref.GetSequenceAt(std::get<0>(repeat_hits[0])) << ", " << read.GetSequenceAt(std::get<1>(repeat_hits[0])) << ") = " << std::get<2>(repeat_hits[0]) << std::endl;
                  // std::cout << "x: " << x << " y: " << y  << "\t" << ref.GetSequenceCommentAt(std::get<0>(triple)) << "\tedit_distance(" << ref.GetSequenceAt(std::get<0>(triple)) << ", " << read.GetSequenceAt(std::get<1>(triple)) << ") = " << std::get<2>(triple) << "\n" <<std::endl;
                  }
                else {
                  //isMapped = false;
                  break;
                }
                  }else {
                    std::cout << "sscanf failed" << std::endl;
                  }
            }
          }else isMapped = true; // no repeat

    if (!isMapped) invalid_mapping_count_++;
    if (isMapped) valid_mapping_count_++;
    return 3;
  }
  return -1;
}