#include "edlib.h"
#include "align.h"
bool Align::cmp(const std::pair<uint32_t, int>& left, const std::pair<uint32_t, int>& right) {
    return left.second < right.second; 
}

// struct myComp {
//     constexpr bool operator()(
//         std::pair<uint32_t, int> const& a,
//         std::pair<uint32_t, int> const& b)
//         const noexcept
//     {
//         return a.second < b.second;
//     }
// };

int Align::CandidateRef(const MappingMetadata &mapping_metadata, const SequenceBatch &ref, SequenceBatch &read, uint32_t &i){
    int numRefHits = mapping_metadata.GetNumRefHits();
    bool isMapped = false;
    std::cout << "numRefHits: " << numRefHits << std::endl;
    if(numRefHits == 0) {
        unmapped_reads_.emplace_back(read.GetSequenceAt(i));
        invalid_mapping_count_++;
        return 0;
    }
    std::vector<EdlibAlignResult> alignments;
    alignments.reserve(numRefHits);
    EdlibAlignResult result;
    if(numRefHits > 0 && numRefHits <= 5) {
        for(const auto &pair: mapping_metadata.reference_hit_count_){
            //std::cout << "reference_id: " << pair.first  << " reference_seq: " << ref.GetSequenceAt(pair.first)<< " count: " << pair.second << "\n" << std::endl;
            result = edlibAlign(ref.GetSequenceAt(pair.first), 30, read.GetSequenceAt(i), 30, edlibDefaultAlignConfig());
            if (result.status == EDLIB_STATUS_OK) {
                printf("edit_distance(%s, %s) = %d\n", ref.GetSequenceAt(pair.first), read.GetSequenceAt(i), result.editDistance);
                if(result.editDistance <= 3) {         
                    isMapped = true;
                }
            }
            edlibFreeAlignResult(result);
        }
        if(!isMapped) {
            unmapped_reads_.emplace_back(read.GetSequenceAt(i));
            invalid_mapping_count_++;
        }
        if(isMapped) valid_mapping_count_++;
        return numRefHits;
    } else if(numRefHits > 5){
        int originalSize = mapping_metadata.reference_hit_count_.size();
        std::priority_queue<std::pair<int, uint32_t>, std::vector<std::pair<int, uint32_t>>, decltype(&Align::cmp)> topHits(Align::cmp);
        //std::priority_queue<std::pair<int, uint32_t>, std::vector<std::pair<int, uint32_t>>, myComp> topHits;
        for(const auto &pair: mapping_metadata.reference_hit_count_){
            //std::cout << "first: " << pair.first << " second: " << pair.second << std::endl;
            if(pair.second >5) topHits.push({pair.first, pair.second});            
            // if (topHits.size() > 15) {
            //     topHits.pop();
            // }
        } 
        //uint32_t ref_count = 0;
        while (!topHits.empty()) {
            const auto& [ref_id, count_] = topHits.top();
            //std::cout << "ref_id: " << ref_id << " count: " << count_ << std::endl;
            result = edlibAlign(ref.GetSequenceAt(ref_id), 30, read.GetSequenceAt(i), 30, edlibDefaultAlignConfig());
            if (result.status == EDLIB_STATUS_OK) {
                //std::cout<<"ref id: "<<ref_id<<"read id: "<<i<<std::endl;
                printf("edit_distance(%s, %s) = %d\n", ref.GetSequenceAt(ref_id), read.GetSequenceAt(i), result.editDistance);
                if(result.editDistance <= 3) {
                    isMapped = true;
                }
            }
            edlibFreeAlignResult(result);
            topHits.pop();
            // ref_count++;
            // if(ref_count == 8) break;
        }
        if(!isMapped) {
            unmapped_reads_.emplace_back(read.GetSequenceAt(i));
            invalid_mapping_count_++;
        }
        if(isMapped) valid_mapping_count_++;
        return 3;
    }
    return -1;
}