#include "edlib.h"
#include "align.h"
bool Align::cmp(const std::pair<uint32_t, int>& left, const std::pair<uint32_t, int>& right) {
    return left.second < right.second; 
}

int Align::CandidateRef(MappingMetadata &mapping_metadata, SequenceBatch &ref, SequenceBatch &read, uint32_t &i){
    int numRefHits = mapping_metadata.GetNumRefHits();
    bool isMapped = false;
    std::cout << "numRefHits: " << numRefHits << std::endl;
    if(numRefHits == 0) {
        unmapped_reads_.emplace_back(read.GetSequenceAt(i));
        return 0;
    }
    std::vector<EdlibAlignResult> alignments;
    alignments.reserve(numRefHits);
    EdlibAlignResult result;
    // if(numRefHits > 0 && numRefHits <=3) {
    //     for(const auto &pair: mapping_metadata.reference_hit_count_){
    //         //std::cout << "reference_id: " << pair.first  << " reference_seq: " << ref.GetSequenceAt(pair.first)<< " count: " << pair.second << "\n" << std::endl;
    //         result = edlibAlign(ref.GetSequenceAt(pair.first), 30, read.GetSequenceAt(i), 30, edlibDefaultAlignConfig());
    //         if (result.status == EDLIB_STATUS_OK) {
    //             printf("edit_distance(%s, %s) = %d\n", ref.GetSequenceAt(pair.first), read.GetSequenceAt(i), result.editDistance);
    //             if(result.editDistance <= 3) {         
    //                 isMapped = true;
    //             }
    //         }
    //         edlibFreeAlignResult(result);
    //     }
    //     if(!isMapped) unmapped_reads_.emplace_back(read.GetSequenceAt(i));
    //     if(isMapped) valid_mapping_count_++;
    //     return numRefHits;
    // }
    // // Only keeps 3 hits with largest counts
    // if(numRefHits > 3){
    //     std::priority_queue<std::pair<int, uint32_t>, std::vector<std::pair<int, uint32_t>>, decltype(&Align::cmp)> topHits(Align::cmp);
    //     for(const auto &pair: mapping_metadata.reference_hit_count_){
    //         topHits.push({pair.first, pair.second});            
    //         if (topHits.size() > 3) {
    //             topHits.pop();
    //         }
    //     }
    //     while (!topHits.empty()) {
    //         const auto& [ref_id, count_] = topHits.top();
    //         topHits.pop();
    //         result = edlibAlign(ref.GetSequenceAt(ref_id), 30, read.GetSequenceAt(i), 30, edlibDefaultAlignConfig());
    //         if (result.status == EDLIB_STATUS_OK) {
    //             printf("edit_distance(%s, %s) = %d\n", ref.GetSequenceAt(ref_id), read.GetSequenceAt(i), result.editDistance);
    //             if(result.editDistance <= 3) {
    //                 isMapped = true;
    //             }
    //         }
    //         edlibFreeAlignResult(result);
    //     }
    //     if(!isMapped) unmapped_reads_.emplace_back(read.GetSequenceAt(i));
    //     if(isMapped) valid_mapping_count_++;
    //     return 3;
    // }
    // return -1;
    if(numRefHits > 0) {
        for(const auto &pair: mapping_metadata.reference_hit_count_){
            //std::cout << "reference_id: " << pair.first  << " reference_seq: " << ref.GetSequenceAt(pair.first)<< " count: " << pair.second << "\n" << std::endl;
            result = edlibAlign(ref.GetSequenceAt(pair.first), 30, read.GetSequenceAt(i), 30, edlibDefaultAlignConfig());
            if (result.status == EDLIB_STATUS_OK) {
                printf("edit_distance(%s, %s) = %d\n", ref.GetSequenceAt(pair.first), read.GetSequenceAt(i), result.editDistance);
                if(result.editDistance <= 1) {         
                    isMapped = true;
                    break;
                }
                if(result.editDistance <= 2) {         
                    isMapped = true;
                }
            }
            edlibFreeAlignResult(result);
        }
        if(!isMapped) unmapped_reads_.emplace_back(read.GetSequenceAt(i));
        if(isMapped) valid_mapping_count_++;
        return numRefHits;
    }
}