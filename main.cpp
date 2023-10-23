#include "seed.h"
#include "sequence_batch.h"
#include <iostream>
#include <fstream>
#include "index.h"
#include "mapping_metadata.h"
#include "align.h"
// void printBinary(uint64_t num) {
//     for (uint64_t i = 63; i != static_cast<uint64_t>(-1); --i) {
//         if (i == 32) {
//             std::cout << " ";
//         }
//         std::cout << ((num & (1ULL << i)) ? '1' : '0');
//     }
//     std::cout << std::endl;
// }
enum{
    kmer_size = 9,
    offset = 0,
};

int main(int argc, char * argv[]) {
    // --------------------------------------------construct index for reference---------------------------------------//
    IndexParameters index_parameters = {kmer_size, offset, 1, "wl.fa", "test.idx"};
    Index index(index_parameters);

    SequenceBatch ref_batch;
    ref_batch.InitializeLoading(index_parameters.reference_file_path);
    ref_batch.LoadAllSequences();

    index.Construct(ref_batch.GetNumSequences(), ref_batch);
    //index.CheckIndex(ref_batch.GetNumSequences(), ref_batch);
    ref_batch.FinalizeLoading();
    //exit(1);
     
    // --------------------------------------------construct index for reference---------------------------------------//
    
    // --------------------------------------------process reads---------------------------------------//
    SequenceBatch read_batch;
    read_batch.InitializeLoading("bctax.fq");
    read_batch.LoadAllSequences();
    MappingMetadata read_metadata;
    Align align;
    align.unmapped_reads_.reserve(read_batch.GetNumSequences());
    // int map = 0;
    for(uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
        read_metadata.PrepareForMappingNextRead(20);
        SeedGenerator seed_generator(offset, kmer_size, 30);
        seed_generator.GenerateSeeds(read_batch, i, read_metadata.seed_);
        int numCandidatePositions = index.GenerateCandidatePositions(read_metadata, ref_batch, read_batch, i);
        align.CandidateRef(read_metadata, ref_batch, read_batch, i);
    }
    read_batch.FinalizeLoading();
    std::cout << "map: " << align.valid_mapping_count_ << std::endl;
    std::cout << "unmap: " << align.invalid_mapping_count_ << std::endl;
    // write unmapped seq
    std::ofstream unmapped("unmapped.txt");
    for(const auto &read_name: align.unmapped_reads_) {
        unmapped << read_name << "\n";
    }
    unmapped.close();
    // --------------------------------------------process reads---------------------------------------//
    return 0;
}
