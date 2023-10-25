#include "align.h"
#include "index.h"
#include "mapping_metadata.h"
#include "seed.h"
#include "sequence_batch.h"
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <omp.h>
#define OPTIONS "r:b:o:k:w:s:m:e:t:h"

typedef struct {
  std::string whitelist;
  std::string barcode;
  std::string output;
  uint32_t kmer;
  uint32_t winsize;
  uint32_t offset;
  uint32_t margin;
  uint32_t editD;
  uint32_t threads;
} opt_t;

constexpr uint32_t kmer_size = 13;
constexpr uint32_t winsize = 0;
constexpr uint32_t offset = 3;
constexpr uint32_t margin = 5;
constexpr uint32_t edit_distance = 4;
constexpr uint32_t threads = 8;

void usage(char *exec) {
  fprintf(stderr,
          "SYNOPSIS\n"
          " Mapped barcodes to whitelist and correct barcodes\n"
          "\n"
          "USAGE\n"
          " %s [-r whitelist] [-b barcode] [-o output] [-k kmer_size] [-w "
          "winsize] [-s offset] [-m margin] [-e edit_distance] [-h]\n"
          "\n"
          "OPTIONS\n"
          " -r input          Specify whitelist.fa\n"
          " -b barcode        Specify barcode.fq\n"
          " -o output         Specify output file or location\n"
          " -k kmer_size      Specify kmer size (default: 13)\n"
          " -w winsize        Specify window size (default: 0)\n"
          " -s offset         Specify offset value (default: 3)\n"
          " -m margin         Specify margin value (default: 5)\n"
          " -e edit_distance  Specify edit distance value (default: 4)\n"
          " -t threads        Specify number of threads (default: 8)\n"
          " -h                Display program help and usage\n",
          exec);
}

int main(int argc, char *argv[]) {
  opt_t opt = {"", "", "", kmer_size, winsize, offset, margin, edit_distance, threads};
  int opt_;
  if (argc == 1) {
    usage(argv[0]);
    exit(1);
  }
  while ((opt_ = getopt(argc, argv, OPTIONS)) != -1) {
    switch (opt_) {
    case 'r':
      opt.whitelist = optarg;
      break;
    case 'b':
      opt.barcode = optarg;
      break;
    case 'o':
      opt.output = optarg;
      break;
    case 'k':
      opt.kmer = atoi(optarg);
      break;
    case 'w':
      opt.winsize = atoi(optarg);
      break;
    case 's':
      opt.offset = atoi(optarg);
      break;
    case 'm':
      opt.margin = atoi(optarg);
      break;
    case 'e':
      opt.editD = atoi(optarg);
      break;
    case 't':
      opt.threads = atoi(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
    default:
      usage(argv[0]);
      exit(1);
    }
  }
  // --------------------------------------------construct index for
  // reference---------------------------------------//
  IndexParameters index_parameters = {opt.kmer, opt.offset, opt.winsize, opt.whitelist,
                                      "test.idx"};
  Index index(index_parameters);

  SequenceBatch ref_batch;
  ref_batch.InitializeLoading(index_parameters.reference_file_path);
  ref_batch.LoadAllSequences();

  index.Construct(ref_batch.GetNumSequences(), ref_batch);
  // index.CheckIndex(ref_batch.GetNumSequences(), ref_batch);
  ref_batch.FinalizeLoading();
  // exit(1);

  // --------------------------------------------construct index for
  // reference---------------------------------------//

  // --------------------------------------------process
  // reads---------------------------------------//
  SequenceBatch read_batch;
  read_batch.InitializeLoading(opt.barcode);
  read_batch.LoadAllSequences();
  //MappingMetadata read_metadata(opt.margin);
  // MappingMetadata read_metadata;
  // Align align(opt.editD);
  // align.unmapped_reads_.reserve(read_batch.GetNumSequences());
  // SeedGenerator seed_generator(opt.offset, opt.kmer, 30); 
  // for (uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
  //   read_metadata.PrepareForMappingNextRead(20);
  //   seed_generator.GenerateSeeds(read_batch, i, read_metadata.seed_);
  //   int numCandidatePositions = index.GenerateCandidatePositions(
  //       read_metadata, ref_batch, read_batch, i);
  //   std::cout << "numCandidatePositions: " << numCandidatePositions
  //             << std::endl;
  //    //align.CandidateRef(read_metadata, ref_batch, read_batch, i); 
  // }

  // MappingMetadata read_metadata;
  //   Align align(4);
  //   align.unmapped_reads_.reserve(read_batch.GetNumSequences());
  //   // int map = 0;
  //   for(uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
  //       read_metadata.PrepareForMappingNextRead(20);
  //       SeedGenerator seed_generator(opt.offset, opt.kmer, 30);
  //       // std::cout << "kmer_size: " << kmer_size << std::endl;
  //       // std::cout << "offset: " << offset << std::endl;
        
  //       seed_generator.GenerateSeeds(read_batch, i, read_metadata.seed_);
  //       //std::cout << "seed num: "<<read_metadata.seed_.size()<< std::endl;
  //       int numCandidatePositions = index.GenerateCandidatePositions(read_metadata, ref_batch, read_batch, i);
  //       std::cout << "numCandidatePositions: " << numCandidatePositions << std::endl;
  //       align.CandidateRef(read_metadata, ref_batch, read_batch, i);
  //   }
     
     
omp_set_num_threads(opt.threads);
MappingMetadata read_metadata(opt.margin);
SeedGenerator seed_generator(opt.offset, opt.kmer, opt.winsize, 30);
Align align(opt.editD);
align.unmapped_reads_.reserve(read_batch.GetNumSequences());
#pragma omp parallel for firstprivate(seed_generator, read_metadata)           \
    shared(opt, ref_batch, read_batch)
  for (uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
    read_metadata.PrepareForMappingNextRead(20);
    seed_generator.GenerateSeeds(read_batch, i, read_metadata.seed_);
    int numCandidatePositions = index.GenerateCandidatePositions(
        read_metadata, ref_batch, read_batch, i);
    std::cout << "numCandidatePositions: " << numCandidatePositions
              << std::endl;
#pragma omp critical
    { 
      align.CandidateRef(read_metadata, ref_batch, read_batch, i); 
      }
  }
  read_batch.FinalizeLoading();
  std::cout << "map: " << align.valid_mapping_count_ << std::endl;
  std::cout << "unmap: " << align.invalid_mapping_count_ << std::endl;
  std::cout << "valid_bc_p: "
            << static_cast<float>(align.valid_mapping_count_) /
                   (align.valid_mapping_count_ + align.invalid_mapping_count_)
            << std::endl;
  // write unmapped seq
  std::ofstream unmapped("unmapped.txt");
  for (const auto &read_name : align.unmapped_reads_) {
    unmapped << read_name << "\n";
  }
  unmapped.close();
  // --------------------------------------------process
  // reads---------------------------------------//
  return 0;
}
