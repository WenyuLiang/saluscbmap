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
constexpr uint32_t edit_distance = 3;
constexpr uint32_t threads = 16;

void usage(char *exec) {
  fprintf(stderr,
          "SYNOPSIS\n"
          " Mapped barcodes to spatial.fa and correct barcodes\n"
          "\n"
          "USAGE\n"
          " %s [-r spatial.fa] [-b barcode] [-o output] [-k kmer_size] [-w "
          "winsize] [-s offset] [-m margin] [-e edit_distance] [-h]\n"
          "\n"
          "OPTIONS\n"
          " -r input          Specify spatial.fa\n"
          " -b barcode        Specify barcode.fq\n"
          " -o output         Specify output file or location\n"
          " -k kmer_size      Specify kmer size (default: 13)\n"
          " -w winsize        Specify window size (default: 0)\n"
          " -s offset         Specify offset value (default: 3)\n"
          " -m margin         Specify margin value (default: 5)\n"
          " -e edit_distance  Specify edit distance value (default: 3)\n"
          " -t threads        Specify number of threads (default: 16)\n"
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
  omp_set_num_threads(opt.threads);
  IndexParameters index_parameters = {opt.kmer, opt.offset, opt.winsize + 1, opt.threads, opt.whitelist};
  Index index(index_parameters);

  SequenceBatch ref_batch;
  ref_batch.InitializeLoading(index_parameters.reference_file_path);
  ref_batch.LoadAllSequences();
  
  index.Construct(ref_batch.GetNumSequences(), ref_batch);
  ref_batch.FinalizeLoading();
  
  // exit(1);

  // --------------------------------------------construct index for
  // reference---------------------------------------//

  // --------------------------------------------process
  // reads---------------------------------------//
  SequenceBatch read_batch;
  read_batch.InitializeLoading(opt.barcode);
  read_batch.LoadAllSequences();

MappingMetadata read_metadata(opt.margin);
SeedGenerator seed_generator(opt.offset, opt.kmer, opt.winsize + 1, 30); // for loop update seed_generator

// Local variables for reduction
uint32_t local_invalid_count = 0;
uint32_t local_valid_count = 0;

#pragma omp parallel for firstprivate(seed_generator, read_metadata) \
    reduction(+:local_invalid_count, local_valid_count) \
    shared(opt, ref_batch, read_batch)
for (uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
    Align local_align(opt.editD);  // Create a private instance of Align for each thread
    read_metadata.PrepareForMappingNextRead(20);
    seed_generator.GenerateSeeds(read_batch, i, read_metadata.seed_);
    int numCandidatePositions = index.GenerateCandidatePositions(
        read_metadata, ref_batch, read_batch, i);
    local_align.CandidateRef(read_metadata, ref_batch, read_batch, i);
    
    // Update local counters
    local_invalid_count += local_align.invalid_mapping_count_;
    local_valid_count += local_align.valid_mapping_count_;
}

//Update the align object with the aggregated results
// align.invalid_mapping_count_ += local_invalid_count;
// align.valid_mapping_count_ += local_valid_count;

read_batch.FinalizeLoading();
// local_invalid_count and local_valid_count are now the total counts
std::cout << "map: " << local_valid_count << std::endl;
std::cout << "unmap: " << local_invalid_count << std::endl;
std::cout << "valid_bc_p: "
          << static_cast<float>(local_valid_count) /
                 (local_valid_count + local_invalid_count)
          << std::endl;

  // write unmapped seq
  // std::ofstream unmapped("unmapped.txt");
  // for (const auto &read_name : align.unmapped_reads_) {
  //   unmapped << read_name << "\n";
  // }
  // unmapped.close();
  // --------------------------------------------process
  // reads---------------------------------------//
  return 0;
}
