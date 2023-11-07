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
#define OPTIONS "r:b:o:k:w:s:m:e:t:z:h"

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
  uint32_t batch_size;
} opt_t;

constexpr uint32_t kmer_size = 11;
constexpr uint32_t winsize = 5;
constexpr uint32_t offset = 3;
constexpr uint32_t margin = 5;
constexpr uint32_t edit_distance = 3;
constexpr uint32_t threads = 16;
constexpr uint32_t batch_size = 100000;

void usage(char *exec) {
  fprintf(stderr,
          "SYNOPSIS\n"
          " Mapped barcodes to spatial.fa and correct barcodes\n"
          "\n"
          "USAGE\n"
          " %s [-r spatial.fa] [-b barcode] [-o output] [-k kmer_size] [-w "
          "winsize] [-s offset] [-m margin] [-e edit_distance] [-t threads] "
          "[-z batch_size] [-h]\n"
          "\n"
          "OPTIONS\n"
          " -r input          Specify spatial.fa\n"
          " -b barcode        Specify barcode.fq\n"
          " -o output         Specify output file prefix or location\n"
          " -k kmer_size      Specify kmer size (default: 11)\n"
          " -w winsize        Specify window size (default: 5)\n"
          " -s offset         Specify offset value (default: 3)\n"
          " -m margin         Specify margin value (default: 5)\n"
          " -e edit_distance  Specify edit distance value (default: 3)\n"
          " -t threads        Specify number of threads (default: 16)\n"
          " -z batch_size     Specify batch size (default: 100000)\n"
          " -h                Display program help and usage\n",
          exec);
}

int main(int argc, char *argv[]) {
  opt_t opt = {"",      "",        "",     kmer_size,
               winsize, offset,    margin, edit_distance,
               threads, batch_size};
  int opt_;
  if (argc == 1) {
    usage(argv[0]);
    exit(0);
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
    case 'z':
      opt.batch_size = atoi(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(0);
    default:
      usage(argv[0]);
      exit(0);
    }
  }
  // --------------------------------------------construct index for
  // reference---------------------------------------//
  omp_set_num_threads(opt.threads);
  IndexParameters index_parameters = {opt.kmer,        opt.offset,
                                      opt.winsize + 1, opt.threads,
                                      opt.whitelist,   opt.batch_size};
  Index index(index_parameters);

  SequenceBatch ref_batch;
  ref_batch.InitializeLoading(index_parameters.reference_file_path);
  ref_batch.LoadAllRefSequences();

  index.Construct(ref_batch.GetNumSequences(), ref_batch);
  ref_batch.FinalizeLoading();
  //  --------------------------------------------construct index for
  //  reference---------------------------------------//
  std::string outbc = opt.output + "_bc.fq";
  std::string outwl = opt.output + "_wl.txt";

  std::ofstream bc(outbc);
  std::ofstream wl(outwl);
  std::ofstream spatial(opt.output + "_spatial.txt");
  std::ofstream stat(opt.output + "_stat.txt");
  // --------------------------------------------process
  // reads---------------------------------------//
  SequenceBatch read_batch(opt.batch_size);
  read_batch.InitializeLoading(opt.barcode);
  // read_batch.LoadAllReadSequences();

  MappingMetadata read_metadata(opt.margin);
  SeedGenerator seed_generator(
      opt.offset, opt.kmer, opt.winsize + 1); // for loop update seed_generator

  // Local variables for reduction
  uint32_t local_invalid_count = 0;
  uint32_t local_valid_count = 0;
  bool doneLoading = false;

  // // Get the stream buffer
  // std::streambuf* pBuffer = bc.rdbuf();

  //   // Define the new buffer size
  // const size_t newBufferSize = 1024 * 1024; // 1 MB, for example

  // // Allocate memory for the new buffer
  // char* buffer = new char[newBufferSize];

  // // Set the new buffer
  // pBuffer->pubsetbuf(buffer, newBufferSize);

  for (;;) {
    doneLoading = read_batch.LoadBatchReadSequences();
#pragma omp parallel for firstprivate(seed_generator, read_metadata)           \
    reduction(+ : local_invalid_count, local_valid_count)                      \
    shared(opt, ref_batch, read_batch)
    for (uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
      Align local_align(
          opt.editD); // Create a private instance of Align for each thread
      read_metadata.PrepareForMappingNextRead(20);
      seed_generator.GenerateSeeds(read_batch, i, read_metadata.seed_);
      int numCandidatePositions = index.GenerateCandidatePositions(
          read_metadata, ref_batch, read_batch, i);
      local_align.CandidateRef(read_metadata, ref_batch, read_batch, i);

      // Update local counters
      local_invalid_count += local_align.invalid_mapping_count_;
      local_valid_count += local_align.valid_mapping_count_;
    }
    for (uint32_t i = 0; i < read_batch.GetNumSequences(); ++i) {
      bc << ">" << read_batch.GetSequenceNameAt(i) << "\n"
         << read_batch.GetSequenceAt(i) << "\n+\n"
         << read_batch.GetSequenceQualAt(i) << "\n";
    }
    if (doneLoading) [[unlikely]]
      break;
  }
  read_batch.FinalizeLoading();
  // delete[] buffer;

  // local_invalid_count and local_valid_count are now the total counts
  stat << "map: " << local_valid_count << std::endl;
  stat << "unmap: " << local_invalid_count << std::endl;
  stat << "valid_bc_p: "
       << static_cast<float>(local_valid_count) /
              (local_valid_count + local_invalid_count)
       << std::endl;

  // std::cerr << "Writing output files..." << std::endl;

  for (uint32_t i = 0; i < ref_batch.GetNumSequences(); ++i) {
    if (ref_batch.ref_sequence_keep_[i]) {
      spatial << ref_batch.GetSequenceAt(i) << " "
              << ref_batch.GetSequenceCommentAt(i) << "\n";
      wl << ref_batch.GetSequenceAt(i) << "\n";
    }
  }
  bc.close();
  wl.close();
  spatial.close();
  stat.close();
  std::cerr << "Done Alignment!!!" << std::endl;
  std::cerr << "Releasing memory..." << std::endl;
  return 0;
}
