#include "CommandLines.h"
#include "ketopt.h"
static ko_longopt_t long_options[] = {
    {"kmer", ko_required_argument, 300},
    {"winsize", ko_required_argument, 301},
    {"offset", ko_required_argument, 302},
    {"editD", ko_required_argument, 303}} opt_t;
