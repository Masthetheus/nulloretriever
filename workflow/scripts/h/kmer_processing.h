#ifndef KMER_PROCESSING_H_
#define KMER_PROCESSING_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>


void process_kmers(const char* seq, int seqlen, int k, char* seen, int bytes_per_sequence, int* tot);

#endif // KMER_PROCESSING_H_
