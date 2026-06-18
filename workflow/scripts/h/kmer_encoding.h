#ifndef KMER_ENCODING_H_
#define KMER_ENCODING_H_

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

uint64_t encode_kmer(const char *seq, int k);
uint64_t encode_rev(const char *seq, int k, int seqlen, int cont);
char *decode_kmer(uint64_t val, int k);

#endif // KMER_ENCODING_H_
