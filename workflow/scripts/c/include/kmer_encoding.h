#ifndef KMER_ENCODING_H_
#define KMER_ENCODING_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

uint64_t encode_kmer(const char* seq, int k);
char* decode_kmer(uint64_t val, int k);

#endif // KMER_ENCODING_H_
