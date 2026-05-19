#ifndef KMERIO_H_
#define KMERIO_H_

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void write_decoded_kmer(FILE *f, uint64_t val, int k);
void flush_output_buffer(int *buff_tot);
void print_packed_binary(uint64_t val, int bytes_per_half, int bytes_per_kmer,
                         int half_k);
void print_packed_binary_test(uint64_t val, int bytes_per_sequence);
void print_binary_bytes(uint64_t val, int k);

#endif // KMERIO_H_
