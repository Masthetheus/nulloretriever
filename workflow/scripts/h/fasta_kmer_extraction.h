#ifndef FASTA_KMER_EXTRACTION_H_
#define FASTA_KMER_EXTRACTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

void flush_output_buffer(int* buff_tot);
void print_packed_binary(uint64_t val, int bytes_per_half, int bytes_per_kmer, int l);
void print_packed_binary_test(uint64_t val, int bytes_per_sequence);
uint64_t encode_kmer(const char* seq, int k);
void print_binary_bytes(uint64_t val, int k);
char* decode_kmer(uint64_t val, int k);
void write_decoded_kmer(FILE *f, uint64_t val, int k);
void process_kmers(const char* seq, int seqlen, int k, char* seen, int bytes_per_sequence, int* tot);
int main(int argc, char* argv[]);

#endif // FASTA_KMER_EXTRACTION_H_
