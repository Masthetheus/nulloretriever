#ifndef FASTA_KMER_EXTRACTION_H_
#define FASTA_KMER_EXTRACTION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

void process_kmers(const char* seq, int seqlen, int k, char* seen, int bytes_per_sequence, int* tot);
int main(int argc, char* argv[]);

#endif // FASTA_KMER_EXTRACTION_H_
