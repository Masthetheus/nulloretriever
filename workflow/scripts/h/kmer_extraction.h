#ifndef KMER_EXTRACTION_H_
#define KMER_EXTRACTION_H_

#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *generate_revcomp_seq(int seqlen, char *seq);
void process_kmers(const char *seq, int seqlen, int k, char *seen,
                   int bytes_per_sequence, int *tot);
//void process_rev(const char *seq, int seqlen, int k, char *seen,
//                   int bytes_per_sequence, int *tot);
char *read_fasta_body(char *line, char *seq, int *seqlen);

#endif // KMER_EXTRACTION_H_
