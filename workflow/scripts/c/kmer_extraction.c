#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "../h/kmer_extraction.h"

#define MAX_SEQ 40000000

char* generate_revcomp_seq(int seqlen, char* seq){
  char* revcomp_seq = malloc(seqlen + 1);
  for (int i = 0; i < seqlen; i++){
    char b = seq[seqlen - 1 - i];
    switch(b){
      case 'A': revcomp_seq[i] = 'T'; break;
      case 'C': revcomp_seq[i] = 'G'; break;
      case 'G': revcomp_seq[i] = 'C'; break;
      case 'T': revcomp_seq[i] = 'A'; break;
      default: revcomp_seq[i] = 'N';
      }
    }
  revcomp_seq[seqlen] = '\0';
  return revcomp_seq;
  }

void process_chromossome(){}

char* read_fasta_body(char* line, char* seq, int seqlen){
  char* p = line;
  while (*p && *p != '\n' && *p != '\r') {
    seq[seqlen++] = *p++;
  return seq;
}

void extract_fasta_sequences(char* line, FILE* f, int* line_count){
  char* seq = malloc(MAX_SEQ);
  int seqlen = 0;
  while (fgets(line, sizeof(line), f)) {
    *line_count++;
    if (line[0] == '>') {
      fprintf(stderr, "DEBUG: Found header at line %d: %.50s\n", line_count, line);
      if (seqlen > 0) {
        // Process forward strand
        process_kmers(seq, seqlen, k, seen, bytes_per_sequence, &tot);
        // Generate reverse complement
        char* revcomp_seq = generate_revcomp_seq(seqlen, seq);
        process_kmers(revcomp_seq, seqlen, k, seen, bytes_per_sequence, &tot);
        free(revcomp_seq);
        seqlen = 0;
        memset(seen, 0, (total + 7) / 8);
      }
    } else {
      char* seq = read_fasta_body(line, seq, seqlen);
      }
    }
  }

  // Process last sequence
  if (seqlen > 0) {
    process_kmers(seq, seqlen, k, seen, bytes_per_sequence, &tot);
    char* revcomp_seq = malloc(seqlen + 1);
    for (int i = 0; i < seqlen; i++) {
      char b = seq[seqlen - 1 - i];
      switch(b) {
      case 'A': revcomp_seq[i] = 'T'; break;
      case 'C': revcomp_seq[i] = 'G'; break;
      case 'G': revcomp_seq[i] = 'C'; break;
      case 'T': revcomp_seq[i] = 'A'; break;
      default: revcomp_seq[i] = 'N';
      }
    }
    revcomp_seq[seqlen] = '\0';
    process_kmers(revcomp_seq, seqlen, k, seen, bytes_per_sequence, &tot);
    free(revcomp_seq);
  }
}
