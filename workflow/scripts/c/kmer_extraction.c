#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "../h/kmer_extraction.h"
#include "../h/kmer_encoding.h"
#include "../h/kmerio.h"

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

void process_kmers(const char* seq, int seqlen, int k, char* seen, int bytes_per_sequence, int* tot) {
    for (int i = 0; i <= seqlen - k; i++) {
        uint64_t idx = encode_kmer(seq + i, k);
        if (idx == UINT64_MAX) {
            continue;
        }

        uint64_t byte = idx / 8, bit = idx % 8;
        if (!(seen[byte] & (1 << bit))) {
            print_packed_binary_test(idx, bytes_per_sequence);
            seen[byte] |= (1 << bit);
            *tot += 1;
            int total = *tot;
            printf("%d", total);
        }
    }
}
char* read_fasta_body(char* line, char* seq, int* seqlen) {
    char* p = line;
    int current_len = *seqlen;

    while (*p && *p != '\n' && *p != '\r') {
        seq[current_len] = *p;
        current_len++;
        p++;
    }
    *seqlen = current_len;  // Update the length
    return seq;
}
