#include "../include/kmer_encoding.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

uint64_t encode_kmer(const char* seq, int k) {
    uint64_t val = 0;
    for (int i = 0; i < k; i++) {
        val <<= 2;
        switch(seq[i]) {
            case 'A': val |= 0; break;
            case 'C': val |= 2; break;
            case 'G': val |= 3; break;
            case 'T': val |= 1; break;
            default:
                return UINT64_MAX;
        }
    }
    return val;
}

char* decode_kmer(uint64_t val, int k) {
    static const char nt[4] = {'A', 'T', 'C', 'G'}; // Adjust to your encoding
    char* seq = malloc(k + 1);
    if (!seq) return NULL;
    for (int i = 0; i < k; i++) {
        int base_bits = (val >> (2 * (k - i - 1))) & 3;
        seq[i] = nt[base_bits];
    }
    seq[k] = '\0';
    return seq;
}
