#include "../h/kmer_encoding.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

uint64_t encode_kmer(const char *seq, int k) {
        uint64_t val = 0;
        for (int i = 0; i < k; i++) {
                val <<= 2;
                val |= (seq[i]>>1) & 3; 
                }
        return val;
}

uint64_t encode_rev(uint64_t val, int k) {
        uint64_t revval = 0;
        for (int i = 0; i < k; i++) {
                uint64_t base = val & 3;
                uint64_t basecomp = base ^ 2;
                revval = (revval << 2) | basecomp; 
                val >>= 2;
                }
        return revval;
}


char *decode_kmer(uint64_t val, int k) {
        static const char nt[4] = {'A', 'C', 'T',
                                   'G'}; // Adjust to your encoding
        char *seq = malloc(k + 1);
        if (!seq)
                return NULL;
        for (int i = 0; i < k; i++) {
                int base_bits = (val >> (2 * (k - i - 1))) & 3;
                seq[i] = nt[base_bits];
        }
        seq[k] = '\0';
        return seq;
}
