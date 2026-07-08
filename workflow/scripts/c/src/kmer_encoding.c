#include "../include/kmer_encoding.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int8_t base_to_bits[256];

void init_base_table(void) {
    // Set all entries to -1 (invalid)
    memset(base_to_bits, -1, sizeof(base_to_bits));

    base_to_bits['A'] = base_to_bits['a'] = 0;
    base_to_bits['C'] = base_to_bits['c'] = 1;
    base_to_bits['T'] = base_to_bits['t'] = 2;
    base_to_bits['G'] = base_to_bits['g'] = 3;
}
uint64_t encode_kmer(const char *seq, int k) {
    init_base_table();
    uint64_t val = 0;
    for (int i = 0; i < k; i++) {
        int bits = base_to_bits[(unsigned char)seq[i]];
        if (bits < 0) {   // invalid character
            return UINT64_MAX;
        }
        val = (val << 2) | bits;
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
