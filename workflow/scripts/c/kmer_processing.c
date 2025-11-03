#include "../h/kmer_processing.h"
#include "../h/kmer_encoding.h"
#include "../h/kmerio.h"


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
        }
    }
}
