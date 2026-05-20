#include "../h/kmerio.h"

#define BUFFER_SIZE 65536

void flush_output_buffer(unsigned char *write_buffer, size_t buffer_position) {
        //fprintf(stderr,"DEBUG: flush no buffer com tot = %ls.\n", buff_tot);
        fwrite(write_buffer, 1, buffer_position, stdout);
        buffer_position = 0;
        //++*buff_tot;
}

void print_binary_bytes(uint64_t val, int k) {
        // Output actual bytes, not ASCII characters
        for (int i = k * 2 - 1; i >= 0; i--) {
                unsigned char bit = (val >> i) & 1;
                fwrite(&bit, 1, 1, stdout);
        }
        fflush(stdout); // Force flush
}

void write_decoded_kmer(FILE *f, uint64_t val, int k) {
        static const char nt[4] = {'A', 'T', 'C', 'G'}; // match your encoding
        char seq[k + 1];
        for (int i = 0; i < k; i++) {
                int base_bits = (val >> (2 * (k - i - 1))) & 3;
                seq[i] = nt[base_bits];
        }
        seq[k] = '\0';
        fprintf(f, "%s\n", seq);
}
