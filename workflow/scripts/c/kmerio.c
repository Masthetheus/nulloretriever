#include "../h/kmerio.h"

#define BUFFER_SIZE 65536

static int buff_tot = 0;
static unsigned char write_buffer[BUFFER_SIZE];
static size_t buffer_position = 0;

void flush_output_buffer(int* buff_tot){
    if (buffer_position > 0){
        fwrite(write_buffer, 1, buffer_position, stdout);
        buffer_position = 0;
        ++ *buff_tot;
    }
}
void print_packed_binary(uint64_t val, int bytes_per_half, int bytes_per_kmer, int half_k){
    if (buffer_position + bytes_per_kmer > BUFFER_SIZE){
        flush_output_buffer(&buff_tot);
    }
    int bits_per_half = half_k * 2;

    for(int j = 1; j >= 0; j--) {
        uint64_t half = (val >> (j * bits_per_half)) & ((1ULL << bits_per_half) - 1);

        for (int i = bytes_per_half - 1; i >= 0; i--) {
            unsigned char byte = (half >> (i * 8)) & 0xFF;
            write_buffer[buffer_position++] = byte;
        }
    }
}
void print_packed_binary_test(uint64_t val, int bytes_per_sequence){
    if (buffer_position + bytes_per_sequence > BUFFER_SIZE){
        flush_output_buffer(&buff_tot);
    }

    for (int i = bytes_per_sequence - 1; i >= 0; i--) {
        unsigned char byte = (val >> (i * 8)) & 0xFF;
        write_buffer[buffer_position++] = byte;
    }
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
    char seq[k+1];
    for (int i = 0; i < k; i++) {
        int base_bits = (val >> (2 * (k - i - 1))) & 3;
        seq[i] = nt[base_bits];
    }
    seq[k] = '\0';
    fprintf(f, "%s\n", seq);
}
