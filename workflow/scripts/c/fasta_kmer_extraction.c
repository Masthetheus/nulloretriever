#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

#define MAX_SEQ 40000000
#define MAX_K 20
#define BUFFER_SIZE 65536

static unsigned char write_buffer[BUFFER_SIZE];
static size_t buffer_position = 0;

void flush_output_buffer(void){

    if (buffer_position > 0){
        fwrite(write_buffer, 1, buffer_position, stdout);
        buffer_position = 0;
    }
}
void print_packed_binary(uint64_t val, int bytes_per_half, int bytes_per_kmer, int l){
    if (buffer_position + bytes_per_kmer > BUFFER_SIZE){
        flush_output_buffer();
    }
    int bits_per_half = l * 2;
    
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
        flush_output_buffer();
    }

    for (int i = bytes_per_sequence - 1; i >= 0; i--) {
        unsigned char byte = (val >> (i * 8)) & 0xFF;
        write_buffer[buffer_position++] = byte;
    }
}
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

void print_binary_bytes(uint64_t val, int k) {
    // Output actual bytes, not ASCII characters
    for (int i = k * 2 - 1; i >= 0; i--) {
        unsigned char bit = (val >> i) & 1;
        fwrite(&bit, 1, 1, stdout);
    }
    fflush(stdout); // Force flush
}

void process_kmers(const char* seq, int seqlen, int k, char* seen, int bytes_per_sequence, int l) {

    for (int i = 0; i <= seqlen - k; i++) {
        uint64_t idx = encode_kmer(seq + i, k);
        if (idx == UINT64_MAX) {
            continue;
        }
        
        uint64_t byte = idx / 8, bit = idx % 8;
        if (!(seen[byte] & (1 << bit))) {
            print_packed_binary_test(idx, bytes_per_sequence);
            seen[byte] |= (1 << bit);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Uso: %s <arquivo_fasta> <k>\n", argv[0]);
        return 1;
    }
    
    fprintf(stderr, "DEBUG: Starting program with file=%s, k=%s\n", argv[1], argv[2]);
    
    FILE* f = fopen(argv[1], "r");
    if (!f) {
        perror("Erro ao abrir arquivo");
        return 1;
    }
    
    int k = atoi(argv[2]);
    fprintf(stderr, "DEBUG: k=%d\n", k);
    int l = k/2;
    uint64_t total = 1ULL << (2 * k);
    char* seen = calloc((total + 7) / 8, 1);
    char* seq = malloc(MAX_SEQ);
    int seqlen = 0;
    char line[1024];
    int line_count = 0;
    int bytes_per_sequence = (k * 2 + 7) / 8;

    while (fgets(line, sizeof(line), f)) {
        line_count++;
        if (line[0] == '>') {
            fprintf(stderr, "DEBUG: Found header at line %d: %.50s\n", line_count, line);
            if (seqlen > 0) {
                // Process forward strand
                process_kmers(seq, seqlen, k, seen, bytes_per_sequence, l);

                // Generate reverse complement
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
                process_kmers(revcomp_seq, seqlen, k, seen, bytes_per_sequence,l);
                free(revcomp_seq);

                seqlen = 0;
                memset(seen, 0, (total + 7) / 8);
            }
        } else {
            char* p = line;
            while (*p && *p != '\n' && *p != '\r') {
                seq[seqlen++] = *p++;
            }
        }
    }
    
    // Process last sequence
    if (seqlen > 0) {
        process_kmers(seq, seqlen, k, seen, bytes_per_sequence, l);
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
        process_kmers(revcomp_seq, seqlen, k, seen, bytes_per_sequence, l);
        free(revcomp_seq);
    }
    flush_output_buffer();
    fprintf(stderr, "DEBUG: Program finished\n");
    free(seen);
    free(seq);
    fclose(f);
    return 0;
}
