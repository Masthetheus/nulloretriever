#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>

#define MAX_SEQ 40000000
#define MAX_K 20

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

void print_packed_binary(uint64_t val, int k, int bytes_needed) {
    // Write in big-endian format
    for (int i = bytes_needed - 1; i >= 0; i--) {
        unsigned char byte = (val >> (i * 8)) & 0xFF;
        fwrite(&byte, 1, 1, stdout);
    }
}

void process_kmers(const char* seq, int seqlen, int k, char* seen) {
    int count = 0;
    int bytes_needed = (k * 2 + 7) / 8;
    int half_bytes_needed = bytes_needed/2;
    for (int i = 0; i <= seqlen - k; i++) {
        uint64_t idx = encode_kmer(seq + i, k);
        if (idx == UINT64_MAX) {
            continue;
        }
        
        uint64_t byte = idx / 8, bit = idx % 8;
        if (!(seen[byte] & (1 << bit))) {
            print_packed_binary(idx, k, bytes_needed);
            seen[byte] |= (1 << bit);
            count++;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Uso: %s <arquivo_fasta> <k> <genome_size>\n", argv[0]);
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
    uint64_t size = strtoull(argv[3], NULL, 10);
    uint64_t possible_kmers_org = (size - k) + 1;
    uint64_t possible_kmers_global = 1ULL << (2 * k);
    uint64_t total_bits;
    if (possible_kmers_global <= possible_kmers_org){
    total_bits = possible_kmers_global;
    } else{
    total_bits = possible_kmers_org;
    }
    char* seen = calloc((total_bits + 7) / 8, 1);
    char* seq = malloc(MAX_SEQ);
    int seqlen = 0;
    char line[1024];
    int line_count = 0;

    while (fgets(line, sizeof(line), f)) {
        line_count++;
        if (line[0] == '>') {
            fprintf(stderr, "DEBUG: Found header at line %d: %.50s\n", line_count, line);
            if (seqlen > 0) {
                // Process forward strand
                process_kmers(seq, seqlen, k, seen);

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
                process_kmers(revcomp_seq, seqlen, k, seen);
                free(revcomp_seq);

                seqlen = 0;
                memset(seen, 0, (total_bits + 7) / 8);
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
        process_kmers(seq, seqlen, k, seen);
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
        process_kmers(revcomp_seq, seqlen, k, seen);
        free(revcomp_seq);
    }
    
    fprintf(stderr, "DEBUG: Program finished\n");
    free(seen);
    free(seq);
    fclose(f);
    return 0;
}