#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAX_SEQ 40000000
#define MAX_K 20

// ...cabeçalho e defines...

// Função para codificar um k-mer linear em índice
uint64_t encode_kmer(const char* seq, int k) {
    uint64_t val = 0;
    for (int i = 0; i < k; i++) {
        val <<= 2;
        switch(seq[i]) {
            case 'A': val |= 0; break;
            case 'C': val |= 1; break;
            case 'G': val |= 2; break;
            case 'T': val |= 3; break;
            default: return UINT64_MAX;
        }
    }
    return val;
}
// Função para imprimir k-mer no formato numérico
void print_kmer_numeric(const char* seq, int k) {
    for (int j = 0; j < k; j++) {
        switch(seq[j]) {
            case 'A': putchar('0'); break;
            case 'C': putchar('1'); break;
            case 'G': putchar('2'); break;
            case 'T': putchar('3'); break;
            default: putchar('N'); // ou pule, se preferir
        }
        if (j < k - 1) putchar(',');
    }
}

// Função para processar todos os k-mers de uma sequência linear
void process_kmers(const char* seq, int seqlen, int k, char* seen) {
    for (int i = 0; i <= seqlen - k; i++) {
        uint64_t idx = encode_kmer(seq + i, k);
        if (idx == UINT64_MAX) continue;
        uint64_t byte = idx / 8, bit = idx % 8;
        if (!(seen[byte] & (1 << bit))) {
            print_kmer_numeric(seq + i, k); // <-- imprime no formato numérico
            putchar(' ');
            seen[byte] |= (1 << bit);
        }
    }
}

int main(int argc, char* argv[]) {
    // ...parsing de argumentos, abertura de arquivo, etc...
    if (argc != 3) {
        fprintf(stderr, "Uso: %s <arquivo_fasta> <k>\n", argv[0]);
        return 1;
    }
    FILE* f = fopen(argv[1], "r");
    if (!f) {
        perror("Erro ao abrir arquivo");
        return 1;
    }
    int k = atoi(argv[2]);
    uint64_t total = 1ULL << (2 * k);
    char* seen = calloc((total + 7) / 8, 1);
    char* seq = malloc(MAX_SEQ);
    int seqlen = 0;
    char line[1024];

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '>') {
            if (seqlen > 0) {
                // Processa k-mers da fita original
                process_kmers(seq, seqlen, k, seen);

                // Gera complementar reversa
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
                // Processa k-mers da complementar reversa
                process_kmers(revcomp_seq, seqlen, k, seen);
                free(revcomp_seq);

                seqlen = 0;
                memset(seen, 0, (total + 7) / 8);
            }
        } else {
            char* p = line;
            while (*p && *p != '\n' && *p != '\r') seq[seqlen++] = *p++;
        }
    }
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
    
    free(seen);
    free(seq);
    fclose(f);
    return 0;
}