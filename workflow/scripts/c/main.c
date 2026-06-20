#include "../h/kmer_extraction.h"
#include "../h/kmerio.h"
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ 8000000000
#define MAX_K 20

int main(int argc, char *argv[]) {
        if (argc != 3) {
                fprintf(stderr, "Usage: %s <fasta_file> <k>\n", argv[0]);
                return 1;
        }

        fprintf(stderr, "DEBUG: Starting program with file=%s, k=%s\n", argv[1],
                argv[2]);

        FILE *f = fopen(argv[1], "r");
        if (!f) {
                perror("Error opening file!");
                return 1;
        }

        int k = atoi(argv[2]);
        fprintf(stderr, "DEBUG: k=%d\n", k);
        uint64_t total = 1ULL << (2 * k);
        char *seen = calloc((total + 7) / 8, 1);
        char *seq = malloc(MAX_SEQ);
        int seqlen = 0;
        char line[1024];
        int line_count = 0;
        int bytes_per_sequence = (k * 2 + 7) / 8;
        int tot = 0;

        while (fgets(line, sizeof(line), f)) {
                line_count++;
                if (line[0] == '>') {
                        fprintf(stderr,
                                "DEBUG: Found header at line %d: %.50s\n",
                                line_count, line);
                        if (seqlen > 0) {
                                // Process forward strand
                                process_kmers(seq, seqlen, k, seen,
                                              bytes_per_sequence, &tot);
                                //process_rev(seq, seqlen, k, seen,
                                //              bytes_per_sequence, &tot);

                                // Generate and process reverse complement
                                //char *revcomp_seq = generate_revcomp_seq(seqlen, seq);
                                //process_kmers(revcomp_seq, seqlen, k, seen,
                                //              bytes_per_sequence, &tot);

                                //free(revcomp_seq);
                                seqlen = 0;
                        }
                } else {
                        char *p = line;
                        while (*p && *p != '\n' && *p != '\r') {
                                seq[seqlen++] = *p++;
                        }
                }
        }

        // Process last sequence
        if (seqlen > 0) {
                process_kmers(seq, seqlen, k, seen, bytes_per_sequence, &tot);
                //char *revcomp_seq = generate_revcomp_seq(seqlen, seq);
                //process_kmers(revcomp_seq, seqlen, k, seen,
                //              bytes_per_sequence, &tot);
                //free(revcomp_seq);
        }
        fprintf(stderr, "DEBUG: Program finished\n");
        free(seen);
        free(seq);
        fclose(f);
        fprintf(stderr, "DEBUG: Total of %d sequences inserted\n", tot);
        return 0;
}
