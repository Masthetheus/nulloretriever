#include "../include/kmer_extraction.h"
#include "../include/chained_bitarray.h"
#include "../include/kmer_encoding.h"
#include "../include/kmerio.h"
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ 40000000
#define BUFFER_SIZE 65536

static unsigned char write_buffer[BUFFER_SIZE];
static size_t  buffer_position = 0;
static int buffers_sent=0;

char *generate_revcomp_seq(int seqlen, char *seq) {
        char *revcomp_seq = malloc(seqlen + 1);
        for (int i = 0; i < seqlen; i++) {
                char b = seq[seqlen - 1 - i];
                switch (b) {
                case 'A':
                        revcomp_seq[i] = 'T';
                        break;
                case 'C':
                        revcomp_seq[i] = 'G';
                        break;
                case 'G':
                        revcomp_seq[i] = 'C';
                        break;
                case 'T':
                        revcomp_seq[i] = 'A';
                        break;
                default:
                        revcomp_seq[i] = 'N';
                }
        }
        revcomp_seq[seqlen] = '\0';
        return revcomp_seq;
}

void process_kmers(const char *seq, int seqlen, int k, char *seen,
                   int bytes_per_sequence, int *tot) {


        for (int i = 0; i <= seqlen - k; i++) {
                uint64_t idx = encode_kmer(seq + i, k);
                uint64_t idx_comp = encode_rev(idx, k);
                if (idx == UINT64_MAX) {
                        continue;
                }

                uint64_t byte = idx / 8, bit = idx % 8;
                if (!(seen[byte] & (1 << bit))) {
                        if (buffer_position + bytes_per_sequence > BUFFER_SIZE) {
                                flush_output_buffer(write_buffer,&buffer_position);
                                buffer_position = 0;
                                buffers_sent++;
                        }
                        for (int i = bytes_per_sequence - 1; i >= 0; i--) {
                                unsigned char byte = (idx >> (i * 8)) & 0xFF;
                                write_buffer[buffer_position++] = byte;
                        }

                        seen[byte] |= (1 << bit);
                        *tot += 1;
                        int total = *tot;
                        fprintf(stderr, "%d", total);
                }
                byte = idx_comp / 8;
                bit = idx_comp % 8;
                if (!(seen[byte] & (1 << bit))) {
                        if (buffer_position + bytes_per_sequence > BUFFER_SIZE) {
                                flush_output_buffer(write_buffer,&buffer_position);
                                buffer_position = 0;
                                buffers_sent++;
                        }
                        for (int i = bytes_per_sequence - 1; i >= 0; i--) {
                                unsigned char byte = (idx_comp >> (i * 8)) & 0xFF;
                                write_buffer[buffer_position++] = byte;
                        }

                        seen[byte] |= (1 << bit);
                        *tot += 1;
                        int total = *tot;
                        fprintf(stderr, "%d", total);
                }
        }
        if (buffer_position > 0){
                flush_output_buffer(write_buffer, &buffer_position);
                buffer_position = 0;
                buffers_sent++;
        }
        fprintf(stderr, "DEBUG: até o momento foram enviados %d buffers.\n\n", buffers_sent);

}

//void process_rev(const char *seq, int seqlen, int k, char *seen,
//                   int bytes_per_sequence, int *tot) {
//
//        int cont = 0;
//        for (int i = 0; i <= seqlen - k; i++) {
//                uint64_t idx = encode_rev(seq, k, seqlen, cont++);
//                if (idx == UINT64_MAX) {
//                        continue;
//                }
//
//                uint64_t byte = idx / 8, bit = idx % 8;
//                if (!(seen[byte] & (1 << bit))) {
//                        if (buffer_position + bytes_per_sequence > BUFFER_SIZE) {
//                                flush_output_buffer(write_buffer,&buffer_position);
//                                buffer_position = 0;
//                                buffers_sent++;
//                        }
//                        for (int i = bytes_per_sequence - 1; i >= 0; i--) {
//                                unsigned char byte = (idx >> (i * 8)) & 0xFF;
//                                write_buffer[buffer_position++] = byte;
//                        }
//
//                        seen[byte] |= (1 << bit);
//                        *tot += 1;
//                        int total = *tot;
//                        fprintf(stderr, "%d", total);
//                }
//        }
//        if (buffer_position > 0){
//                flush_output_buffer(write_buffer, &buffer_position);
//                buffer_position = 0;
//                buffers_sent++;
//        }
//        fprintf(stderr, "DEBUG: até o momento foram enviados %d buffers.\n\n", buffers_sent);
//
//}
