#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define MAX_SEQ 36
#define EMPTY_SLOT UINT64_MAX
#define INITIAL_CAPACITY 100000

static int8_t base_to_bits[256];

void init_base_table(void) {
    // Set all entries to -1 (invalid)
    memset(base_to_bits, -1, sizeof(base_to_bits));

    base_to_bits['A'] = base_to_bits['a'] = 0;
    base_to_bits['C'] = base_to_bits['c'] = 1;
    base_to_bits['T'] = base_to_bits['t'] = 2;
    base_to_bits['G'] = base_to_bits['g'] = 3;
}

static const char BITS_TO_BASE[4] = {'A', 'C', 'T', 'G'};

uint64_t encode_kmer(const char *seq, int k) {
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

void decode_kmer(uint64_t idx, int k, char *seq) {
    for (int i = 0; i < k; i++) {
        uint64_t base = (idx >> ((k-i-1)*2)) & 3;
        seq[i] = BITS_TO_BASE[base];
    }
    seq[k] = '\0';
}

static void generate(int k, int d, int start, int *current, int depth, int *data, int *idx){
        if (depth == d){
                for (int i = 0; i < d; i++){
                        data[(*idx) * d + i] = current[i];
                }
                (*idx)++;
                return;
        }
        for (int i = start; i <= k - (d - depth); i++){
                current[depth] = i;
                generate(k, d, i + 1, current, depth + 1, data, idx);
        }
}

void precompute_combinations(int k, int d, int **combos, int *num_combos){
        long long num = 1;
        for (int i = 1; i <= d; i++){
                num = num * (k - d + i) / i;
        }
        *num_combos = (int)num;

        int *data = malloc((size_t)(*num_combos) * d * sizeof(int));
        if (!data){
                *num_combos = 0;
                *combos = NULL;
                return;
        }

        int *current = malloc((size_t)d * sizeof(int));
        if(!current){
                free(data);
                *num_combos = 0;
                *combos = NULL;
                return;
        }

        int idx = 0;
        generate(k, d, 0, current, 0, data, &idx);

        free(current);
        *combos = data;
}

bool seen_countains(uint64_t idx, char *seen){
        size_t byte = idx/8;
        int bit = idx%8;
        if ((seen[byte]&(1 << (bit))) != 0){
                return true;
        }else{
                return false;
        }
}

bool all_neighbors_exist(uint64_t kmer, int d, const int *combos, int num_combos, char *seen){

        for (int c = 0; c < num_combos; c++){
                int posicoes[d];
                for (int p = 0; p < d; p++){
                        posicoes[p] = combos[c * d + p];
                }

                int max_pat = 1 << (2*d);
                for (int pat = 0; pat < max_pat; pat++){
                        uint64_t neighbor = kmer;
                        bool valid = true;

                        for (int p = 0; p < d; p++){
                                int base = (pat >> (2 * p)) & 3;
                                int pos = posicoes[p];
                                int original = (int)((kmer >> (2 * pos)) & 3);

                                if (base == original){
                                        valid = false;
                                        break;
                                }

                                neighbor &= ~(3ULL << (2 * pos));
                                neighbor |= ((uint64_t)base << (2 * pos));
                        }

                        if (!valid) continue;

                        if (!seen_countains(neighbor, seen)){
                                return false;
                        }

                }
        }
        return true;
}

int main(int argc, char *argv[]){
        init_base_table();

        if (argc != 3){
                fprintf(stderr, "Usage: %s <file> <order>\n", argv[0]);
                return 1;
        }

        FILE *f = fopen(argv[1], "rb");
        if (!f){
                perror("Error opening file!");
                return 1;
        }
        int order = atoi(argv[2]);

        uint64_t *nullomers = NULL;
        size_t count = 0;
        size_t capacity = 0;

        fseek(f, 6, SEEK_SET);

        uint16_t parameters[4];
        fread(parameters,2,2, f);

        int k = parameters[0];
        int half_k = parameters[1];

        uint8_t codes[2];
        fread(codes, 1, 2, f);
        unsigned char byte_size = codes[0];
        unsigned char counter_size = codes[1];

        uint64_t total = 1ULL << (2*k);
        char *seen = calloc((total+7) / 8, 1);

        uint64_t v1 = 0;
        while (1) {
                v1 = 0;
                if (fread(&v1, byte_size, 1, f) != 1) break;

                uint64_t counter = 0;
                fread(&counter, counter_size, 1, f);
                uint64_t null_count = counter;

                uint8_t *null_collection = malloc(null_count * byte_size);
                fread(null_collection, byte_size, null_count, f);

                for (size_t i = 0; i < null_count; i++) {
                        uint64_t v2 = 0;
                        memcpy(&v2, null_collection + i * byte_size, byte_size);
                        uint64_t idx = (v1 << ((uint64_t) half_k * 2)) | v2;

                        size_t byte = idx/8;
                        int bit = idx%8;
                        seen[byte] |= (1 << bit);

                        if (count == capacity){
                                size_t new_capacity = (capacity==0) ? INITIAL_CAPACITY : capacity * 2;
                                uint64_t *tmp = realloc(nullomers, new_capacity * sizeof(uint64_t));
                                if (!tmp) {
                                        perror("Realloc failed");
                                        free(nullomers);
                                        fclose(f);
                                        return 1;
                                }
                                nullomers = tmp;
                                capacity = new_capacity;
                        }
                        nullomers[count++] = idx;
                }
                free(null_collection);
        }

        printf("%zu nullomers allocated initially\n", count);

        uint64_t *candidates = malloc(count * sizeof(uint64_t));
        memcpy(candidates, nullomers, count * sizeof(uint64_t));
        size_t cand_count = count;

        free(nullomers);
        nullomers = NULL;

        for (int i = 1; i <= order; i++){
                int *combos_d = NULL;
                int num_combos_d = 0;
                precompute_combinations(k, i, &combos_d, &num_combos_d);

                uint64_t *passed = NULL;
                size_t passed_count = 0;
                size_t passed_cap = 0;

                for (size_t j = 0; j < cand_count; j++){
                        if (all_neighbors_exist(candidates[j], i, combos_d, num_combos_d, seen)){
                                if (passed_count == passed_cap){
                                        passed_cap = (passed_cap ==0) ? 256 : passed_cap * 2;
                                        uint64_t *tmp = realloc(passed, passed_cap * sizeof(uint64_t));
                                        if (!tmp) {exit(1);}
                                        passed = tmp;
                                }
                                passed[passed_count++] = candidates[j];
                        }
                }
                
                free(combos_d);
                combos_d = NULL;
                free(candidates);

                if (passed_count == 0){
                        printf("No higher order null\n");
                        candidates = NULL;
                        break;
                }

                candidates = passed;
                cand_count = passed_count;
                printf("%zu passed for order %d\n", cand_count, i);
                

                FILE *out;
                char buffer[30];

                snprintf(buffer, sizeof(buffer), "organism_order_%d.txt", i);

                out = fopen(buffer,"w");
                char *seq = malloc((k+1) * sizeof(char));
                for(size_t m = 0; m < cand_count; m++){
                        decode_kmer(candidates[m], k, seq);
                        fprintf(out, "%s\n", seq);
                }
                free(seq);
                fclose(out);
        }

        fclose(f);
        
        printf("Loaded %zu nullomers.\n", count);

        free(candidates);
        free(seen);
        return 0;
        
}
