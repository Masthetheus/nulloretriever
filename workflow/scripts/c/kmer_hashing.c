#include "../h/kmer_hashing.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Node {
        char *key;
        char *value;
        struct Node *next;
} Node;

typedef struct HashItem {
        uint64_t *key;
        uint64_t value;
} HashItem;

typedef struct HashTable {
        HashItem **items;
        int size;
        int count;
} HashTable;

static uint64_t hash_funct(const uint64_t *idx){uint64_t hash_value = idx >> 2}
