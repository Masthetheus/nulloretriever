#include "../h/kmer_deduplication.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

char* bitarray_deduplication(int k){
  uint64_t total = 1ULL << (2*k);
  char* seen = calloc((total+7) / 8, 1);
  return seen;
}

char* deduplication_method(int k){
  if(k>=8 && k<=16){
    return bitarray_deduplication(k);
  }
  else if(k <= 20)
    chainhash_deduplication()
  else
    printf("Invalid k value!")
}
