#include "../h/kmer_deduplication.h"

char* bitarray_deduplication(int k){
  uint64_t total = 1ULL << (2*k);
  char* seen = calloc((total+7) / 8, 1);
  return seen;
}

void deduplication_method(int k){
  if(k>=8 && k<=16)
    bitarray_deduplication(k);
  else if(k <= 20)
    chainhash_deduplication()
  else
    printf("Invalid k value!")
}
