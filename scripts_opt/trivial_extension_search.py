"""Script for trivial extension search."""

import bitarray

from nulloretriever.analysis.processing import mount_trie_from_bitfile

smaller_file = "novoresult_doze_result"
bigger_file = "novoresult_treze_result"

smaller_trie = mount_trie_from_bitfile(smaller_file)
bigger_trie = mount_trie_from_bitfile(bigger_file)
#temq ser na bigger

def obtain_v1_ext(v1, half_k, k):
        odd = k % 2
        m = 4**(half_k+ (1-odd))
        ext_v1s = bitarray.bitarray(m)
        #fazer um bitarray e ir incrementing
        for i in range(4):
            b = i
            extended_v1 = (b << (half_k*2)) | v1
            ext_v1s[extended_v1] = 1
        return ext_v1s

def search_trivial_test(trie):
    trivial = bigger_trie.find_trivial_ext(smaller_trie)
    count=bigger_trie.count_kmers()
    perc = (trivial/count) *100
    print(perc)
    #print(v1s_small)
    #print(ext_v1s)
    #print(found)


search_trivial_test(smaller_trie)
