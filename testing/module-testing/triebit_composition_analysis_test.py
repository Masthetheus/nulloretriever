from nulloretriever.analysis.composition import *

file = "../data/trie_binary_test"
print(f"Starting GC mean calculation test on the file {file}")
gc_mean = nullomers_gc_mean(file)
print(f"The tested organism has a GC% of {gc_mean}!")
print("GC mean functions module tested!")