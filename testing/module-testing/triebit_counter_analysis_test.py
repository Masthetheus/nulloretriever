from nulloretriever.analysis.counter import *

file = "data/trie_binary_test"
print(f"Starting the counter test on the file {file}")
count = quick_nullomer_count(file)
print(f"Counted {count} nullomers on {file}.")
print(f"Counter working!")
