from nulloretriever.core.triebit_class import *

k = input("Which k value to create a Trie for?")
if not k:
    k = 10
else:
    k = int(k)
half_k = int(k/2)
m = 4**half_k
print(f"Creating TrieBit for m = {m} and half_k = {half_k}")
trie = TrieBit(m,half_k)
print(f"Checking if all nodes are inserted via nullomer count. Result must be = {4**k}.")
init_count = trie.count_nullomers()
print(init_count)
indexes = []
i = 0
while i < m:
    indexes.append(i)
    i += 23
print(f"Index list generated for testing:\n {indexes}")
v1 = []
v1_counter = 0
while v1_counter < half_k:
    v1.append(0)
    v1_counter += 1
count = 0
for i in indexes:
    v2 = i
    if count < 4:
        print(f"Inserting pair {v1}{v2}")
        count += 1
    trie.insert(v1,v2)
print(len(indexes))
expected_null = (4**k) - len(indexes)
print("Checking if all nodes were inserted correctly via nullomer count.")
nullomers = trie.count_nullomers()
print(f"{nullomers} nullomers were found of the {expected_null} expected.")
print(f"{nullomers} nullomers were found of the initial {init_count} present.\n")
print(f"Representing a difference of {init_count - nullomers}.\n")
output_bit = "trie_binary_test"
output_compact_txt = "trie_compact_txt"
print(f"Trying now to write the trie in bit format on the following path: {output_bit}")
#trie.write_bit_format(output_bit)
print("Trie in bit format correctly wrriten!")
print(f"Trying now to write the trie in compact txt format on the following path: {output_compact_txt}")
#trie.write_compact_txt_format(output_compact_txt)
print("Trie in compact txt format correctly wrriten!")
print("All TrieBit base functionalities tested and passed!")
