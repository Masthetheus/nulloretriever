from nulloretriever.core.triebit_class import *

k = input("Which k value to create a Trie for?")
if not k:
    k = 10
else:
    k = int(k)
l = int(k/2)
m = 4**l
print(f"Creating TrieBit for m = {m} and l = {l}")
trie = TrieBit(m,l)
print(f"Checking if all nodes are inserted via nullomer count. Result must be = {4**k}.")
count = trie.count_nullomers()
print(count)
indexes = []
i = 0
while i < m:
    indexes.append(i)
    i += 23
print(f"Index list generated for testing:\n {indexes}")
v1 = []
v1_counter = 0
while v1_counter < l:
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
