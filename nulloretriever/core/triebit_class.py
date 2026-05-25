"""Defines the TrieBit class for nullomer retrieving pipeline"""

from bitarray import bitarray
import struct


class TrieBitNode:
    def __init__(self):
        self.children = [None] * 4  

    def iterate(self, path=None):
        if path is None:
            path = []
        yield path, self
        for i, child in enumerate(self.children):
            if child is not None:
                yield from child.iterate(path + [i])

class TrieBitLeaf:
    def __init__(self, m):
        self.v2_set = bitarray(m)
        self.v2_set.setall(0)

class TrieBit:
    def __init__(self, m, half_k):
        self.root = TrieBitNode()
        self.m = m
        self.half_k = half_k

        if m <= 255:
            self.index_format = 'B'
            self.format_code = 1
        elif m <= 65535:
            self.index_format = 'H'
            self.format_code = 2
        elif m <= 4294967295:
            self.index_format = 'I'
            self.format_code = 4
        else:
            self.index_format = 'Q'
            self.format_code = 8

        def build(node, depth):
            print(f"Depth atual {depth}")
            if depth == half_k:
                for i in range(4):
                    print(f"DEBUG: LOOP {i}")
                    if node.children[i] is None:
                        print(f"Leaf adicionada no idx {i} na depth {depth}")
                        node.children[i] = TrieBitLeaf(self.m)
            for i in range(4):
                print("CHILDRENS")
                print(node.children)
                if node.children[i] is None:
                    node.children[i] = TrieBitNode()
                    print(f"Criado node de idx {i} na depth {depth}")
                    build(node.children[i], depth + 1)
                print(f"DEBUG: terminando loop maior de idx {i}")
        build(self.root, 0)

    def insert(self, v1, v2):
        node = self.root
        for value in v1:
            node = node.children[value]
        node.v2_set[v2] = 1

    def iterate(self):
        yield from self.root.iterate([])

    def count_nullomers(self, target_length=None):
        half_k = target_length or self.half_k

        def dfs(node, depth):
            counter = 0
            if depth == half_k:
                print(f"Returning {node.v2_set.count(0)} null counted .")
                return node.v2_set.count(0)

            return sum(dfs(child, depth + 1)
                       for child in node.children
                       if child is not None)
        return dfs(self.root, 0)

    def write_bit_format(self, output):
        """
            Saves TrieBit to a compact binary format.
            Format: [header][nodes...]
            Header: b'TRIE'[4] + version(2) + half_k(2) + format_code(1)
            Args:
                self: TrieBit object
                output(str): path to save the file
                target_leght(int): size of each halve of the nullomers.
            Returns:
                file: all nullomers sequences in binary format, where:
                    v1: index of the first half of the sequence, with size half_k (k/2)
                    v2_set size: total of v2 for the given v1
                    v2: index of the second half of the sequence, with size half_k (k/2), that are directly connected to the previous v1 value
        """
        with open(output, 'wb') as f:
            f.write(b'TRIE')  # Magic number
            version = 1
            # version, half_k, byte_size
            f.write(struct.pack('<HHB', version, self.half_k, self.format_code))

            def collect_nodes(node, path):
                if len(path) == self.half_k:
                    nullomers = [i for i, bit in enumerate(
                        node.v2_set) if not bit]
                    if nullomers:
                        # Compute lexicographic index for v1 path
                        index = sum(base * (4 ** (self.half_k - i - 1))
                                    for i, base in enumerate(path))
                        f.write(struct.pack(f'<{self.index_format}', index))
                        f.write(struct.pack(
                            f'<{self.index_format}', len(nullomers)))
                        for v2_index in nullomers:
                            f.write(struct.pack(
                                f'<{self.index_format}', v2_index))
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        collect_nodes(child_node, path + [child_value])
            collect_nodes(self.root, [])

    def write_compact_txt_format(self, output):
        """Writes a trie paths and relative v2 values in a compact txt format
        Args:
            self: TrieBit object
            output(str): path where the compact txt shall be stored
        Returns:
            file: compact .txt file as below:
                >v1_index(int)
                v2_values(arr): all v2 index values for the previous v1 value
        """
        with open(output, 'w') as f:
            def dfs(node, path):
                if len(path) == self.half_k:
                    nullomers = [i for i, bit in enumerate(
                        node.v2_set) if not bit]
                    if nullomers:
                        v1_index = sum(base * (4 ** (self.half_k - i - 1))
                                       for i, base in enumerate(path))
                        f.write(f">{v1_index}\n")
                        v2_values = ",".join(str(i) for i in nullomers)
                        f.write(f"{v2_values}\n")
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        dfs(child_node, path + [child_value])
            dfs(self.root, [])
