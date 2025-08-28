"""Defines the TrieBit class for nullomer retrieving pipeline"""

from bitarray import bitarray
import struct

class TrieBitNode:
    def __init__(self, m):
        self.children = [None] * 4  # 0:A, 1:T, 2:C, 3:G
        self.v2_set = bitarray(m)
        self.v2_set.setall(0)

    def iterate(self, path=None):
        if path is None:
            path = []
        yield path, self
        for i, child in enumerate(self.children):
            if child is not None:
                yield from child.iterate(path + [i])

class TrieBit:
    def __init__(self, m, l):
        self.root = TrieBitNode(m)
        self.m = m
        self.l = l

        if m <= 256:
            self.index_format = 'B'
            self.format_code = 1
        elif m <= 65536:
            self.index_format = 'H' 
            self.format_code = 2
        elif m <= 4294967296:
            self.index_format = 'I'
            self.format_code = 4
        else:
            self.index_format = 'Q'
            self.format_code = 8

        def build(node, depth):
            if depth == l:
                return
            for i in range(4):
                if node.children[i] is None:
                    node.children[i] = TrieBitNode(m)
                build(node.children[i], depth + 1)
        build(self.root, 0)

    def insert(self, v1, v2):
        node = self.root
        for value in v1:
            if node.children[value] is None:
                node.children[value] = TrieBitNode(self.m)
            node = node.children[value]
        node.v2_set[v2] = 1

    def iterate(self):
        yield from self.root.iterate([])

    def count_nullomers(self, target_length=None):
        l = target_length or self.l
        def dfs(node, depth):
            if depth == l:
                return sum(1 for genome_id in range(len(node.v2_set)) 
                        if not node.v2_set[genome_id])
            
            return sum(dfs(child, depth + 1) 
                    for child in node.children 
                    if child is not None)
        return dfs(self.root, 0)
    
    def write_bit_format(self, output):
        """
            Saves TrieBit to a compact binary format.
            Format: [header][nodes...]
            Header: b'TRIE'[4] + version(2) + l(2) + format_code(1)
            Args:
                self: TrieBit object
                output(str): path to save the file
                target_leght(int): size of each halve of the nullomers.
            Returns:
                file: all nullomers sequences in binary format, where:
                    v1: index of the first half of the sequence, with size l (k/2)
                    v2_set size: total of v2 for the given v1
                    v2: index of the second half of the sequence, with size l (k/2), that are directly connected to the previous v1 value
        """
        with open(output, 'wb') as f:
            f.write(b'TRIE')  # Magic number
            version = 1
            f.write(struct.pack('<HHB', version, self.l, self.format_code))  # version, l, byte_size
            def collect_nodes(node, path):
                if len(path) == self.l:
                    nullomers = [i for i, bit in enumerate(node.v2_set) if not bit]
                    if nullomers:
                        # Compute lexicographic index for v1 path
                        index = sum(base * (4 ** (self.l - i - 1)) for i, base in enumerate(path))
                        f.write(struct.pack(f'<{self.index_format}', index))
                        f.write(struct.pack(f'<{self.index_format}', len(nullomers)))
                        for v2_index in nullomers:
                            f.write(struct.pack(f'<{self.index_format}', v2_index))
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
                if len(path) == self.l:
                    nullomers = [i for i, bit in enumerate(node.v2_set) if not bit]
                    if nullomers:
                        v1_index = sum(base * (4 ** (self.l - i - 1)) for i, base in enumerate(path))
                        f.write(f">{v1_index}\n")
                        v2_values = ",".join(str(i) for i in nullomers)
                        f.write(f"{v2_values}\n")
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        dfs(child_node, path + [child_value])
            dfs(self.root, [])

