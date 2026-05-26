"""Defines the TrieBit class for nullomer retrieving pipeline"""

from bitarray import bitarray
import struct


class TrieBitNode:
    """Class for non-terminal nodes.
    Regular node of a Bit Trie class, representing one of the four nucleotides.
    Each node is composed of up to four referenced children nodes.
    """
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
    """Class for terminal nodes.
    Leaf node of a Bit Trie class, representing the last base of given half_k
    path.
    It's main difference consists on it's "children", that being an bitarray
    that represents the possible sequences indexes for given k value.
    All bits start at 0, representing absence of such sequences.
    """
    def __init__(self, m):
        self.v2_set = bitarray(m)
        self.v2_set.setall(0)

class TrieBit:
    """Main class for nullomer optimized search and retrieval.
    Composite class of multiple BitNodes objects with bit arrays attached to
    terminal BitLeaf nodes at depth half_k or half_k + 1.
    Lessen the impact of redundancy in the search by splitting each sequence in
    two main indexes, the first one for the main path and the second one being
    represented via bit setting in the respective bit array.
    Args:
        m(int): Number of possible indexes for the second split half. Obtained
        in it's general form by 4**(half_k).
        half_k(int): Half of the k value. For odd k values, half_k must be
        (k/2)+1, since, by own convention, the first half balances the division.
    Return:
        triebit(trie)S: An initialized TrieBit object with half_k depth.
    """
    def __init__(self, m, half_k):
        self.root = TrieBitNode()
        self.m = m
        self.half_k = half_k

        # m-1 aiming to adjust to the 0 index, since m = 4**half_k, m will
        # represent the total possible indexes, not accounting 0. 
        if m-1 <= 255:
            self.index_format = 'B'
            self.format_code = 1
        elif m-1 <= 65535:
            self.index_format = 'H'
            self.format_code = 2
        elif m-1 <= 4294967295:
            self.index_format = 'I'
            self.format_code = 4
        else:
            self.index_format = 'Q'
            self.format_code = 8

    def insert(self, v1, v2):
        node = self.root
        count = 0
        for value in v1:
            count += 1
            if node.children[value] is None:
                if count != self.half_k:
                    node.children[value] = TrieBitNode()
                else:
                    node.children[value] = TrieBitLeaf(self.m)
            if count == self.half_k:
                node.children[value].v2_set[v2] = 1
            else:
                node = node.children[value]

    def iterate(self):
        yield from self.root.iterate([])

    def count_kmers(self, target_length=None):
        half_k = target_length or self.half_k

        def dfs(node, depth):
            counter = 0
            leafs = 0
            if depth == half_k:
                return node.v2_set.count(1)

            return sum(dfs(child, depth + 1)
                       for child in node.children
                       if child is not None)
        return dfs(self.root, 0)

    def write_bit_format(self, output):
        """Saves TrieBit to a compact binary format.
            Format: [header][nodes...]
            Header: b'TRIE'[4] + version(2) + half_k(2) + format_code(1)
            Args:
                self(TrieBit): TrieBit object to be saved in compact binary
                output(str): path to save the file
            Returns:
                file: all nullomers sequences in binary format, where:
                    v1(int): index of the first half of the sequence, with size
                    half_k (k/2).
                    v2_set size(int): count of v2 for the given v1
                    v2(array): indexes of the second half of the sequence,
                        calculated with half_k (k/2), that are directly connected
                        to the previous v1 value
        """
        with open(output, 'wb') as f:
            f.write(b'TRIE')  # Magic number
            version = 1
            # version, half_k, byte_size
            f.write(struct.pack('<HHB', version, self.half_k, self.format_code))

            def collect_nodes(node, path):
                if len(path) == self.half_k:
                    nullomers = node.v2_set.search(bitarray('0'))
                    null_count = node.v2_set.count()
                    if nullomers:
                        buffer = bytearray()
                        index = sum(base * (4 ** (self.half_k - i - 1))
                                    for i, base in enumerate(path))
                        buffer.extend(struct.pack(f'{self.index_format}',
                                                  index))
                        buffer.extend(struct.pack(f'{self.index_format}',
                                                  null_count))
                        for v2_index in nullomers:
                            buffer.extend(struct.pack(f'{self.index_format}',
                                          v2_index))
                        f.write(buffer)
                    return
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        collect_nodes(child_node, path + [child_value])
            collect_nodes(self.root, [])

    def write_txt_format(self, output):
        """Writes a trie paths and relative v2 values in a compact txt format
            Args:
                self(TrieBit): TrieBit object to be saved in compact binary
                output(str): path to save the file
            Returns:
                file: compact .txt file as below:
                    >(char): v1 delimiter, for further automation of file
                    reading and processing
                    v1_index(int)
                    v2_values(array): comma separated v2 index values for the
                                    previous v1
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
                    return
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        dfs(child_node, path + [child_value])
            dfs(self.root, [])
