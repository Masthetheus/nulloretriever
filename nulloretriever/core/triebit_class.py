"""Defines the TrieBit class for nullomer retrieving pipeline"""

from bitarray import bitarray
import struct

from nulloretriever.analysis.composition import (
    generate_gc_dict
)
from nulloretriever.analysis.motifs import (
    generate_cpg_dict,
    generate_complement_index_dict,
    generate_homopolymer_array,
    retrieve_nullomers_cpg_stats,
)

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
    def __init__(self, m, k, half_k):
        self.root = TrieBitNode()
        self.m = m
        self.k = k
        self.half_k = half_k

        byte_to_format = {1:'B', 2:'H', 4:'I', 8:'Q'}
        if half_k < 4:
            idx_sz, cnt_sz = 1, 1
        elif half_k == 4:
            idx_sz, cnt_sz = 1, 2
        elif half_k < 8:
            idx_sz, cnt_sz = 2, 2
        elif half_k == 8:
            idx_sz, cnt_sz = 2, 4
        else:
            idx_sz, cnt_sz = 4, 8

        self.format_code, self.counter_code = idx_sz, cnt_sz
        self.index_format, self.counter_format = byte_to_format[idx_sz], byte_to_format[cnt_sz]

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

    def insert_from_bit(self, v1, v2s):
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
                list(map(node.children[value].v2_set.__setitem__, v2s, [1] *
                         len(v2s)))
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

    def traverse_till_half_k(self, callback):
        def walk(node, depth, idx_acc):
            if depth == self.half_k:
                callback(node, idx_acc)
                return
            for child_value, child_node in enumerate(node.children):
                if child_node is not None:
                    next_idx = (idx_acc << 2) | child_value
                    walk(child_node, depth + 1, next_idx)
        walk(self.root, 0, 0)

    def traverse_till_custom(self, callback, target_idx,v2=None):
        def walk(node, depth, idx_acc):
            i = 1
            if depth == self.half_k:
                callback(node, idx_acc,v2)
                return
            else:
                next_idx = (target_idx >> (2*(self.half_k - depth - 1)) & 3)
                idx_acc = (idx_acc << 2)|next_idx 
                walk(node.children[next_idx], depth+1, idx_acc)
        walk(self.root, 0, 0)

    def count_gc(self):
        """Count percentage of GC of organism nullomers.
            Args:
                self(TrieBit): TrieBit object to be saved in compact binary
                output(str): path to save the file
            Returns:
                gc_percent(float): GC counted/number of nullomer bases
        """
        half_k = self.k//2
        gc_dict = generate_gc_dict(half_k)
        gc_tot = 0
        null_count = 0
        def count_gc(node, v1_idx):
            nonlocal gc_tot, null_count
            v2_idxs = node.v2_set.search(1)
            for v2 in v2_idxs:
                gc_tot += gc_dict[v2]
            v2_count = node.v2_set.count(bitarray('1'))
            null_count += v2_count
            v1_gc = 0
            for base in range(self.half_k):
                v1_gc += (v1_idx >> (2*1) & 1)
            v1_gc = gc_dict[v1_idx]
            gc_tot += (v1_gc*v2_count)
            return
        self.traverse_till_half_k(callback=count_gc)

        tot_bases = null_count * self.k
        if tot_bases > 0:
            gc_percent = (gc_tot / tot_bases) * 100
        else:
            gc_percent = 0
        return gc_percent

    def retrieve_nullomers_cpg_stats(self):
        half_k = self.k // 2
        null_count = 0
        cpg_dict = generate_cpg_dict(half_k)
        cpg_tot = 0
        null_with_cpg = 0
        v2_size = half_k
        v2_shift = (half_k - 1) * 2
        v1_cpg = 0
        def count_cpg(node, idx_acc):
            nonlocal v1_cpg, cpg_tot, null_count, null_with_cpg
            for i in range (self.half_k -1):
                shift = (self.half_k - i - 2) * 2
                pair = (idx_acc >> shift) & 15
                if pair == 6:
                    v1_cpg += 1
            v1_last = idx_acc & 3
            v2_idxs = node.v2_set.search(1)
            for v2 in v2_idxs:
                v2_first = v2 >> v2_shift
                has_cpg = False
                if cpg_dict[v2] != 0:
                    cpg_tot += cpg_dict[v2]
                    has_cpg = True
                if v1_last == 1 and v2_first == 3:
                    cpg_tot +=1
                    has_cpg = True
                if has_cpg:
                    null_with_cpg += 1
            v2_count = node.v2_set.count(bitarray('1'))
            null_count += v2_count
            cpg_tot += (v1_cpg * v2_count)
            return
        self.traverse_till_half_k(callback=count_cpg)
        try:
            cpg_count_mean = cpg_tot/null_with_cpg
        except ZeroDivisionError:
            cpg_count_mean = 0
        try:
            cpg_global_mean = (null_with_cpg/null_count)*100
        except ZeroDivisionError:
            cpg_global_mean = 0
        cpg_stats = {
            "total": cpg_tot,
            "nullomers_with_cpg": null_with_cpg,
            "global_mean": cpg_global_mean,
            "mean_nullomers_with_cpg": cpg_count_mean
        }
        return cpg_stats

    def retrieve_palindrome_stats(self):
        half_k = self.k//2
        palindrome_count = 0
        total_null = 0
        complement_index_dict = generate_complement_index_dict(half_k)
        def check_palindromy(node, idx_acc):
            nonlocal palindrome_count, total_null
            v1_adjusted = idx_acc >> 2
            v1_comp = complement_index_dict[v1_adjusted]
            total_null += node.v2_set.count(bitarray('1'))
            try:
                if node.v2_set[v1_comp]:
                    palindrome_count += 1
            except Exception as e:
                print(e)
                print(f"v2 set len {len(node.v2_set)}")
                print(f"k {self.k} and half {self.half_k}")
            return
        self.traverse_till_half_k(callback=check_palindromy)
        try:
            palindrome_relative = (palindrome_count/total_null) * 100
        except ZeroDivisionError:
            palindrome_relative = 0
        palindrome_stats = {
            "count": palindrome_count,
            "relative_fraction": palindrome_relative
        }
        return palindrome_stats

    def retrieve_homopolymer_stats(self):
        half_k = self.k//2
        found_homopolymers = []
        homopolymers = set(generate_homopolymer_array(half_k))
        homopolymers_v1 = set(generate_homopolymer_array(self.half_k))
        def gather_homopolymer(node, target_idx):
            nonlocal found_homopolymers
            if self.half_k != half_k:
                target_idx = target_idx >> 2
            if node.v2_set[target_idx]:
                found_homopolymers.append(target_idx)
            return
        for homopolymer in homopolymers_v1:
            self.traverse_till_custom(target_idx = homopolymer, callback = gather_homopolymer)
        return found_homopolymers

    def retrieve_v1_list(self):
        m = 4**self.half_k
        v1s = bitarray(m)
        def gather_v1s(node, idx_acc):
            nonlocal v1s
            v1s[idx_acc] = 1
            return
        self.traverse_till_half_k(callback=gather_v1s)
        return v1s

    def retrieve_v2_list(self, target_idx):
        v2s = bitarray(self.m)
        def gather_v2s(node,idx_acc):
            nonlocal v2s
            v2s = node.v2_set
            return
        self.traverse_till_custom(target_idx=target_idx, callback=gather_v2s)
        return v2s

    def retrieve_prime_null(self, second_trie):
        v1s = self.retrieve_v1_list()
        second_v1s = second_trie.retrieve_v1_list()
        common = v1s & second_v1s
        common_idxs = list(common.search(bitarray('1')))
        common_v2s = {}
        count = 0
        for v1 in common_idxs:
            v2s = self.retrieve_v2_list(v1)
            second_v2s = self.retrieve_v2_list(v1)
            common_v2s_array = v2s & second_v2s
            common_v2s_idxs = list(common_v2s_array.search(bitarray('1')))
            common_v2s[v1] = common_v2s_idxs
            count += 1
        return common_v2s

    def find_trivial_ext(self, small_trie):
        trivial = 0
        half_k = self.k//2
        smaller_k= self.k - 1
        is_odd = smaller_k%2
        smaller_hk = (smaller_k//2) + is_odd
        smaller_mask = (1 << ((smaller_hk-is_odd)*2)) - 1
        mask = (1 << ((self.k *2)-2)) - 1
        def is_in_trie(node, target_idx, v2):
            nonlocal trivial
            if node.v2_set[v2]:
                trivial += 1
            return
        def search_trivial(node, v1_idx):
            nonlocal trivial, smaller_mask
            ext_idxs = set()
            v2_idxs = node.v2_set.search(1)
            for v2 in v2_idxs:
                idx = (v1_idx << (self.half_k*2))|v2
                ext_idxs.add(idx & mask)
                ext_idxs.add(idx >> 2)
            for idx in ext_idxs:
                v1 = idx >> ((smaller_hk-is_odd)*2)
                v2 = idx & smaller_mask
                small_trie.traverse_till_custom(target_idx=v1, callback=is_in_trie, v2 = v2)

        self.traverse_till_half_k(callback=search_trivial)
        print(trivial)
        return trivial

    def find_root_v2(self):
        v2s = bitarray(self.m)
        def v2_minus_last(node,idx_acc):
            nonlocal v2s
            v2_idxs = list(node.v2_set.search(1))
            for v2 in v2_idxs:
                print(v2)
                v2_minus = v2 >> 2
                print(f"v2 {v2} and minud {v2_minus}")
                v2s[v2_minus] = 1
                return
        self.traverse_till_half_k(callback=v2_minus_last)
        return v2s
    def missing_path_idx(self, path):
        print(path)
        init_idx = sum(base*(4**(self.half_k - i - 1)) for i, base in
                       enumerate(path))
        abs_idx = 4**(self.half_k - len(path))
        print(f"abs: {abs_idx}")
        return range(init_idx, init_idx + abs_idx)

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
            f.write(struct.pack('<HHHBB', version, self.k, self.half_k, self.format_code, self.counter_code))
            def collect_nodes(node, path):
                if len(path) == self.half_k:
                    nullomers = node.v2_set.search(bitarray('0'))
                    null_count = node.v2_set.count(bitarray('0'))
                    if nullomers and null_count > 0:
                        buffer = bytearray()
                        index = sum(base * (4 ** (self.half_k - i - 1))
                                    for i, base in enumerate(path))
                        buffer.extend(struct.pack(f'{self.index_format}',
                                                  index))
                        buffer.extend(struct.pack(f'{self.counter_format}',
                                                  null_count))
                        for v2_index in node.v2_set.search(bitarray('0')):
                            buffer.extend(struct.pack(f'{self.index_format}',
                                          v2_index))
                        f.write(buffer)
                    return
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        collect_nodes(child_node, path + [child_value])
                    else:
                        for missing_idx in self.missing_path_idx(path+[child_value]):
                            buffer = bytearray()
                            buffer.extend(struct.pack(f'{self.index_format}',
                                                      missing_idx))
                            buffer.extend(struct.pack(f'{self.counter_format}',
                                                      0))
                            f.write(buffer)
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
