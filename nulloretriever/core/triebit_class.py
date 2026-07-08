"""Defines the TrieBit class for nullomer retrieving pipeline"""

from bitarray import bitarray
import array
import struct

from nulloretriever.analysis.composition import generate_gc_dict
from nulloretriever.analysis.motifs import (
    generate_cpg_dict,
    generate_complement_index_dict,
    generate_homopolymer_array,
)


class TrieBitNode:
    """Class for non-terminal nodes.
    Regular node of a Bit Trie class, representing one of the four nucleotides.
    Each node is composed of up to four referenced children nodes.
    """

    __slots__ = ('children',)

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

    __slots__ = ('v2_set',)

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

        byte_to_format = {1: "B", 2: "H", 4: "I", 8: "Q"}
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
        self.index_format, self.counter_format = (
            byte_to_format[idx_sz],
            byte_to_format[cnt_sz],
        )

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

    def insert(self, v1, v2):
        node = self.root
        half_k = self.half_k
        for shift in range(half_k - 1, 0, -1):
            base = (v1 >> (shift * 2)) & 3
            child = node.children[base]
            if child is None:
                child = TrieBitNode()
                node.children[base] = child
            node = child
        base = v1 & 3
        leaf = node.children[base]
        if leaf is None:
            leaf = TrieBitLeaf(self.m)
            node.children[base] = leaf
        leaf.v2_set[v2] = 1

    def iterate(self):
        yield from self.root.iterate([])

    def count_kmers(self, target_length=None):
        half_k = target_length or self.half_k

        def dfs(node, depth):
            if depth == half_k:
                return node.v2_set.count(1)

            return sum(
                dfs(child, depth + 1) for child in node.children if child is not None
            )

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

    def traverse_till_custom(self, callback, target_idx, v2=None):
        def walk(node, depth, idx_acc):
            if depth == self.half_k:
                if v2:
                    return callback(node, idx_acc, v2)
                else:
                    return callback(node, idx_acc)
            else:
                next_idx = target_idx >> (2 * (self.half_k - depth - 1)) & 3
                idx_acc = (idx_acc << 2) | next_idx
                if node.children[next_idx]:
                    return walk(node.children[next_idx], depth + 1, idx_acc)

        return walk(self.root, 0, 0)

    def count_gc(self):
        """Count percentage of GC of organism nullomers.
        Args:
            self(TrieBit): TrieBit object to be saved in compact binary
            output(str): path to save the file
        Returns:
            gc_percent(float): GC counted/number of nullomer bases
        """
        gc_dict = generate_gc_dict(self.half_k)
        gc_tot = 0
        null_count = 0

        def count_gc(node, v1_idx):
            nonlocal gc_tot, null_count
            v2_idxs = node.v2_set.search(1)
            for v2 in v2_idxs:
                gc_tot += gc_dict[v2]
            v2_count = node.v2_set.count(bitarray("1"))
            null_count += v2_count
            v1_gc = 0
            for _ in range(self.half_k):
                v1_gc += v1_idx >> (2 * 1) & 1
            gc_tot += v1_gc * v2_count
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
        cpg_dict = generate_cpg_dict(half_k)
        cpg_tot = 0
        null_with_cpg = 0
        v2_shift = (half_k - 1) * 2

        def count_cpg(node, idx_acc):
            nonlocal cpg_tot, null_with_cpg
            v1_cpg = 0
            for i in range(self.half_k - 1):
                shift = (self.half_k - i - 2) * 2
                pair = (idx_acc >> shift) & 15
                if pair == 6:
                    v1_cpg += 1
            v1_last = idx_acc & 3
            v2_idxs = node.v2_set.search(1)
            v2_count = 0
            for v2 in v2_idxs:
                v2_first = v2 >> v2_shift
                has_cpg = False
                if cpg_dict[v2] != 0:
                    cpg_tot += cpg_dict[v2]
                    has_cpg = True
                    v2_count += 1
                if v1_last == 1 and v2_first == 3:
                    cpg_tot += 1
                    has_cpg = True
                    if cpg_dict[v2] == 0 and v1_cpg == 0:
                        v2_count += 1
                if has_cpg:
                    null_with_cpg += 1
            cpg_tot += v1_cpg * v2_count
            return

        self.traverse_till_half_k(callback=count_cpg)
        null_count = self.count_kmers()
        try:
            cpg_count_mean = cpg_tot / null_with_cpg
        except ZeroDivisionError:
            cpg_count_mean = 0
        try:
            cpg_global_mean = (null_with_cpg / null_count) * 100
        except ZeroDivisionError:
            cpg_global_mean = 0
        cpg_stats = {
            "total": cpg_tot,
            "nullomers_with_cpg": null_with_cpg,
            "global_mean": cpg_global_mean,
            "mean_nullomers_with_cpg": cpg_count_mean,
        }
        return cpg_stats

    def retrieve_palindrome_stats(self):
        half_k = self.k // 2
        palindrome_count = 0
        total_null = 0
        complement_index_dict = generate_complement_index_dict(half_k)

        def check_palindromy(node, idx_acc):
            nonlocal palindrome_count, total_null
            v1_adjusted = idx_acc >> 2
            v1_comp = complement_index_dict[v1_adjusted]
            total_null += node.v2_set.count(bitarray("1"))
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
            palindrome_relative = (palindrome_count / total_null) * 100
        except ZeroDivisionError:
            palindrome_relative = 0
        palindrome_stats = {
            "count": palindrome_count,
            "relative_fraction": palindrome_relative,
        }
        return palindrome_stats

    def retrieve_homopolymer_stats(self):
        half_k = self.k // 2
        found_homopolymers = []
        homopolymers_v1 = set(generate_homopolymer_array(self.half_k))

        def gather_homopolymer(node, target_idx):
            nonlocal found_homopolymers, half_k
            if self.half_k != half_k:
                target_idx = target_idx >> 2
            if node.v2_set[target_idx]:
                found_homopolymers.append(target_idx)
            return

        for homopolymer in homopolymers_v1:
            self.traverse_till_custom(
                target_idx=homopolymer, callback=gather_homopolymer
            )
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

        def gather_v2s(node, idx_acc, v2=None):
            nonlocal v2s
            v2s = node.v2_set
            return

        self.traverse_till_custom(target_idx=target_idx, callback=gather_v2s)
        return v2s

    def retrieve_prime_null(self, second_trie):
        v1s = self.retrieve_v1_list()
        second_v1s = second_trie.retrieve_v1_list()
        common = v1s & second_v1s
        v1_only_self = v1s & ~second_v1s
        common_idxs = list(common.search(bitarray("1")))

        def change_v2_set(node, idx_acc, v2=None):
            nonlocal second_v2s
            node.v2_set = second_v2s & node.v2_set
            return

        def zero_v2_set(node, idx_acc, v2=None):
            node.v2_set = bitarray(len(node.v2_set))
            return

        for v1 in common_idxs:
            second_v2s = second_trie.retrieve_v2_list(v1)
            self.traverse_till_custom(target_idx=v1, callback=change_v2_set)
        v1_only_self_idxs = list(v1_only_self.search(bitarray("1")))
        for v1 in v1_only_self_idxs:
            self.traverse_till_custom(target_idx=v1, callback=zero_v2_set)
        return

    def find_trivial_ext(self, small_trie):
        trivial = 0
        found = 0
        smaller_k = self.k - 1
        is_odd = smaller_k % 2
        smaller_hk = (smaller_k // 2) + is_odd
        smaller_mask = (1 << ((smaller_hk - is_odd) * 2)) - 1
        mask = (1 << ((self.k * 2) - 2)) - 1

        def is_in_trie(node, target_idx, v2):
            return node.v2_set[v2]

        def search_trivial(node, v1_idx):
            nonlocal trivial, smaller_mask, found
            v2_idxs = node.v2_set.search(1)
            for v2 in v2_idxs:
                idx = (v1_idx << (self.half_k * 2)) | v2
                for possibility in (idx >> 2, idx & mask):
                    v1 = possibility >> ((smaller_hk - is_odd) * 2)
                    v2_loop = possibility & smaller_mask
                    found = small_trie.traverse_till_custom(
                        target_idx=v1, callback=is_in_trie, v2=v2_loop
                    )
                    if found:
                        trivial += 1
                        break

        self.traverse_till_half_k(callback=search_trivial)
        return trivial

    def find_root_v2(self):
        v2s = bitarray(self.m)

        def v2_minus_last(node, idx_acc):
            nonlocal v2s
            v2_idxs = list(node.v2_set.search(1))
            for v2 in v2_idxs:
                v2_minus = v2 >> 2
                v2s[v2_minus] = 1
                return

        self.traverse_till_half_k(callback=v2_minus_last)
        return v2s

    def missing_path_idx(self, path):
        print(path)
        init_idx = sum(
            base * (4 ** (self.half_k - i - 1)) for i, base in enumerate(path)
        )
        abs_idx = 4 ** (self.half_k - len(path))
        print(f"abs: {abs_idx}")
        return range(init_idx, init_idx + abs_idx)

    def write_bit_format(self, output):
        with open(output, "wb") as f:
            f.write(b"TRIE")
            version = 1
            f.write(struct.pack("<HHHBB", version, self.k, self.half_k,
                               self.format_code, self.counter_code))

            # Map your format codes to array types
            type_map = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
            arr_type = type_map[self.format_code]
            cnt_type = type_map[self.counter_code]

            def collect_nodes(node, path):
                if len(path) == self.half_k:
                    nullomers = node.v2_set.search(bitarray("0"))
                    if nullomers:
                        null_count = len(nullomers)
                        index = sum(base * (4 ** (self.half_k - i - 1))
                                   for i, base in enumerate(path))

                        # Use array.array for fast C-level conversion
                        result = array.array(arr_type, [index])
                        count_arr = array.array(cnt_type, [null_count])
                        v2_arr = array.array(arr_type, nullomers)

                        f.write(result.tobytes())
                        f.write(count_arr.tobytes())
                        f.write(v2_arr.tobytes())
                    return

                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        collect_nodes(child_node, path + [child_value])
                    else:
                        for missing_idx in self.missing_path_idx(path + [child_value]):
                            result = array.array(arr_type, [missing_idx])
                            count_arr = array.array(cnt_type, [0])
                            f.write(result.tobytes())
                            f.write(count_arr.tobytes())

    def idx_to_seq(self, idx, k):
        DECODE = {0: "A", 1: "C", 2: "T", 3: "G"}
        bases = []
        for _ in range(k):
            bases.append(DECODE[idx & 3])
            idx >>= 2
        return "".join(reversed(bases))

    def write_sequences(self, filepath):
        results = []

        def collect(node, v1_idx):
            for v2 in node.v2_set.search(1):
                seq_v1 = self.idx_to_seq(v1_idx, self.half_k)
                seq_v2 = self.idx_to_seq(v2, self.k - self.half_k)
                results.append(seq_v1 + seq_v2)

        self.traverse_till_half_k(callback=collect)
        with open(filepath, "w") as f:
            f.write("\n".join(results))

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
        with open(output, "w") as f:

            def dfs(node, path):
                if len(path) == self.half_k:
                    nullomers = [i for i, bit in enumerate(node.v2_set) if bit]
                    if nullomers:
                        v1_index = sum(
                            base * (4 ** (self.half_k - i - 1))
                            for i, base in enumerate(path)
                        )
                        f.write(f">{v1_index}\n")
                        v2_values = ",".join(str(i) for i in nullomers)
                        f.write(f"{v2_values}\n")
                    return
                for child_value, child_node in enumerate(node.children):
                    if child_node is not None:
                        dfs(child_node, path + [child_value])

            dfs(self.root, [])
