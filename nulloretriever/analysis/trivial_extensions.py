"""Defines functions aimed at analysis of trivial nullomer extensions."""

import struct
from collections import defaultdict

from nulloretriever.analysis.parsings import obtain_trie_infos


def obtain_sequence_list_from_v1(v1, half_k):
    """Receives a v1 value and expand into it's singular value form."""
    v1_bits = []
    i = half_k - 1
    while i >= 0:
        v1_bits.append(v1 >> (i*2) & 3)
        i -= 1
    return v1_bits


def obtain_new_v1_from_v2(v1_bits, v2, half_k):
    """Given a v2 from a list obtain all possible related v1 for k+1."""
    indexv1 = 0
    j = half_k
    extended_v1 = v1_bits + [(v2 >> ((half_k-1)*2) & 3)]
    indexv1 = 0
    for base in extended_v1:
        indexv1 += base*(4**j)
        j -= 1
    return indexv1


def obtain_v2_extensions(v2, v2_extensions, half_k):
    """Given certain v2 value, returns it's possible 4 base extensions."""
    mask = (1 << ((half_k*2)-2)) - 1
    v2_diminished = v2 & mask
    for i in range(4):
        v2_extensions.append((v2_diminished << 2) | i)
    return v2_extensions


def compare_null_tries(file1, file2):
    """Search for trivial extension nullomers on file 1.

    Receives two files of nullomers from an organism where k'-k" = +-1. Then,
    searches the bigger k value file for nullomers that are trivial extensions
    of nullomers from the smaller k value file, that being one's that consists
    of a nullomer already found with one new base extension in either "side".
    Args:
        smallerkfile(str): path to the smaller k TrieBit bit file location
        biggerkfile(str): path to the bigger k TrieBit bit file location
    Returns:
        found(str): extended nullomers found
        count(int): total of extended nullomers found
    """

    return True


def v1_possible_extensions(v1, half_k):
    """Retrieves possible trivial extensions from a v1."""
    v1_extensions = []
    for i in range(4):
        v1_extensions.append(v1+(i*(4**half_k)))
    return v1_extensions


def assign_v2_to_v1(v2s, half_k):
    """Assign to each v2 extension it's related v1 value."""
    v2s_extensions = {0: [], 1: [], 2: [], 3: []}
    for v2 in v2s:
        v2s_extended = []
        obtain_v2_extensions(v2, v2s_extended, half_k)
        first_base = v2 >> (((half_k-1)*2)) & 3
        # print(f"v2: {v2}, first base: {first_base}")
        for value in v2s_extended:
            v2s_extensions[first_base].append(value)
    return v2s_extensions


def retrieve_trivial_extensions(v1, v2s, v1_indexes_file2, file2, half_k,
                                trivial_extensions, odd):
    """Count the number of trivial nullomer extensions in a bit file."""
    if odd:
        v1_adjusted = v1 >> 2
        v1_extensions = v1_possible_extensions(v1_adjusted, half_k)
    else:
        v1_extensions = v1_possible_extensions(v1, half_k)
    orig_v2_extensions = assign_v2_to_v1(v2s, half_k)
    v2s_bigger_k = set()
    byte_size, byte_format = obtain_trie_infos(file2)
    new_counter = 0
    for v1 in v1_extensions:
        counter = 0
        v1_first = v1 >> (((half_k-1)*2) - 2) & 3
        with open(file2, 'rb') as f:
            f.seek(v1_indexes_file2[v1], 0)
            try:
                nullomer_count_bytes = f.read(byte_size)
                nullomer_count = struct.unpack(
                    f'<{byte_format}', nullomer_count_bytes)[0]
                while counter < nullomer_count:
                    current_v2_bytes = f.read(byte_size)
                    v2s_bigger_k.add(struct.unpack(
                        f'<{byte_format}', current_v2_bytes)[0])
                    counter += 1
                set_v2s = set(orig_v2_extensions[v1_first])
                set_v2s_orig = set(v2s)
                extensions1 = v2s_bigger_k & set_v2s
                extensions2 = v2s_bigger_k & set_v2s_orig
                if not extensions1 and not extensions2:
                    new_counter = 0
                    continue
                else:
                    new_counter += len(extensions1) + len(extensions2)
                    trivial_extensions[v1].update(extensions1)
                    trivial_extensions[v1].update(extensions2)
            except (struct.error, OSError):
                pass
    return new_counter
