"""Defines function that parses files aiming downstream analysis."""

import struct


def obtain_trie_infos(filename):
    """Opens an trie object and obtain it's byte format."""
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(8)  # Skip magic(4) + version(2) + half_k_val(2)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]

    return byte_size, byte_format


def read_per_v1(filename, skip_counter, total_null_count):
    """Estimates and counts trivial nullomeric extensions of k+1.

    Given an triebit file filename, gathers only one of it's v1 indexes and
    it's correspondent v2 values. It's intent is to lower RAM usage during
    comparisons between two or more bit files.
    Args:
        filename(str): path to the TrieBit bit file location
        skip_counter(int): number of bytes to skip in the file, to get to the
            wanted v1, allowing also iteration. If equals 0, no bytes are
            'jumped'.
        total_null_count(int): Counts the amount of null obtained after all
            iterations. Used as metrics for iterations.
    Returns:
        v1(int): v1 index gathered.
        v2s(array): Correspondent v2 values to the v1 gathered.
        counter(int): Count of how many nullomers were extracted.
        total_null_count(int): Cumulative nullomer counter, for loop purposes.
        skip_counter(int): .tell() value, to skip v1 values as last stop.
    """
    byte_size, byte_format = obtain_trie_infos(filename)
    v2s = []
    v1 = 0
    counter = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(9)  # Skip magic(4) + version(2) + half_k_val(2)
        if skip_counter > 0:
            f.seek(skip_counter)
        try:
            bytes_v1 = f.read(byte_size)
            v1 = struct.unpack(f'<{byte_format}', bytes_v1)[0]
            nullomer_count_bytes = f.read(byte_size)
            nullomer_count = struct.unpack(
                f'<{byte_format}', nullomer_count_bytes)[0]
            while counter < nullomer_count:
                v2_bytes = f.read(byte_size)
                v2s.append(struct.unpack(f'<{byte_format}', v2_bytes)[0])
                counter += 1
            total_null_count += counter
            skip_counter = f.tell()

        except (struct.error, OSError):
            pass
    return v1, v2s, counter, total_null_count, skip_counter


def sequence_length_conversion(half_k, v1=False, v2=False):
    """Transform a index from size k to k-1."""
    if v1:
        v1_extensions = []
        for i in range(4):
            print(i)
            v1_extensions.append(v1+(i*(4**half_k)))
        return v1_extensions
    if v2:
        half_k
