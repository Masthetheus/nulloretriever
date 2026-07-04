"""Defines function that parses files aiming downstream analysis."""

import struct
from collections import defaultdict


def obtain_trie_infos(filename):
    """Opens an trie object and obtain it's byte format."""
    byte_to_format = {1: "B", 2: "H", 4: "I", 8: "Q"}
    with open(filename, "rb") as f:
        # Skip header
        f.seek(8)  # Skip magic(4) + version(2) + half_k_val(2)
        byte_size = struct.unpack("<B", f.read(1))[0]
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
    with open(filename, "rb") as f:
        # Skip header
        f.seek(9)  # Skip magic(4) + version(2) + half_k_val(2)
        if skip_counter > 0:
            f.seek(skip_counter)
        try:
            bytes_v1 = f.read(byte_size)
            v1 = struct.unpack(f"<{byte_format}", bytes_v1)[0]
            nullomer_count_bytes = f.read(byte_size)
            nullomer_count = struct.unpack(f"<{byte_format}", nullomer_count_bytes)[0]
            while counter < nullomer_count:
                v2_bytes = f.read(byte_size)
                v2s.append(struct.unpack(f"<{byte_format}", v2_bytes)[0])
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
            v1_extensions.append(v1 + (i * (4**half_k)))
        return v1_extensions
    if v2:
        half_k


def obtain_v1_positions(file):
    """Retrieves the location of the v1 indexes inside a bit trie file."""
    byte_to_format = {1: "B", 2: "H", 4: "I", 8: "Q"}
    v1_locations = defaultdict(set)
    with open(file, "rb") as f:
        # Skip header
        f.seek(8)
        byte_size = struct.unpack("<B", f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        try:
            while True:
                bytes_v1 = f.read(byte_size)
                if len(bytes_v1) < byte_size:
                    break
                v1_location = f.tell()
                v1 = struct.unpack(f"<{byte_format}", bytes_v1)[0]
                v1_locations[v1] = v1_location
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f"<{byte_format}", nullomer_count_bytes)[
                    0
                ]
                f.seek(nullomer_count * byte_size, 1)
        except (struct.error, OSError) as e:
            print(e)
            pass
    return v1_locations


def gather_v1_related_data(v1_bit_index, file):
    """Gather informations from given v1 value from a bit nullomer file."""
    byte_to_format = {1: "B", 2: "H", 4: "I", 8: "Q"}
    with open(file, "rb") as f:
        f.seek(8)
        byte_size = struct.unpack("<B", f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        f.seek(v1_bit_index, 0)
        counter = 0
        v2s = set()
        try:
            nullomer_count_bytes = f.read(byte_size)
            nullomer_count = struct.unpack(f"<{byte_format}", nullomer_count_bytes)[0]
            while counter < nullomer_count:
                v2_bytes = f.read(byte_size)
                v2s.add(struct.unpack(f"<{byte_format}", v2_bytes)[0])
                counter += 1
        except (struct.error, OSError):
            pass
    return v2s


def obtain_nullomer_set(file, v1_positions):
    """Retrieves the collection of v1 and v2 from a nullomer bit file."""
    byte_size, byte_format = obtain_trie_infos(file)
    v1_count = 0
    nullomer_set = defaultdict(set)
    with open(file, "rb") as f:
        for v1, position in v1_positions.items():
            f.seek(position, 0)
            counter = 0
            try:
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f"<{byte_format}", nullomer_count_bytes)[
                    0
                ]
                while counter < nullomer_count:
                    v2_bytes = f.read(byte_size)
                    nullomer_set[v1].add(struct.unpack(f"<{byte_format}", v2_bytes)[0])
                    counter += 1
            except (struct.error, OSError) as e:
                print(e)
                pass
            v1_count += 1
    return nullomer_set
