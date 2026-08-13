"""Functions for nullomer data encoding to bit format."""

import struct
import array

def generate_sub_indexes(idx, k):

    is_odd = k%2
    half_k = (k//2) + is_odd
    half_k_v2 = k//2
    m = 4**(k//2)

    v1 = (idx >> (half_k_v2*2)) & ((1 << ((2*half_k))) - 1)
    v2 = idx & ((1 << (half_k_v2 * 2)) - 1)

    return v1, v2

def encode_to_bit(filename, k):

    nullomers = {}

    ENCODING = {"A": 0, "C": 1, "T": 2, "G": 3}

    with open(filename,"r") as f:
        for line in f:
            sequence = line.rstrip("\n")
            count = 0
            idx = 0
            while count < k:
                idx = (idx << 2) | ENCODING[sequence[count]]
                count += 1

            v1,v2 = generate_sub_indexes(idx,k)
            nullomers.setdefault(v1, []).append(v2)

    return nullomers

def write_bit_file(filename, output, k):

    is_odd = k%2
    half_k = (k//2) + is_odd
    m = 4**(k//2)

    nullomers = encode_to_bit(filename, k)
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

    format_code, counter_code = idx_sz, cnt_sz
    index_format, counter_format = (
        byte_to_format[idx_sz],
        byte_to_format[cnt_sz],
    )

    with open(output, "wb") as f:
        f.write(b"TRIE")
        version = 1
        f.write(struct.pack("<HHHBB", version, k, half_k, format_code, counter_code))
        type_map = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        arr_type = type_map[format_code]
        cnt_type = type_map[counter_code]

        for v1, v2s in nullomers.items():

            result = array.array(arr_type, [v1])
            count_arr = array.array(cnt_type,[len(v2s)])
            v2_arr = array.array(arr_type, v2s)

            f.write(result.tobytes())
            f.write(count_arr.tobytes())
            f.write(v2_arr.tobytes())

    return












