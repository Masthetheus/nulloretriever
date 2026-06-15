"""Functions for nullomer data processing aiming further analysis."""

import pandas as pd
import struct

from nulloretriever.core.triebit_class import TrieBit

def json_to_csv_mapping(df_json, df_csv):
    """Add all json columns to a csv df."""
    columns_to_add = list(df_json)
    print(columns_to_add)
    print(df_json.index())
    for column in columns_to_add:
        try:
            mapping = df_json.set_index('e')[column].replace(".","_")
            print("\n")
            print(mapping)
            column = column.lower().replace(" ", "_")
            df_csv[column] = df_csv['organism'].map(mapping)
        except Exception:
            print(f"Colum {column} Nao foi")
            continue
    return df_csv


def mount_trie_from_bitfile(filename):
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        k = struct.unpack('<H', l_bytes)[0]
        l_bytes = f.read(2)
        half_k = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        counter_size = struct.unpack('<B', f.read(1))[0]
        counter_byte_format = byte_to_format[counter_size]
        if k % 2 == 0:
            m = 4**half_k
            k_mask = (2**k) - 1
            v1_size = v2_size = half_k
        else:
            m = 4**(half_k-1)
            k_mask = (2**(k-1)) - 1
            v1_size = half_k
            v2_size = k - half_k
        trie = TrieBit(m, k, half_k)
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1= struct.unpack(f'<{byte_format}', index_bytes)[0]
                i = half_k - 1
                v1_bits = []
                while i >= 0:
                    v1_bits.append(v1 >> (i*2) & 3)
                    i -= 1
                nullomer_count_bytes = f.read(counter_size)
                if len(nullomer_count_bytes) < counter_size:
                    break
                nullomer_count = struct.unpack(f'<{counter_byte_format}', nullomer_count_bytes)[0]
                total_bytes = nullomer_count * byte_size
                v2s = []
                i = 0
                while i < nullomer_count:
                    nullomer_byte = f.read(byte_size)
                    total_bytes -= byte_size
                    v2s.append(struct.unpack(f'<{byte_format}',
                                             nullomer_byte)[0])
                    i += 1
                trie.insert_from_bit(tuple(v1_bits), v2s)
        except (struct.error, OSError):
            pass
    return trie
