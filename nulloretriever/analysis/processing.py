"""Functions for nullomer data processing aiming further analysis."""

import pandas as pd
import struct

from nulloretriever.core.triebit_class import TrieBit

def json_to_csv_mapping(df_json, df_csv):
    """Add all json columns to a csv df."""
    df_json.index = df_json.index.str.replace('.', '_', regex = False)
    print(df_json)
    df_json.columns = df_json.columns.str.lower()
    df_json.columns = df_json.columns.str.replace(' ', '_')
    columns_to_add = list(df_json)
    print(f"COLUMNS TO ADD:")
    print(columns_to_add)
    print(f"COLUMNS JSON:")
    print(df_json.columns)
    print(f"COLUMNS CSV:")
    print(df_csv.columns)
    for column in columns_to_add:
        print(column)
        print("========== INICIO ==========")
        print(df_json.columns)
        print(df_csv.columns)
        try:
            mapping = df_json[column]
            print("\n")
            df_csv[column] = df_csv['organism'].map(mapping)
            print(df_csv[column])
        except Exception as e:
            print(f"Column {column} Nao foi")
            print(f"Razão {e} e motivo {type(e).__name__}")
            print('\n')
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
                shift_ranges = list(range(half_k - 1, -1, -1))
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                v1_bits = tuple((v1 >> (i * 2)) & 3 for i in shift_ranges)
                i = half_k - 1
                nullomer_count_bytes = f.read(counter_size)
                if len(nullomer_count_bytes) < counter_size:
                    break
                nullomer_count = struct.unpack(f'<{counter_byte_format}', nullomer_count_bytes)[0]
                if nullomer_count == 0:
                    v2s=list(range(m))
                else:
                    total_bytes = nullomer_count * byte_size
                    null_bytes = f.read(total_bytes)
                    if len(null_bytes) < total_bytes:
                        break
                    v2s = struct.unpack(f'<{nullomer_count}{byte_format}', null_bytes)
                trie.insert_from_bit(tuple(v1_bits), v2s)
        except (struct.error, OSError):
            pass
    return trie
