"""Functions aimed at comparing bit tries searching for prime nullomers."""

import struct
from collections import defaultdict


def search_prime_nullomers(common_null_set, common_v1, file2_v1_loc, file2):
    """Search for prime nullomers between two files."""
    new_common_null_set = defaultdict(set)
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    with open(file2, 'rb') as f:
        f.seek(8)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        for v1 in common_v1:
            counter = 0
            current_v2_values = set()
            f.seek(file2_v1_loc[v1], 0)
            try:
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(
                    f'<{byte_format}', nullomer_count_bytes)[0]
                while counter < nullomer_count:
                    v2_bytes = f.read(byte_size)
                    current_v2_values.add(struct.unpack(
                        f'<{byte_format}', v2_bytes)[0])
                    counter += 1
                new_common_null_set[v1] = current_v2_values & common_null_set[v1]
            except (struct.error, OSError) as e:
                print(e)
                pass
    return new_common_null_set
