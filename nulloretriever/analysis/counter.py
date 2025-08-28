"""Defines the function used for fast nullomer counting from a TrieBit bit output file"""

import struct

def quick_nullomer_count(filename):
    """Fast counting of nullomeric sequences on TrieBit bit output
    Uses the v2 occurence counter already available for each v1 value.
    Args:
        filename(str): path to the TrieBit bit file location
    Returns:
        count(int): total of nullomeric sequences found in given genome for a specific k value
    """
    count = 0
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l = struct.unpack('<H', f.read(2))[0]
        l = int(l)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                # Skip v2 indices (nullomer IDs)
                f.seek(nullomer_count * byte_size, 1)
                
                # Add the number of nullomers for this v1
                count += nullomer_count
                
        except (struct.error, OSError):
            pass
    
    return count

