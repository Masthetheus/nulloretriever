"""Functions aimed at composition analysis of nullomeric sequences from bit files"""

import struct

def calculate_gc_index(index, l):
    """Calculate nullomeric gc content for each relative index
    Args:
        index(int): index related to a k-mer sequence
        l(int): size of the k-mer original sequence
    Returns:
        gc(int): number of occurences of G or C nucleotides on the sequence representend by the given index
    """
    gc = 0
    for i in range(l):
        base = (index // (4 ** (l - i - 1))) % 4
        if base == 2 or base == 3:
            gc += 1
    return gc

def generate_gc_dict(l):
    """Generate a dict pairing each possible index to it's total GC count
    Args:
        l(int): size of the k-mer original sequence
    Returns:
        gc_dict(dict): dictionary with pairs index:gc_count for all possible k-mers sequence of lenght l
    """
    gc_dict = {}
    for index in range(4**l):
        gc = calculate_gc_index(index, l)
        gc_dict[index] = gc
    return gc_dict

def nullomers_gc_mean(filename):
    """Calculate mean GC% of all nullomeric sequences on given organism
    Args:
        filename(str): TrieBit bit file
    Returns:
        gc_percent(float): mean of GC presence in all nullomeric sequences of given organism for given k value
    """
    count = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        l = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        gc_dict = generate_gc_dict(l)
        gc_tot = 0
        v1_count = 0
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                v1_count += 1
                gc_tot += gc_dict[v1]
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                i = 0
                while i < nullomer_count:
                    nullomer_byte = f.read(byte_size)
                    nullomer_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                    gc_tot += gc_dict[nullomer_index]
                    i += 1
                count += nullomer_count    
        except (struct.error, OSError):
            pass
    total_bases = (v1_count*l)+(count*l)
    gc_percent = (gc_tot/total_bases)*100
    return gc_percent