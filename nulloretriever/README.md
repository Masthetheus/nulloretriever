# Nulloretriever Python Package

This package contains the shared Python code used by the Nulloretriever workflow and optional scripts. It provides the trie representation for nullomer bitfiles, genome download helpers, analysis helpers, and common utilities such as path and integrity checks.

## Structure
```
nulloretriever/
├── README.md
├── analysis
│   ├── composition.py
│   ├── motifs.py
│   └── processing.py
├── core
│   └── triebit_class.py
├── data
│   ├── ncbiapidata.py
│   └── ncbidownload.py
└── utils
    ├── bit_encoding.py
    ├── integrity.py
    ├── paths.py
    ├── progress_bar.py
    └── validation.py
```
## Installation

The package is installed automatically with the Conda environment. For development:
```
pip install -e .
```

This installs the package in editable mode, so changes to the source code are reflected immediately.

### Core Module:

- **triebit_class**

The main class is TrieBit, which implements a compressed trie for k-mer storage.

**Examples**

Basic scripts operation consists of:
```
from nulloretriever.core.triebit_class import TrieBit
# Initialize a TrieBit object
# Variables are calculated related to k size, further
# information can be found in the main README file
trie = TrieBit(m, k, half_k)

# Insert a kmer split in v1/v2 into the trie
trie.insert(v1,v2)

# Insert a v1 with multiple v2s
trie.insert_v2_list(v1,v2s)

# Write final trie into a compact binary file
trie.write_bit_format(ouput_path)

# If full sequence output is desired, use:
# Where identifier relates to the origin of the trie
# If from genome, use 0, if mounted from bit file, use 1
trie.write_sequences(output_path, identifier)
```

The retrieval of statistiscs or operation on the saved bit file can be made with:
```
# Load a bit file into a trie object
# Trie parameters are inside the bit file
mount_trie_from_bitfile(null_bit_file)

# Get number of: nullomers if runned on top of mounted trie or
# the number of unique kmers found on the original genome
# when on top of genome mounted trie
count = trie.count_kmers()

# Compute GC content
gc_stats = trie.count_gc()

# Retrieve motif statistics
motifs_results = {
"cpg": trie.retrieve_nullomers_cpg_stats(),
"palindromy": trie.retrieve_palindrome_stats(),
"homopolymers": trie.retrieve_homopolymer_stats(),
}
return motifs_results
```

### Analysis Modules:

1. **composition.py**

* calculate_gc_index(index, half_k): returns GC count for given sequence index.
* generate_dc_dict(half_k): returns python dictionary with pairs index: gc count.

2. **motifs.py**

* calculate_cpg_index(index, k): returns count of C directly followed by G occurences.
* generate_cpg_dict(k): returns python dictionary with pairs index: CpG count.
* generate_complement_index_dict(half_k): returns python dictionary with pairs idx: idx to form palindrome.
* generate_homopolymer_array(half_k): returns array of homopolymer indexes for given k value.

3. **processing.py**

* mount_trie_from_bitfile(filename): creates a TrieBit object from one bit nullomer file.

### Data Modules:

1. **ncbiapidata.py**

* read_accession_txt(filepath): reads NCBI accessions from simple txt format.
* read_accession_list(filepath, column=None): reads NCBI accessions from a csv, tsv or txt file(invoking the function above).
* get_accession_summary_data(acc): returns xml summary data for given NCBI accession code.
* get_genome_download_link(accessions): gather the list of Entrez genome links for given assembly accessions.
* get_genome_metadata(accessions, params=None): gather metadata for one accession code.
* get_genome_length(metadata): searches organism metadata for it's length.
* get_taxonomy_metadata(metadata): searches for organism TaxID informations on NCBI db.

2. **ncbidownload.py**

* download_genome_bioentrez(accessions, out): download the accessions genomes from NCBI with Bio.Entrez.

### Utils Modules:

1. **integrity.py**

* unzip_fasta_file(gzip_path): auxiliary for genome decompressing after API download.
* casplock_file(target_file): capitalizes the entire file.
* check_genome_integrity(genome): checks if the genome is comprised only by A, T, C and G.
* check_multiple_genomes_integrity(genomes): checks if an array of genomes passes integrity test.

2. **paths.py**

* gather_files_names(location): gather all files names from given directory.
* gather_files_paths(location): gather the specific path of all files inside location.

3. **progress_bar.py**: progress_bar(current, total, start, size=30): custom progress_bar for genome download.

4. **validation.py**

* get_valid_email(): prompts for NCBI-registered email.
* get_valid_tool(): prompts for tool name.

5. **bit_encoding.py**

* generate_sub_indexes(idx, k): Calculates v1 and v2 from given idx.
* encode_to_bit(filename, k): Reads a sequence nullomer file and retrieve the nullomer set in bit format.
* write_bit_file(filename, output, k): Writes the retrieved nullomer set in bit format to output.

### Usage in Snakemake

Most functions are designed to be called from Snakemake and optional scripts, receiving input and output files during script execution. Even though not being designed for standalone use, each module provides functions that accept explicit file paths and allow such.

**All dependencies are listed in the root requirements.txt.**

## Contributing

This package is part of the Nulloretriever project. For changes, please follow the existing code style and update the relevant docstrings.

## License

GNU General Public License v3.0 (see LICENSE in the repository root).

Last updated: 2026-07-06
