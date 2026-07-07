# C Source Code - kmer_extractor

This directory contains the C implementation of the k-mer extractor used by the Nulloretriever pipeline. The code is compiled into a single binary (kmer_extractor) that reads a FASTA file, extracts all k-mers (forward and reverse complement), and outputs their packed indices in binary format.

Structure

scripts/c/
├── src/                            # base scripts for the kmer_extractor
│   ├── main.c
│   ├── kmer_encoding.c             # sequence encoding functions
│   ├── kmer_extraction.c           # general logic of sequencing processing
│   └── kmerio.c                    # concentrates I/O operations
├── include/                        # necessary headers
│   ├── fasta_kmer_extraction.h
│   ├── kmer_encoding.h
│   ├── kmer_extraction.h
│   └── kmerio.h
├── build/                    # (object files, generated during compilation)
├── bin/                      # (compiled binary: kmer_extractor)
└── Makefile

## File Descriptions (src/)

**main.c**
Entry point and command-line interface. Parses arguments (<fasta_file> <k>), reads the FASTA file line by line, and calls process_kmers() for each sequence. Currently processes only both strands and simultaneously, saving reading and parsing time.
Handles memory allocation and cleanup.

**kmer_encoding.c**
Implements 2-bit encoding of DNA bases (A, C, G, T) into packed integers.
- encode_kmer(seq, k): converts a DNA string to a uint64_t using bitwise operations, reading one base at a time as it's ASCII correspondent.
- encode_rev(val, k): returns the reverse complement of a packed k-mer.
- decode_kmer(val, k): converts a packed integer back to a DNA string (for debugging).

**kmer_extraction.c**
Core extraction logic in process_kmers():
- For each position in the sequence, encodes the k-mer and its reverse complement.
- Checks if each has been seen using a bitarray.
- If not seen, appends the packed value to an output buffer.
- Periodically flushes the buffer to stdout.

**kmerio.c**
Manages binary output:
- flush_output_buffer(): writes the accumulated buffer to stdout.
- write_decoded_kmer(): prints a human-readable k-mer (debugging only).

## Headers (include/)

Each .c file has a corresponding .h file in include/ with function prototypes and macros.
- fasta_kmer_extraction.h: aggregates all headers and declares main functions.
- kmer_encoding.h: encoding functions.
- kmer_extraction.h: extraction functions and helpers.
- kmerio.h: output functions.

## Compilation

The Makefile compiles all source files in src/. Flags used:
- -fPIE and -pie for Position Independent Executable support.
- -O3 for optimization.
- -Wall -Wextra for warnings.

Run make in this directory to build the binary. The output by default is placed in bin/kmer_extractor.

## Usage

From the command line:
```
./kmer_extractor genome.fasta 15 > output.bit
```
The binary writes binary data to stdout, so you must redirect the output to a file. Keep in mind the intended use circumvent this by reading and processing directly the stdout via auxiliary functions inside a Python script.

## Data Flow

1. FASTA file is read by main.c.
2. For each sequence, process_kmers() extracts all k-mers and their reverse complements.
3. Each k-mer is encoded as a uint64_t.
4. The value is checked against a bitarray to avoid duplicates.
5. If not seen, the packed index is appended to an output buffer.
6. The buffer is flushed periodically to stdout.
7. The resulting .bit file contains the packed indices for downstream Python processing.

## Dependencies

- GCC (or any C99-compliant compiler)
- Standard C library (no external dependencies)

## Troubleshooting

- "Usage: ..." error: provide exactly two arguments: FASTA file and k.
- Compilation errors: ensure all files are present and include paths are correct.
- Permission denied on binary: run chmod +x bin/kmer_extractor after compilation.

## Maintenance Notes

- The bitarray allocated is calloc((total + 7) / 8, 1) where total = 1ULL << (2*k). This can be memory-intensive for large k (>20). The code currently limits k to 20 via MAX_K. Future work for better duplicated resolving is on the development path.
- The MAX_SEQ buffer by default is set to 8 GB – adjust if needed for very large genomes.

## License

This code is part of the Nulloretriever project and is licensed under GNU GPL v3.0.

Last updated: 2026-07-07
