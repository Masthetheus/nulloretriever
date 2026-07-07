# Python Scripts – Nulloretriever Workflow

This directory contains the Python scripts called by the Snakemake rules in the Nulloretriever pipeline. They handle genome integrity checking, k-mer extraction, statistics retrieval, and summary aggregation. All scripts are designed to run within the Snakemake environment and receive inputs/outputs via the snakemake object.

## Scripts Overview

1. nullomer_extraction.py
2. nullomer_statistics_retrieval.py
3. nullomer_statistics_summary.py
4. genome_checking.py
5. path_checking.py


1. nullomer_extraction.py

Purpose:
- Executes the C binary (kmer_extractor) on a given genome and k value.
- Reads the binary output stream, decodes each packed k-mer, and inserts it into a TrieBit structure.
- Writes the resulting TrieBit to a .bit file for downstream processing.

Inputs (from Snakemake):
- snakemake.input.bin: path to the C binary (kmer_extractor).
- snakemake.input.genome: path to the genome FASTA file.
- snakemake.params.k_val: k-mer size (integer).
- snakemake.output[0]: path to the output .bit file.

Algorithm:
- Spawns the C binary via subprocess, redirecting stdout to a pipe.
- Reads bytes_per_sequence = ceil((2*k)/8) bytes per k-mer.
- Decodes each packed index, splits into v1 and v2 based on k parity. For odd values, |v1|>|v2|.
- Inserts into the TrieBit using trie.insert(v1_bits_tuple, v2).
- Writes the trie to the output .bit file.

Dependencies:
- subprocess (standard library)
- nulloretriever.core.triebit_class.TrieBit
- snakemake.script.snakemake (for access to Snakemake objects)

Integration with Snakemake:
- Called from the rule extract_nullomers.
- Uses wildcards {k} and {organism} to generate dynamic paths.

Notes:
- The script expects the C binary to output packed indices in big-endian byte order.
- The number of bytes per k-mer is computed as (2*k + 7) // 8.
- The TrieBit is constructed with parameters m (number of nodes), k, and half_k.
- The reverse complement is already handled by the C binary; this script only needs to process the forward index stream.

2. nullomer_statistics_retrieval.py

Purpose:
- Loads a TrieBit from a .bit file for a given k and organism.
- Computes specified statistics: composition (GC content), counter (number of k-mers), motifs (CPG, palindromy, homopolymers), and optionally trivial extensions.
- Writes the statistics to a CSV file.

Inputs (from Snakemake):
- snakemake.input.current: path to the .bit file for the current k.
- snakemake.input.previous (optional): path to the .bit file for k-1 (used for trivial extensions).
- snakemake.params.stats: list of statistics to compute (e.g., ["composition", "counter", "motifs"]).
- snakemake.wildcards.organism: organism identifier.
- snakemake.wildcards.k: k-mer size (as a string).
- snakemake.output[0]: path to the output CSV file.

Statistics computed:
- composition: GC content of nullomers (via trie.count_gc()).
- counter: number of nullomers (via trie.count_kmers()).
- motifs: CPG, palindromy, homopolymer statistics (via motif_wrapper()).
- trivial (only if previous .bit exists): trivial extensions (via trie.find_trivial_ext()).

Output CSV columns:
- organism, k, counter, composition_gc_mean, composition_gc_std, etc.
- For motifs, columns like motifs_cpg_*, motifs_palindromy_*, motifs_homopolymers_*.
- If trivial is computed, columns like trivial_* are added.

Dependencies:
- csv (standard library)
- nulloretriever.analysis.processing.mount_trie_from_bitfile
- snakemake.script.snakemake

Integration with Snakemake:
- Called from the rule retrieve_nullomer_statistics.
- Uses wildcards {k} and {organism}.
- The previous input is conditionally provided by the get_previous_k_input() function in the Snakefile.

Notes:
- The script flattens nested dictionaries (e.g., motifs) using dict_flattener().
- If a statistic is not available (e.g., trivial without previous file), it is **skipped with a warning.**
- The organism and k values are extracted from Snakemake wildcards and added to the output.


3. nullomer_statistics_summary.py

Purpose:
- Aggregates all individual nullomer statistics CSVs (one per organism per k) into a single summary CSV.

Inputs (from Snakemake):
- snakemake.input: list of all nullomer_statistics.csv files (one for each combination of k and organism).
- snakemake.output[0]: path to the aggregated CSV file.

Algorithm:
- Reads each input CSV with pandas.
- Concatenates all dataframes row-wise.
- Writes the combined dataframe to the output CSV.

Dependencies:
- pandas
- snakemake.script.snakemake

Integration with Snakemake:
- Called from the rule summarize_nullomer_statistics.
- The input list is generated using expand() in the Snakefile.

Notes:
- The script ignores missing files (prints a warning and continues) – this is useful if some runs failed.
- The output CSV has the same columns as the individual CSVs, with rows from all organisms and k values.

**All dependencies are listed in the root environment.yaml and pyproject.toml.**

## Integration with Snakemake Workflow

- nullomer_extraction.py is called by the extract_nullomers rule.
- nullomer_statistics_retrieval.py is called by the retrieve_nullomer_statistics rule.
- nullomer_statistics_summary.py is called by the summarize_nullomer_statistics rule.

The scripts rely on the snakemake object, which is automatically injected by Snakemake when using the script directive. This object provides input, output, params, and wildcards attributes.

Maintenance Notes

- If you add a new statistic, update the STATISTICS list in the config.yaml and ensure the dispatch_table in nullomer_statistics_retrieval.py includes it.
- The trivial extension logic requires a previous .bit file (k-1).
- The nullomer_extraction.py script assumes the C binary outputs exactly bytes_per_sequence bytes per k-mer. Any deviation will cause decoding errors.
- For large genomes and high k, the TrieBit structure may use significant memory, being capped by the number of unique sequences present, that being the number of v1's that do happen on the genome.

License

This code is part of the Nulloretriever project and is licensed under GNU GPL v3.0.

Last updated: 2026-07-07
