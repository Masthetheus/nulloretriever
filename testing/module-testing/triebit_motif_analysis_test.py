from nulloretriever.analysis.motifs import *

file = "../data/teste_test_run"
print(f"Starting motifs analysis test on the file {file}")
print("=== Gathering CpG stats ===")
cpg_stats = retrieve_nullomers_cpg_stats(file)
print("=== Stats gathered! ===")
print(f"Total cpg occurrences = {cpg_stats[0]}\n")
print(f"Nullomers with at least 1 cpg = {cpg_stats[1]}\n")
print(f"CpG organism mean = {cpg_stats[2]}\n")
print(f"Number of CpG's inside nullomers mean = {cpg_stats[3]}")
print("=== Gathering Palindrome stats ===")
palindrome_count, palindrome_relative = retrieve_palindrome_stats(file)
print("=== Stats gathered! ===")
print(f"{palindrome_count} palindromes found on the organism. This represents {palindrome_relative:.2f}% of the nullomers of the current organism!")
print("=== Gathering Homopolymer stats ===")
found_homopolymers = retrieve_homopolymer_stats(file)
print("=== Stats gathered! ===")
print("The following homopolymers were found for the current organism:")
for i in found_homopolymers:
    print(i)