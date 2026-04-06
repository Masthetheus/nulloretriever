"""Receives bit nullomer files for two consecutive k values and compare."""

from nulloretriever.analysis.parsings import read_per_v1, sequence_length_conversion
from nulloretriever.analysis.counter import quick_nullomer_count
from nulloretriever.analysis.trivial_extensions import first_half_extensions, second_half_extensions

file1 = "/home/leveduras/integranteslab/matheus/nulloretriever/workflow/results/k12/GCA_000412225_2/null_bit_format"
file2 = "/home/leveduras/integranteslab/matheus/nulloretriever/workflow/results/k13/GCA_000412225_2/null_bit_format"
filetestk8 = "/home/leveduras/integranteslab/matheus/nulloretriever/workflow/results/k8/GCA_050947625_1/null_bit_format"
skip_counter = 0
nullomer_count = 0
count = quick_nullomer_count(file1)
print(count)
print(skip_counter)
half_k = 6
extensions = {}
while nullomer_count < count:
    v2_extensions = set()
    v1, v2s, counter, nullomer_count, skip_counter = read_per_v1(
        file1, skip_counter, nullomer_count)
    extensions.update(first_half_extensions(v1, v2s, half_k, file2))
    i = half_k - 1
    found = second_half_extensions(v1, v2s, half_k, file2)
    v1_extensions = set()
    checker = 10
    indexv1 = 0
print(len(v1_extensions))
print(len(found))
