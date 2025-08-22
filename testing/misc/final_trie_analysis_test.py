from funcoes_mestrado import *
import os
import csv
import pandas as pd

# config_file = "../../config/config.yaml"
# k = input("Which k value was used?")
# if not k:
#     k = 10
# org = input("Which organism to perform the final analysis?")
# if not org:
#     null_file = "../data/teste_test_run"
#     org = "teste"
# else:
#     null_file = "../data/" + org + "_test_run"
# csv_out = "../data/" + org + "_data.csv"
# final_analysis_testing(csv_out,null_file,k,org,config_file)
l = 4
org = input("Which genome to test? (default = teste)")
if not org:
    genome = "../data/teste"
    org = "teste"
else:
    genome = "../data/" + org
out = "../data/"+org+"_test_run"
tot = quick_nullomer_count(out)
print(tot)
gc = nullomers_gc_mean(out)
print(gc)
cpg_stats = nullomers_cpg_presence_mean(out)
print(cpg_stats)