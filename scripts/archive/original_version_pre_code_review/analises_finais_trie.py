from funcoes_mestrado import *
import os
import csv
import pandas as pd

csv_total_null = '/home/leveduras/integranteslab/matheus/resultados_novos/total_nulomeros.csv'
csv_db_orgs = "/home/leveduras/integranteslab/matheus/db/teste.tsv"
base_path = '/home/leveduras/integranteslab/matheus/resultados_novos/'
config_file = '/home/leveduras/integranteslab/matheus/Mestradoteste/Scripts/snakemakefluxo_c_trie_bit/config/config.yaml'

processar_null_batch(csv_total_null, base_path, config_file)