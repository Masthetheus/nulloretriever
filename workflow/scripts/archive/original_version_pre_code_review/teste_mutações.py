from funcoes_mestrado import *
import os
import csv
import pandas as pd
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

base_path = '/home/leveduras/integranteslab/matheus/resultados_novos/'
config_file = '/home/leveduras/integranteslab/matheus/Mestradoteste/Scripts/snakemakefluxo_c_trie_bit/config/config.yaml'
x = int(input("Digite o número de processos a serem utilizados (default 4): ")) or int(4)
paths, erros = path_nulomeros_gerados(base_path, config_file)
print(f'Não foram gerados arquivos de nulômeros para {len(erros)} organismos.')

for k, lista_paths in paths.items():
    l = int(k / 2)
    m = 4 ** l
    if k > 10:
        x = 10
    csv_path = f"{base_path}{k}/mutacoesk_{k}.csv"
    print(f"Iniciando o processamento de {len(lista_paths)} organismos para k = {k}")
    inicio = time.time()
    orgs_com_nulomeros = [(g, a) for (g, a) in lista_paths if total_null_txt(a, l) > 0]
    tot = len(orgs_com_nulomeros)
    tasks = [(genoma, arquivo_txt, l, m, k) for (genoma, arquivo_txt) in orgs_com_nulomeros]
    results = []
    with ProcessPoolExecutor(max_workers=x) as executor:
        for i, res in enumerate(executor.map(processa_organismo, tasks), 1):
            barra_progresso_tempo(i, tot, inicio)
            results.extend(res)
    if results:
        print(f"Total de nulômeros encontrados para k = {k}: {len(results)}")
        print(f"Salvando resultados em {csv_path}...")
        with open(csv_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["kmer", "nulomeros_proximos"])
            for row in results:
                writer.writerow(row)
        print(f"Resultados salvos em {csv_path}")
    else:
        print(f"Nenhum nulômero de mutação encontrado para k = {k}.")
    print(f"Foram importados com sucesso os arquivos de nulômeros para k = {k}")
    