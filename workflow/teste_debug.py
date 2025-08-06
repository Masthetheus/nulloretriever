#!/usr/bin/env python3
import sys
import os
sys.path.append('/home/leveduras/integranteslab/matheus/Mestradoteste/Scripts/nullomer_c_trie_bit/workflow/scripts/python')

from funcoes_mestrado import *
import re,os,sys,glob
import subprocess

genoma = "/home/leveduras/integranteslab/matheus/Mestradoteste/Scripts/nullomer_c_trie_bit/workflow/data/genomes/GCA_000009045_1"
output = "/tmp/teste_nullomers.txt"
org = genoma.split("/")[-1]
k = 16
l = int(k / 2)
m = 4**l
base_values = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
saida_kmers = f"/tmp/{org}_kmers.txt"

print(f"DEBUG: Parâmetros - genoma={genoma}, k={k}, l={l}, m={m}")

# Inicializa a trie
print("DEBUG: Inicializando trie...")
trie = inicializar_triebit_teste(l, m)
print(f"DEBUG: Trie inicializada com {m} posições")

# Executa o programa C e lê os k-mers do stdout
print("DEBUG: Executando programa C...")
proc = subprocess.Popen(
    ['/home/leveduras/integranteslab/matheus/Mestradoteste/./fasta_kmers_novo', genoma, str(k)],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

kmers_inseridos = 0
linha_count = 0

print("DEBUG: Lendo k-mers do programa C...")
for line in proc.stdout:
    linha_count += 1
    if linha_count <= 5:
        print(f"DEBUG: Linha {linha_count}: {line.strip()[:100]}")
    
    line = line.strip()
    if line:
        # Dividir por espaço para obter k-mers individuais
        kmers_na_linha = line.split(' ')
        for kmer_str in kmers_na_linha:
            if kmer_str.strip():
                try:
                    kmer_array = [int(x) for x in kmer_str.split(',')]
                    if len(kmer_array) == k:
                        v1 = kmer_array[:l]
                        v2 = kmer_array[l:]
                        
                        inserir_na_trie_bit_teste(trie, v1, v2, l)
                        kmers_inseridos += 1
                        if kmers_inseridos <= 5:
                            print(f"DEBUG: K-mer {kmers_inseridos} inserido - v1={v1}, v2={v2}")
                except Exception as e:
                    print(f"DEBUG: Erro ao processar k-mer '{kmer_str}': {e}")
                    continue

proc.wait()
print(f"DEBUG: Programa C terminou. Total linhas lidas: {linha_count}")
print(f"DEBUG: Total k-mers inseridos: {kmers_inseridos}")

# Escreve nullomers
print("DEBUG: Escrevendo nullomers...")
escrever_trie_em_txt_bitarray(trie, output, l, m)

print(f"DEBUG: Arquivo gerado: {output}")
print("DEBUG: Primeiras linhas do arquivo:")
with open(output, 'r') as f:
    for i, line in enumerate(f):
        if i < 10:
            print(f"  {line.strip()}")
        else:
            break

# Contar total de nullomers
with open(output, 'r') as f:
    total_nullomers = sum(1 for line in f)
print(f"DEBUG: Total de nullomers encontrados: {total_nullomers}")
print(f"DEBUG: Esperado máximo de nullomers: {m}")
