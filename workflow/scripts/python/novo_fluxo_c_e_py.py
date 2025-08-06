from funcoes_mestrado import *
import re,os,sys,glob
import subprocess


genoma = snakemake.input.genome
output = snakemake.output.nullomers
org = genoma.split("/")[-1]
k = snakemake.params.k
l = int(k / 2)
m = 4**l
base_values = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
saida_kmers = f"/tmp/{org}_kmers.txt"

# Inicializa a trie
trie = inicializar_triebit_teste(l, m)

# Executa o programa C e lê os k-mers do stdout
proc = subprocess.Popen(
    ['/home/leveduras/integranteslab/matheus/Mestradoteste/./fasta_kmers_novo', genoma, str(k)],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True
)

kmers_inseridos = 0
print("DEBUG: Lendo k-mers do programa C...")
for line in proc.stdout:
    line = line.strip()
    if line:
        # Dividir por espaço para obter k-mers individuais
        kmers_na_linha = line.split(' ')
        for kmer_str in kmers_na_linha:
            if kmer_str.strip():
                try:
                    kmer_array = [int(x) for x in kmer_str.split(',')]
                    if len(kmer_array) == k:
                        v1 = tuple(kmer_array[:l])
                        v2_array = kmer_array[l:]
                        
                        # Converte v2_array em um índice único
                        v2 = 0
                        for base in v2_array:
                            v2 = (v2 << 2) | base
                        
                        trie.insert(v1, v2)
                        kmers_inseridos += 1
                        if kmers_inseridos <= 5:
                            print(f"DEBUG: K-mer {kmers_inseridos} inserido - v1={v1}, v2={v2}")
                except Exception as e:
                    print(f"DEBUG: Erro ao processar k-mer '{kmer_str}': {e}")
                    continue

proc.wait()
print(f"Total de k-mers inseridos na trie: {kmers_inseridos}")
print(f"Valores calculados: k={k}, l={l}, m={m}")

print("Iniciando a escrita da trie em arquivo...")
escrever_trie_em_txt_bitarray(trie, output, l)
tot = contar_nulomeros_trie_bit_novo(trie, l)
print("Trie escrita com sucesso!")
print("Contagem de nulômeros concluída com sucesso!")
print(f"Total de nulômeros encontrados: {tot}")
