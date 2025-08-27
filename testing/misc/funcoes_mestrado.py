from itertools import product
import random
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import datetime as dt
from bitarray import bitarray
import csv
import yaml
import re
import pandas as pd
import os
from collections import Counter, defaultdict
import time
import subprocess
import psutil
from memory_profiler import memory_usage
import struct

#----------------------------------------------------------------------#
# Funções em teste
#----------------------------------------------------------------------#
def bit_file_stat_retrieve_template(filename):
    """
    Fastest way to get the total count of nullomer k-mers.
    Sums all nullomer counts across all v1 blocks in the file.
    """
    count = 0
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l = struct.unpack('<H', f.read(2))[0]
        l = int(l)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                # Skip v2 indices (nullomer IDs)
                f.seek(nullomer_count * byte_size, 1)
                
                # Add the number of nullomers for this v1
                count += nullomer_count
                
        except (struct.error, OSError):
            pass
    
    return count


#----------------------------------------------------------------------#
# Funções ocasionais
#----------------------------------------------------------------------#

def genome_random(lengthsG, genomaaleatorio): #apenas para testes
    genoma=""
    for count in range(lengthsG):
        genoma+=random.choice("ATCG")
    with open(genomaaleatorio, "w") as datafile:
        for i in genoma:
            datafile.write(str(i))
        datafile.close()

def gerar_kmers(l):
    """
    Gera todos os k-mers possíveis de comprimento `l` usando as bases A, T, C, G.
    :param l: Comprimento do k-mer.
    :return: Lista de k-mers.
    """
    from itertools import product
    bases = ['A', 'T', 'C', 'G']
    return [''.join(kmer) for kmer in product(bases, repeat=l)]

def ajustar_nome_orgs(nome_org):
    nome_org = nome_org.replace(" ", "_")  # Substitui espaços por underscores
    return re.sub(r'[^a-zA-Z0-9_]', '_', nome_org)  # Remove caracteres especiais, mantendo apenas letras, números e underscores

def adicionar_genome_size(csv_total_null, csv_db_orgs):
    df_null_cont = pd.read_csv(csv_total_null, sep=",")
    df_orgs = pd.read_csv(csv_db_orgs, sep="\t")

    df_orgs['Organism Name Normalizado'] = df_orgs['Organism Name'].apply(normalizar_nome)
    df_null_cont['Organismo Normalizado'] = df_null_cont['Organismo'].apply(normalizar_nome)

    df_null_cont = pd.merge(
        df_null_cont,
        df_orgs[['Organism Name Normalizado', 'Genome Size']],
        left_on='Organismo Normalizado',
        right_on='Organism Name Normalizado',
        how='left'
    )

    df_null_cont = df_null_cont.drop(columns=['Organism Name Normalizado', 'Organismo Normalizado'])
    df_null_cont.to_csv(csv_total_null, sep=",", index=False)
    if 'Genome Size' in df_null_cont.columns:
        df_null_cont = df_null_cont.drop(columns=['Genome Size'])

def porcent_null_genoma(arquivo):
    """
    Calcula a porcentagem de nulômeros em relação ao tamanho do genoma.
    :param arquivo: Caminho para o arquivo CSV contendo os dados de nulômeros e tamanho do genoma.
    """
    df = pd.read_csv(arquivo, sep=",")
    df['Porcentagem Nulômeros'] = (df['Total Nulômeros'] / df['Genome Size']) * 100
    df['Porcentagem Nulômeros'] = df['Porcentagem Nulômeros'].round(2)  # Arredonda para duas casas decimais
    df.to_csv(arquivo, sep=",", index=False)

def adicionar_grupo_por_organismo(csv_principal, csv_grupos, saida_csv):
    """
    Adiciona o identificador de grupo ao csv_principal, de acordo com o nome do organismo.
    :param csv_principal: Caminho do CSV principal.
    :param csv_grupos: Caminho do CSV com colunas 'Organismo' e 'Grupo'.
    :param saida_csv: Caminho do novo CSV com a coluna de grupo adicionada.
    """
    # Carrega os dois CSVs
    df_main = pd.read_csv(csv_principal)
    df_grupos = pd.read_csv(csv_grupos, sep="\t")
    print(list(df_main.columns))
    print(list(df_grupos.columns))
    df_main['Organismo_norm'] = df_main['Organismo'].apply(ajustar_nome_orgs)
    df_grupos['Organismo_norm'] = df_grupos['Organism Name'].apply(ajustar_nome_orgs)

    # Faz o merge pelo nome normalizado
    df_merged = pd.merge(df_main, df_grupos[['Organismo_norm', 'Group_ID']], on='Organismo_norm', how='left')

    # Remove coluna auxiliar
    df_merged = df_merged.drop(columns=['Organismo_norm'])

    # Salva o resultado
    df_merged.to_csv(saida_csv, index=False)

def obter_ks_analisados(base_path):
    """
    Obtém os valores de k analisados a partir dos arquivos no diretório base_path.
    :param base_path: Caminho base onde os arquivos estão localizados.
    :return: Lista de valores de k analisados.
    """
    ks = []
    for nome in os.listdir(base_path):
        caminho = os.path.join(base_path, nome)
        print(f"Verificando caminho: {caminho}")
        if os.path.isdir(caminho) and nome.isdigit():
            ks.append(int(nome))
    return sorted(ks)

def path_nulomeros_gerados(base_path, config_file):
    orgs_bruto = obter_lista_org(config_file)
    orgs = [ajustar_nome_orgs(org) for org in orgs_bruto]
    erros = []
    k_values = obter_ks_analisados(base_path)
    paths_existentes = {}
    for k in k_values:
        l = int(k // 2)
        k_path = base_path + str(k) + '/'
        print(f"Verificando diretório para k={k}")
        cont = 0
        paths_k = []
        paths_genoma = []
        for org in orgs:
            org_path = k_path + org + '/'
            genoma_path = '/home/leveduras/integranteslab/matheus/Mestradoteste/Scripts/snakemakefluxo_c_trie_bit/dados/genomas/genomas/' + org
            arquivo_txt = org_path + f'nulomerostrie_{org}_{k}.txt'
            if not os.path.exists(arquivo_txt):
                erros.append([org, k, 'Arquivo não encontrado'])
                continue
            if not os.path.exists(genoma_path):
                erros.append([org, k, 'Genoma não encontrado'])
                continue
            # Adiciona a tupla (caminho_genoma, caminho_k) para este organismo
            paths_k.append((genoma_path, arquivo_txt))
        if paths_k:
            paths_existentes[k] = paths_k
    return paths_existentes, erros

def barra_progresso_tempo(atual, total, inicio, tamanho=30):
    agora = time.time()
    elapsed = agora - inicio
    if atual > 0:
        tempo_medio = elapsed / atual
        restante = tempo_medio * (total - atual)
    else:
        restante = 0
    minutos, segundos = divmod(int(restante), 60)
    proporcao = atual / total
    completos = int(proporcao * tamanho)
    barra = "|" + "o" * completos + "-" * (tamanho - completos) + "|"
    print(f"\r{barra} {atual}/{total} - Estimativa: {minutos:02d}:{segundos:02d} restantes", end='', flush=True)
    if atual == total:
        print()  # Garante quebra de linha ao final

#----------------------------------------------------------------------#
# TrieBit approach
#----------------------------------------------------------------------#

# Base definitions

class TrieNodeBitTeste:
    def __init__(self, m):
        self.children = [None] * 4  # 0:A, 1:T, 2:C, 3:G
        self.v2_set = bitarray(m)
        self.v2_set.setall(0)

    def iterate(self, path=None):
        if path is None:
            path = []
        yield path, self
        for i, child in enumerate(self.children):
            if child is not None:
                yield from child.iterate(path + [i])

class TrieBitTeste:
    def __init__(self, m):
        self.root = TrieNodeBitTeste(m)

    def insert(self, v1, v2):
        node = self.root
        for value in v1:
            if node.children[value] is None:
                node.children[value] = TrieNodeBitTeste(len(node.v2_set))
            node = node.children[value]
        node.v2_set[v2] = 1  # Marca o bit correspondente

    def iterate(self):
        yield from self.root.iterate([])

def inicializar_triebit_teste(l, m):
    trie = TrieBitTeste(m)
    def construir(node, depth):
        if depth == l:
            return
        for i in range(4):
            if node.children[i] is None:
                node.children[i] = TrieNodeBitTeste(m)
            construir(node.children[i], depth + 1)
    construir(trie.root, 0)
    return trie

def revcomp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

# Misc

def count_trie_nullomers(trie_bit, l):
    def dfs(node, depth):
        total = 0
        if depth == l:
            # Conta os bits 0 em cada conjunto
            m = len(node.v2_set)
            for v2 in range(m):
                if not (node.v2_set[v2]):
                    total += 1
            return total
        for child in node.children:
            if child is not None:
                total += dfs(child, depth + 1)
        return total
    return dfs(trie_bit.root, 0)

# Output functions (bit or compact txt format)

def write_trie_compact_txt(trie, arquivo_saida, l):
    """
    Escreve os caminhos e valores de uma Trie baseada em bitarray (apenas v2_set) em um arquivo de texto no formato:
    >indice_lexicografico_v1
    v2_1,v2_2,v2_3,...
    """
    with open(arquivo_saida, 'w') as f:
        def dfs(node, path):
            # Verifica se há nullômeros (bits 0) neste nó
            if len(path) == l:  # Só processa nós folha (profundidade l)
                nullomeros = [i for i, bit in enumerate(node.v2_set) if not bit]
                if nullomeros:  # Se há nullômeros para escrever
                    # Calcula o índice lexicográfico de v1 (o caminho atual)
                    indice_lexicografico = sum(base * (4 ** (l - i - 1)) for i, base in enumerate(path))
                    f.write(f">{indice_lexicografico}\n")
                    valores = ",".join(str(i) for i in nullomeros)
                    f.write(f"{valores}\n")
            for child_value, child_node in enumerate(node.children):
                if child_node is not None:
                    dfs(child_node, path + [child_value])
        dfs(trie.root, [])

def write_trie_bit_format(trie, output, l):
    """
    Saves TrieBitTeste to a compact binary format.
    Format: [header][nodes...]
    Header: b'TRIE'[4] + version(2) + l(2) + format_code(1)
    """
    m = 4**l
    if m <= 256:
        index_format = 'B'
        format_code = 1
    elif m <= 65536:
        index_format = 'H' 
        format_code = 2
    elif m <= 4294967296:
        index_format = 'I'
        format_code = 4
    else:
        index_format = 'Q'
        format_code = 8

    with open(output, 'wb') as f:
        f.write(b'TRIE')  # Magic number
        version = 1
        f.write(struct.pack('<HHB', version, l, format_code))  # version, l, byte_size
        def collect_nodes(node, path):
            if len(path) == l:
                nullomers = [i for i, bit in enumerate(node.v2_set) if not bit]
                if nullomers:
                    # Compute lexicographic index for v1 path
                    index = sum(base * (4 ** (l - i - 1)) for i, base in enumerate(path))
                    f.write(struct.pack(f'<{index_format}', index))
                    f.write(struct.pack(f'<{index_format}', len(nullomers)))
                    for v2_index in nullomers:
                        f.write(struct.pack(f'<{index_format}', v2_index))
            for child_value, child_node in enumerate(node.children):
                if child_node is not None:
                    collect_nodes(child_node, path + [child_value])
        
        collect_nodes(trie.root, [])

## Nullomer output analysis (bit output)

# Counter
def quick_nullomer_count(filename):
    """
    Fastest way to get the total count of nullomer k-mers.
    Sums all nullomer counts across all v1 blocks in the file.
    """
    count = 0
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l = struct.unpack('<H', f.read(2))[0]
        l = int(l)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                # Skip v2 indices (nullomer IDs)
                f.seek(nullomer_count * byte_size, 1)
                
                # Add the number of nullomers for this v1
                count += nullomer_count
                
        except (struct.error, OSError):
            pass
    
    return count

# GC counter
def calculate_gc_index(index, l):
    gc = 0
    for i in range(l):
        base = (index // (4 ** (l - i - 1))) % 4
        if base == 2 or base == 3:
            gc += 1
    return gc

def generate_gc_dict(l):
    gc_dict = {}
    for index in range(4**l):
        gc = calculate_gc_index(index, l)
        gc_dict[index] = gc
    return gc_dict

def nullomers_gc_mean(filename):
    """
    Counts GC% of given set of nullomers
    """
    count = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        l = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        gc_dict = generate_gc_dict(l)
        gc_tot = 0
        v1_count = 0
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                v1_count += 1
                gc_tot += gc_dict[v1]
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                i = 0
                while i < nullomer_count:
                    nullomer_byte = f.read(byte_size)
                    nullomer_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                    gc_tot += gc_dict[nullomer_index]
                    i += 1
                count += nullomer_count    
        except (struct.error, OSError):
            pass
    total_bases = (v1_count*l)+(count*l)
    gc_percent = (gc_tot/total_bases)*100
    return gc_percent

# CpG dinucleotides counter and stats calculation
def calculate_gpc_index(index, l):
    cpg = 0
    c = 0
    g = 0
    for i in range(l):
        base = (index // (4 ** (l - i - 1))) % 4
        if i == 0 and base == 3:
            g = 1
        if base == 2:
            c = 1
        elif base == 3 and c:
            cpg += 1
            c = 0
        else:
            c = 0
    return cpg,c,g

def generate_cpg_dict(l):
    cpg_dict = {}
    for index in range(4**l):
        cpg, c, g = calculate_gpc_index(index, l)
        cpg_dict[index] = (cpg,c,g)
    return cpg_dict

def retrieve_nullomers_cpg_stats(filename):
    """
    Counts the percentage of nullomers that have at least 1 CpG dinucleotide.
    """
    count = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        l = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        cpg_dict = generate_cpg_dict(l)
        cpg_tot = 0
        null_with_cpg = 0
        try:
            while True:
                end_c = 0
                cpg_exists = 0
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                if cpg_dict[v1][1] == 1:
                    end_c = 1
                if cpg_dict[v1][0] != 0:
                    cpg_exists = 1
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                i = 0
                while i < nullomer_count:
                    end_g = 0
                    nullomer_byte = f.read(byte_size)
                    nullomer_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                    cpg_tot += cpg_dict[v1][0] + cpg_dict[nullomer_index][0]
                    if cpg_dict[nullomer_index][2] == 1:
                        end_g = 1
                    if end_c and end_g:
                            cpg_tot +=1
                    if cpg_exists:
                        null_with_cpg += 1
                    elif end_c and end_g:
                        null_with_cpg += 1
                    else:
                        if cpg_dict[nullomer_index][0] > 0:
                            null_with_cpg += 1
                    i += 1
                count += nullomer_count     
        except (struct.error, OSError):
            pass
    cpg_count_mean = cpg_tot/null_with_cpg
    cpg_global_mean = (null_with_cpg/count)*100
    cpg_stats = [cpg_tot, null_with_cpg, cpg_global_mean, cpg_count_mean]
    return cpg_stats

# Palindrome stats generation
def generate_complement_index_dict(l):
    complement = {0: 1, 1: 0, 2: 3, 3: 2}
    complement_index_dict = {}
    m = 4**l
    for number in range(m):
        temp = number
        bases = []
        for _ in range(l):
            bases.append(temp % 4)
            temp //= 4
        comp_bases = [complement[base] for base in bases]
        comp_index = 0
        for i, base in enumerate(comp_bases):
            comp_index += base * (4 ** (l - i - 1))
        complement_index_dict[number] = comp_index
    return complement_index_dict

def retrieve_palindrome_stats(filename):
    """

    """
    palindrome_count = 0
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l_bytes = f.read(2)
        l = struct.unpack('<H', l_bytes)[0]
        byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        total_null = 0 
        complement_index_dict = generate_complement_index_dict(l)
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                v1_comp = complement_index_dict[v1]
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                total_null += nullomer_count
                # Skip v2 indices (nullomer IDs)
                total_bytes = nullomer_count * byte_size
                i = 0
                while i < nullomer_count:
                    nullomer_byte = f.read(byte_size)
                    total_bytes -= byte_size
                    v2_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                    if v2_index == v1_comp:
                        palindrome_count += 1
                        f.seek(total_bytes,1)
                        break
                    i += 1                
        except (struct.error, OSError):
            pass
        palindrome_relative = (palindrome_count/total_null) * 100
    return palindrome_count, palindrome_relative

# Homopolymers search
def is_homopolymer(index,l):
    same = 0
    current = 0
    for i in range(l):
        base = (index // (4 ** (l - i - 1))) % 4
        if i == 0:
            current = base
        if base == current:
            same += 1
        else:
            same = 0
            break
    return same

def generate_homopolymer_array(l):
    homopolymer_array = []
    for i in range(4**l):
        homopolymer = is_homopolymer(i, l)
        if homopolymer:
            homopolymer_array.append(i)
    return homopolymer_array

def retrieve_homopolymer_stats(filename):
    """
    """
    byte_to_format = {1: 'B', 2: 'H', 4: 'I', 8: 'Q'}
    found_homopolymers = []
    with open(filename, 'rb') as f:
        # Skip header
        f.seek(6)  # Skip magic(4) + version(2)
        l = struct.unpack('<H', f.read(2))[0]
        l = int(l)
        byte_size = struct.unpack('<B', f.read(1))[0]
        byte_format = byte_to_format[byte_size]
        homopolymers = generate_homopolymer_array(l)
        try:
            while True:
                index_bytes = f.read(byte_size)
                if len(index_bytes) < byte_size:
                    break
                v1 = struct.unpack(f'<{byte_format}', index_bytes)[0]
                nullomer_count_bytes = f.read(byte_size)
                if len(nullomer_count_bytes) < byte_size:
                    break
                nullomer_count = struct.unpack(f'<{byte_format}', nullomer_count_bytes)[0]
                total_bytes = nullomer_count * byte_size
                i = 0
                if v1 in homopolymers:
                    while i < nullomer_count:
                        nullomer_byte = f.read(byte_size)
                        total_bytes -= byte_size
                        v2_index = struct.unpack(f'<{byte_format}', nullomer_byte)[0]
                        if v2_index == v1:
                            found_homopolymers.append(v1)
                            f.seek(total_bytes,1)
                            break
                        i += 1
                else:
                    f.seek(total_bytes,1)
        except (struct.error, OSError):
            pass
    return found_homopolymers

# Output analysis management

#----------------------------------------------------------------------#
# Funções em comum
#----------------------------------------------------------------------#

def checar_genoma(entrada):
    """
    Checa se o arquivo FASTA contém apenas as bases válidas (A, T, C, G).
    Retorna True se o arquivo está correto, False caso contrário.
    """
    bases_validas = set('ATCG')
    with open(entrada, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            if not set(line).issubset(bases_validas):
                # Encontrou base inválida
                return False
    return True

def eh_palindromo(sequence): # checa se a sequência é palíndroma

    for i in range(len(sequence) // 2):
        # Compara o caractere na posição i com o caractere simétrico
        if sequence[i] != sequence[-(i + 1)]:
            return False
    return True

def ler_arquivo_caixa_alta(caminho_arquivo): # retorna o arquivo em caixa alta para leitura adequada do genoma
    with open(caminho_arquivo, 'r') as file:
        conteudo = file.read()
    conteudo_caixa_alta = conteudo.upper()
    with open(caminho_arquivo, 'w') as file:
        file.write(conteudo_caixa_alta)

def ler_organismos(caminho_arquivo): # retorna lista de organismos que serão analisados
    with open(caminho_arquivo, 'r') as file:
        organismos = [linha.strip() for linha in file.readlines()]  # Remove espaços e quebras de linha
    return organismos

def calcular_cpg(sequencia):
    """
    Conta o número de dinucleotídeos CG em uma sequência.
    :param sequencia: Sequência de nucleotídeos.
    :return: Número de dinucleotídeos CG.
    """
    cpg_count = 0
    for i in range(len(sequencia) - 1):
        if sequencia[i] == 'C' and sequencia[i + 1] == 'G':
            cpg_count += 1
    return cpg_count

