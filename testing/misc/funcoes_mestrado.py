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
from collections import Counter
import time
import subprocess
import psutil
from memory_profiler import memory_usage

# ---------------------------------------------------------------------------------------------------------------#
# Funções em teste
# ---------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------# 
# Funções ocasionais
# ---------------------------------------------------------------------------------------------------------------#  

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

def criar_dic_reversos(l):
    """
    Cria um dicionário onde cada índice lexicográfico está relacionado ao índice
    lexicográfico da sua sequência reversa, sem duplicidade (ou seja, só um dos pares é registrado).
    Exemplo: se ATCG (i) pareia com GCTA (j), só (i: j) estará no dicionário, não (j: i).

    :param l: Comprimento da sequência.
    :return: Dicionário {indice_lexico: indice_lexico_reverso}
    """
    dic_reversos = {}
    for i in range(4**l):
        seq = indice_para_seq(i, l)
        seq_rev = seq[::-1]
        j = calc_ind_lexicografico(seq_rev, l)
        dic_reversos[i] = j
    return dic_reversos

def indices_homopolimeros(l):
    """
    Retorna um array com os índices lexicográficos das sequências homopolímeras (AAAA..., TTTT..., CCCC..., GGGG...).
    :param l: Comprimento da sequência.
    :return: Lista de índices lexicográficos.
    """
    bases = ['A', 'T', 'C', 'G']
    return [calc_ind_lexicografico(base * l, l) for base in bases]

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

# ---------------------------------------------------------------------------------------------------------------#
# Abordagem via TrieBit direto do genoma
# ---------------------------------------------------------------------------------------------------------------#

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

def contar_nulomeros_trie_bit_novo(trie_bit, l):
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

def importar_triebit_teste_txt(arquivo_entrada, l, m):
    """
    Lê um arquivo no formato exportado por escrever_trie_em_txt_bitarray e reconstrói uma TrieBitTeste.
    :param arquivo_entrada: Caminho do arquivo .txt.
    :param l: Comprimento de v1.
    :param m: Comprimento de v2 (tamanho do bitarray).
    :return: Instância de TrieBitTeste reconstruída.
    """
    trie = TrieBitTeste(m)
    with open(arquivo_entrada, 'r') as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None and line:
                    # Converte o índice lexicográfico para o caminho v1
                    v1 = []
                    temp = v1_lexico
                    for _ in range(l):
                        div = 4 ** (l - len(v1) - 1)
                        v1.append(temp // div)
                        temp = temp % div
                    for v2 in line.split(','):
                        if v2:
                            trie.insert(v1, int(v2))
    return trie

def media_gc_por_organismo_triebit_txt(arquivo_txt, l, k):
    """
    Calcula a média de GC% para cada organismo a partir de um arquivo trie_bit_txt_novo.
    :param arquivo_txt: Caminho do arquivo .txt.
    :param l: Comprimento de v1 e v2.
    :param k: Tamanho total do k-mer.
    :return: Float com a média de GC%, 4 casas decimais.
    """
    soma_gc = 0
    total_pares = 0
    gc_dict = precomputar_gc_cpg(l)  # {indice_lexico_v1: [total_gc, total_seqs]}
    with open(arquivo_txt, 'r') as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None and line:
                    v2_indices = [int(x) for x in line.split(',') if x]
                    for v2 in v2_indices:
                        soma_gc += gc_dict[v1_lexico][0] + gc_dict[v2][0]  # Soma GC% de v1 e v2
                        total_pares += 1
    media_gc = media_gc = (soma_gc / (total_pares * k)) * 100 if total_pares > 0 else 0

    return round(media_gc,4)

def total_palindromo_organismo_triebit_txt(arquivo_txt, l):
    """
    Conta o total de palíndromos em um arquivo trie_bit_txt_novo.
    :param arquivo_txt: Caminho do arquivo .txt.
    :param l: Comprimento de v1 e v2.
    :return: Inteiro com o total de palíndromos.
    """
    ind_reversos = criar_dic_reversos(l)
    total_palindromos = 0
    with open(arquivo_txt, 'r') as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None and line:
                    v1 = ind_reversos[v1_lexico]  # Obtém o índice reverso de v1
                    v2_indices = set(int(x) for x in line.split(',') if x)
                    if v1 in v2_indices:
                        total_palindromos += 1
    return total_palindromos

def total_homopolimeros_organismo_triebit_txt(arquivo_txt, l):
    """
    Conta o total de homopolímeros em um arquivo trie_bit_txt_novo.
    :param arquivo_txt: Caminho do arquivo .txt.
    :param l: Comprimento de v1 e v2.
    :return: Inteiro com o total de homopolímeros.
    """
    homopolimeros = indices_homopolimeros(l)
    total_homopolimeros = 0
    with open(arquivo_txt, 'r') as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None and line:
                    if v1_lexico in homopolimeros:
                        v2_indices = set(int(x) for x in line.split(',') if x)
                        if v1_lexico in v2_indices:
                            total_homopolimeros += 1
    return total_homopolimeros

def total_null_txt(arquivo_txt, l):
    """
    Conta o total de nulômeros em um arquivo trie_bit_txt_novo.
    :param arquivo_txt: Caminho do arquivo .txt.
    :param l: Comprimento de v1 e v2.
    :return: Inteiro com o total de nulômeros.
    """
    total_nulomeros = 0
    with open(arquivo_txt, 'r') as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None and line:
                    v2_indices = set(int(x) for x in line.split(',') if x)
                    total_nulomeros += len(v2_indices)  # Conta cada v2 como um nulômero
    return total_nulomeros

def total_cpg_organismo_triebit_txt(arquivo_txt, l):
    cpg_dict, cg_dict = precomputar_cpg_terminaC_comecaG(l)
    total_cpg = 0
    with open(arquivo_txt) as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None:
                    v2_indices = [int(x) for x in line.split(',') if x]
                    for v2 in v2_indices:
                        total_cpg += cpg_dict[v1_lexico] + cpg_dict[v2]
                        if cg_dict[v1_lexico][0] and cg_dict[v2][1]:
                            total_cpg += 1
    return total_cpg

def total_cpg_triebit_por_cpg(arquivo_txt, l):
    cpg_dict, cg_dict = precomputar_cpg_terminaC_comecaG(l)
    cpg_counter = Counter()
    with open(arquivo_txt) as f:
        v1_lexico = None
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                v1_lexico = int(line[1:])
            else:
                if v1_lexico is not None:
                    total_cpg = 0
                    v2_indices = [int(x) for x in line.split(',') if x]
                    for v2 in v2_indices:
                        total_cpg += cpg_dict[v1_lexico] + cpg_dict[v2]
                        if cg_dict[v1_lexico][0] and cg_dict[v2][1]:
                            total_cpg += 1
                        cpg_counter[total_cpg] += 1
                        total_cpg = 0
    return dict(cpg_counter)

def parametros_cpg(arquivo_txt, l, total_null):
    total_cpg= total_cpg_organismo_triebit_txt(arquivo_txt, l)
    if total_null != 0:
        cpg_relativo_total = (total_cpg * 100) / total_null
    else:
        cpg_relativo_total = 0
    cpg_counter = total_cpg_triebit_por_cpg(arquivo_txt, l)
    return {
        "total_cpg": total_cpg,
        "cpg_total": cpg_relativo_total,
        "cpg_counter": cpg_counter
    }

def parametros_nulomeros(arquivo_txt, l):
    """
    Calcula os parâmetros de nulômeros a partir de um arquivo trie_bit_txt_novo.
    :param arquivo_txt: Caminho do arquivo .txt.
    :param l: Comprimento de v1 e v2.
    :return: Dicionário com os parâmetros calculados.
    """
    total_nulomeros = total_null_txt(arquivo_txt, l)
    total_palindromos = total_palindromo_organismo_triebit_txt(arquivo_txt, l)
    total_homopolimeros = total_homopolimeros_organismo_triebit_txt(arquivo_txt, l)
    total_gc = media_gc_por_organismo_triebit_txt(arquivo_txt, l, 2*l)  # k = 2*l para o cálculo de GC%
    total_cpg = parametros_cpg(arquivo_txt, l, total_nulomeros)
    return {
        "total_nulomeros": total_nulomeros,
        "total_palindromos": total_palindromos,
        "total_homopolimeros": total_homopolimeros,
        "total_gc": total_gc,
        "total_cpg": total_cpg["total_cpg"],
        "cpg_relativo": total_cpg["cpg_total"],
        "cpg_counter": total_cpg["cpg_counter"]
    }

def processar_null_batch(
        csv_saida,
        base_path,
        config_file
):
    """
    Recebe múltiplos arquivos de texto contendo tries compactadas,
    Remonta elas em memória e conta o total de nulômeros para cada,
    retornando o total de nulômeros encontrados em um arquivo csv.
    :param csv_saida: especifica caminho para gerar resultados.
    :param base_path: especifica caminho base para os arquivos de tries.
    :param k_values: lista de valores de k para os quais as tries foram geradas, referente às subpastas do caminho base_path.
    :param config_file: caminho para o arquivo de configuração YAML contendo a lista de organismos.
    """
    orgs_bruto = obter_lista_org(config_file)
    orgs = [ajustar_nome_orgs(org) for org in orgs_bruto]
    erros = []
    resultados = []
    k_values = obter_ks_analisados(base_path)
    for k in k_values:
        l = int(k // 2)
        k_path = base_path + str(k) + '/'
        if not os.path.exists(k_path):
            print(f"Diretório {k_path} não existe.")
        cont = 0
        for org in orgs:
            org_path = k_path + org + '/'
            arquivo_txt = org_path + f'nulomerostrie_{org}_{k}.txt'
            if not os.path.exists(arquivo_txt):
                print(f"Arquivo {arquivo_txt} não existe. Pulando.")
                cont += 1
                erros.append([org, k, 'Arquivo não encontrado'])
                continue

            if not os.path.exists(org_path):
                print(f"Diretório {org_path} não existe.")
                erros.append([org, k, 'Diretório não encontrado'])
            else:
                cont += 1
                print(f"Processando organismo {org}. {cont} de um total de {len(orgs)} organismos para k = {k}.")
                org_path = org_path + 'nulomerostrie_' + org + '_' + str(k) + '.txt'
                parametros = parametros_nulomeros(org_path, l)
                print(f"Total de nulômeros encontrados para {org} com k={k}: {parametros['total_nulomeros']}")
                linha = {'Organismo': org, 'k': k}
                for chave, valor in parametros.items():
                    linha[chave] = valor
                resultados.append(linha)

    if resultados:
        colunas = ['Organismo', 'k'] + [k for k in resultados[0].keys() if k not in ['Organismo', 'k']]
        with open(csv_saida, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=colunas)
            writer.writeheader()
            writer.writerows(resultados)
        # writer.writerows(media_gc_lista)
    with open(base_path + 'erros_nulomeros.txt', 'w') as f:
        for erro in erros:
            f.write(f"{erro[0]}, {erro[1]}, {erro[2]}\n")
    print(f"Organismos não encontrados: {erros}")

# ---------------------------------------------------------------------------------------------------------------#
# Verificação de possíveis nulômeros (k-mers que podem tornar-se nulômeros por mutação de ponto)
# ---------------------------------------------------------------------------------------------------------------#

def medir_memoria_worker(args):
    return max(memory_usage((processa_organismo, (args,)), interval=0.1))

def calcula_workers_ajustado(mem_por_worker_gb=2, max_workers_user=4):
    # Memória total disponível em GB
    mem_total_gb = psutil.virtual_memory().available / (1024 ** 3)
    # Calcula quantos workers cabem, deixando 2GB livres para o sistema
    max_workers_mem = max(1, int((mem_total_gb - 2) // mem_por_worker_gb))
    # Usa o menor entre o pedido pelo usuário e o limite pela RAM
    return min(max_workers_user, max_workers_mem)

def processa_organismo(args):
    genoma, arquivo_txt, l, m, k = args
    trie = importar_triebit_teste_txt(arquivo_txt, l, m)
    proc = subprocess.Popen(
        ['/home/leveduras/integranteslab/matheus/Mestradoteste/./fasta_kmers_novo', genoma, str(k)],
        stdout=subprocess.PIPE,
        text=True
    )
    results = []
    for line in proc.stdout:
        for kmer in line.strip().split():
            kmer_indices = kmer_str_para_indices(kmer)
            min_mut = busca_min_mutacoes(trie, kmer_indices, 3)
            if min_mut is not None:
                results.append(("".join(kmer), min_mut))
    proc.wait()
    return results

def kmer_str_para_indices(kmer_str):
    return [int(x) for x in kmer_str.split(',')]

def busca_min_mutacoes(trie, kmer, max_mut):
    min_mut = None
    def dfs(node, pos, mismatches):
        nonlocal min_mut
        if mismatches > max_mut:
            return
        if pos == len(kmer):
            for v2 in range(len(node.v2_set)):
                if node.v2_set[v2]:
                    if min_mut is None or mismatches < min_mut:
                        min_mut = mismatches
            return
        base = kmer[pos]
        for i in range(4):
            if node.children[i] is not None:
                dfs(node.children[i], pos+1, mismatches + (i != base))
    dfs(trie.root, 0, 0)
    return min_mut

# ----------------------------------------------------------------------------------------------------------------#
# Funções em comum
# ----------------------------------------------------------------------------------------------------------------#

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

def calcular_gc(sequencia):
    """
    Calcula a porcentagem de conteúdo GC em uma sequência, formatada com 4 casas decimais.
    :param sequencia: Sequência de nucleotídeos.
    :return: Porcentagem de GC com 4 casas decimais.
    """
    gc_count = sequencia.count('G') + sequencia.count('C')
    return round((gc_count / len(sequencia)) * 100, 1)

def calcular_gc_total(sequencia):
    """
    Calcula a porcentagem de conteúdo GC em uma sequência, formatada com 4 casas decimais.
    :param sequencia: Sequência de nucleotídeos.
    :return: Porcentagem de GC com 4 casas decimais.
    """
    gc_count = sequencia.count('G') + sequencia.count('C')
    return gc_count

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

# ---------------------------------------------------------------------------------------------------------------#
# Funções "aposentadas"
# ---------------------------------------------------------------------------------------------------------------#