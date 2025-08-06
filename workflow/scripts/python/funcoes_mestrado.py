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

def contar_seq_abordagem_antiga(arquivo):
    """
    Conta o número total de sequências em um arquivo, onde cada sequência está em uma linha separada.
    :param arquivo: Caminho para o arquivo.
    :return: Número total de sequências.
    """
    with open(arquivo, 'r') as f:
        total_sequencias = sum(1 for linha in f if linha.strip())  # Conta apenas linhas não vazias
    return total_sequencias

def ajustar_nome_orgs(nome_org):
    nome_org = nome_org.replace(" ", "_")  # Substitui espaços por underscores
    return re.sub(r'[^a-zA-Z0-9_]', '_', nome_org)  # Remove caracteres especiais, mantendo apenas letras, números e underscores

#def checar_maw_final

def checar_maw_inicial(trie, caminho):

    node = trie.root
    for base in caminho:
        if node.children[base] is None:
            return False
        node = node.children[base]
    return True

def buscar_trie_menor(trie_k,trie_kmenor,k):

    resultados = []
    def dfs(node, path):
        print(f"DFS: {path}")
        if len(path) == k:
            caminho_kmenor = path[1:]  # remove a primeira base
            existe = checar_maw_inicial(trie_kmenor, caminho_kmenor)
            print(f"Caminho: {path}, Caminho k menor: {caminho_kmenor}, Existe: {existe}")
            resultados.append((tuple(path), tuple(caminho_kmenor), existe))
        for base, child in enumerate(node.children):
            if child is not None:
                dfs(child, path + [base])

    dfs(trie_k.root, [])
    return resultados

def normalizar_nome(nome):
    nome = nome.lower().strip()
    nome = re.sub(r'[^a-z0-9]', '_', nome)
    return nome

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

def precomputar_cpg_terminaC_comecaG(l):
    base_map = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    cpg_dict = {}
    cg_dict = {}
    for idx in range(4**l):
        seq = indice_para_seq(idx, l)
        cpg_dict[idx] = calcular_cpg(seq)
        cg_dict[idx] = (seq.endswith('C'), seq.startswith('G'))
    return cpg_dict, cg_dict

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
# Funções voltadas ao uso em database
# ---------------------------------------------------------------------------------------------------------------#

def gerar_total_null_csv(dicio_dados, arquivo_saida, modo = 'w'):
    """
    Gera um arquivo CSV a partir de um dicionário de dados.
    :param dicio_dados: Dicionário contendo os dados a serem salvos.
    :param arquivo_saida: Nome do arquivo CSV de saída.
    """
    with open(arquivo_saida, modo, newline='') as csvfile:
        fieldnames = ['Sequência', 'gc', 'CpG', 'Termina C', 'Comeca G', 'Sequence', 'l_value']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        if modo == 'w':
            writer.writeheader()
        for sequencia, dados in dicio_dados.items():
            writer.writerow({'Sequência': sequencia, **dados})

def criar_dic_organismo_id(arquivo_csv):

    valor = "ID"
    chave = "Organism Name"
    chave_valor = {}
    with open(arquivo_csv, 'r') as f:
        leitor = csv.DictReader(f)
        for linha in leitor:
            chave_valor[linha[chave]] = linha[valor]
    return chave_valor

def adicionar_org_id(db_org, results_org, org,k):
    dic = criar_dic_organismo_id(db_org)
    id = dic[org]
    # Adicionar coluna de ID de organismo ao arquivo CSV
    with open(results_org, 'r') as f:
        leitor = csv.reader(f)
        cabecalho = next(leitor)
        linhas = list(leitor)

    cabecalho.append("organismo_id")
    cabecalho.append("k_value")
    for linha in linhas:
        linha.append(id)
        linha.append(k)

    with open(results_org, 'w', newline='') as f:
        escritor = csv.writer(f)
        escritor.writerow(cabecalho)
        escritor.writerows(linhas)
    print(f"Coluna de ID de organismo adicionada ao arquivo {results_org} com sucesso!")

def importar_csv_para_postgresql(caminho_csv, tabela, conexao):
    """
    Importa um arquivo CSV para uma tabela no PostgreSQL.
    :param caminho_csv: Caminho para o arquivo CSV.
    :param tabela: Nome da tabela no PostgreSQL.
    :param conexao: Conexão com o banco de dados PostgreSQL.
    """
    with conexao.cursor() as cursor:
        with open(caminho_csv, 'r') as f:
            query = f"""
            COPY {tabela}(v1,v2,ind,org_esp)
            FROM STDIN WITH CSV HEADER DELIMITER ','
            """
            cursor.copy_expert(query, f)
    conexao.commit()

def obter_lista_org(caminho_yaml):
    with open(caminho_yaml, 'r') as f:
        dados = yaml.safe_load(f)
        organismos = dados.get('organismos', [])
        return organismos

def criar_dados_nulomeros(l):
    dicio_dados ={}
    sequencias = gerar_kmers(l)

    for sequencia in sequencias:
        lexi = calc_ind_lexicografico(sequencia, l)
        seqg,seqc = checar_inicio_fim(sequencia)
        dicio_dados[lexi] = {
            "gc": calcular_gc_total(sequencia),
            "CpG": calcular_cpg(sequencia),
            "Termina C": seqc,
            "Comeca G": seqg,
            'Sequence': sequencia,
            'l_value': l
        }
    return dicio_dados

def criar_relacao_k_l(k_list):

    k_l = {}

    for k in k_list:
        if k % 2 != 0:
            l1 = int(k / 2)
            l2 = int(k / 2) + 1
            k_l[k] = [l1, l2]
        else:
            l = int(k / 2)
            k_l[k] = [l]

    return k_l

# ---------------------------------------------------------------------------------------------------------------#
# Abordagem via Tries + sets
# ---------------------------------------------------------------------------------------------------------------#

class TrieNode:
    def __init__(self):
        self.children = [None] * 4  # Lista fixa para as bases A, T, C, G
        self.v2_set = set()  # Conjunto para armazenar os valores de v2
        self.v2_complementar_set = set()  # Conjunto para armazenar os valores complementares de v2
        self.v2_duplicado_set = set()  # Conjunto para armazenar os valores duplicados de v2

class Trie:

    def __init__(self):
        self.root = TrieNode()

    def insert(self, v1, v2):
        """
        Insere `v1` na Trie e adiciona `v2` ao conjunto associado ao último nó de `v1`.
        :param v1: Lista de inteiros representando o caminho na Trie.
        :param v2: Valor inteiro a ser adicionado ao conjunto.
        """
        node = self.root
        for value in v1:
            if node.children[value] is None:
                node.children[value] = TrieNode()
            node = node.children[value]
        # Adiciona o valor `v2` ao conjunto
        node.v2_set.add(v2)

    def find_missing_v2(self, v1, m):
        """
        Retorna os valores de v2 que estão ausentes para um caminho v1.
        :param v1: Lista de inteiros representando o caminho na Trie.
        :param m: Valor máximo de v2 (4^l - 1).
        :return: Lista de valores ausentes.
        """
        node = self.root
        for value in v1:
            if node.children[value] is None:
                return list(range(m))  # Se o caminho não existir, todos os valores estão ausentes
            node = node.children[value]
        # Calcula os valores ausentes comparando o conjunto com todos os índices possíveis
        ausentes = [i for i in range(m) if i not in node.v2_set]
        return ausentes

    def print_trie(self):
        """
        Imprime os caminhos da Trie e os valores de v2 associados.
        """
        def dfs(node, path):
            if node.v2_set:
                print(f"Caminho: {' -> '.join(map(str, path))}")
                print(f"Valores de v2: {sorted(node.v2_set)}")
            for child_value, child_node in node.children.items():
                dfs(child_node, path + [child_value])

        dfs(self.root, [])

def dfs_iterativo_limited(trie, l, operation):
    """
    Percorre todos os caminhos de comprimento `l` na Trie usando DFS iterativo
    e aplica uma operação nos sets anexados aos nós terminais.
    
    :param trie: Instância da Trie.
    :param l: Comprimento dos caminhos a serem considerados.
    :param operation: Função a ser aplicada nos sets dos nós terminais.
    """
    stack = [(trie.root, [])]  # Pilha para armazenar (nó atual, caminho atual)

    while stack:
        node, path = stack.pop()

        # Se o caminho atual tiver comprimento `l`, aplica a operação no set do nó terminal
        if len(path) == l:
            if node.v2_set:
                operation(path, node.v2_set)
            continue  # Não expande mais este caminho

        # Adiciona os filhos do nó atual à pilha
        for child_value, child_node in enumerate(node.children):  # Corrigido para iterar sobre a lista
            if child_node is not None:  # Verifica se o filho existe
                stack.append((child_node, path + [child_value]))

def print_trie_completa(trie):
    """
    Imprime os caminhos da Trie e os valores associados:
    - v2_set (valores originais)
    - v2_complementar_set (valores complementares)
    - v2_duplicado_set (valores duplicados)
    """
    def dfs(node, path):
        # Verifica se o nó atual possui valores em qualquer conjunto
        if node.v2_set or node.v2_complementar_set or node.v2_duplicado_set:
            caminho_formatado = "-".join(map(str, path))  # Formata o caminho como uma string
            
            # Formatar os conjuntos
            valores_originais = ",".join(map(str, sorted(node.v2_set))) if node.v2_set else "[]"
            valores_complementares = ",".join(map(str, sorted(node.v2_complementar_set))) if node.v2_complementar_set else "[]"
            valores_duplicados = ",".join(map(str, sorted(node.v2_duplicado_set))) if node.v2_duplicado_set else "[]"
            
            # Imprimir os dados
            print(f"Caminho: {caminho_formatado}")
            print(f"  Valores originais de v2: {valores_originais}")
            print(f"  Valores complementares de v2: {valores_complementares}")
            print(f"  Valores duplicados de v2: {valores_duplicados}")
            print("-" * 50)  # Separador para melhor visualização

        # Percorrer os filhos do nó atual
        for child_value, child_node in enumerate(node.children):
            if child_node is not None:  # Verifica se o filho existe
                dfs(child_node, path + [child_value])

    # Inicia a busca a partir da raiz da Trie
    dfs(trie.root, [])

def adicionar_caminhos_complementares(trie, imp, l):
    """
    Adiciona os caminhos complementares e os valores de v2 complementares na mesma Trie.
    Se o valor complementar já existir no caminho complementar, ele é movido para o conjunto de duplicados.
    """
    complemento_base = {0:1, 1:0, 2:3, 3:2}  # A=0 ↔ T=1, C=2 ↔ G=3
    
    def calcular_complemento_caminho(caminho):
        return [complemento_base[base] for base in reversed(caminho)]

    def abrir_v2(v2, l):
        bases = []
        for _ in range(l):
            bases.append(v2 % 4)
            v2 //= 4
        return bases
        # bases_complementares = [complemento_base[base] for base in reversed(bases)]
        # indice_complementar = 0
        # for i, base in enumerate(bases_complementares):
        #     indice_complementar += base * (4 ** (l - i - 1))
        # return indice_complementar
    def ind_v2_comp(v2, l):
        indice_complementar = 0
        for i, base in enumerate(v2):
            indice_complementar += base * (4 ** (l - i - 1))
        return indice_complementar

    def processar_caminho(path, v2_set):
        for v2 in list(v2_set):  # Cria uma cópia do conjunto para evitar modificações durante a iteração
            print(f"Processando caminho: {path} com v2: {v2}")
            caminho_completo = calcular_complemento_caminho((path + abrir_v2(v2, l)))
            print(f"Caminho completo: {caminho_completo}")
            v1 = caminho_completo[:l]  # O caminho v1 é o primeiro l elementos do caminho completo
            v2_complementar = ind_v2_comp(caminho_completo[l:],l)  # O valor v2 é o próximo elemento após o caminho v1
            print(f"v1: {v1}, v2: {v2}")
            caminho_complementar = calcular_complemento_caminho(v1)
            print(f"v2 complementar: {v2_complementar}")
            node = trie.root
            for value in caminho_complementar:
                if node.children[value] is None:
                    node.children[value] = TrieNode()
                node = node.children[value]
            
            # Se o valor complementar já existir no caminho complementar
            if v2_complementar in node.v2_set:
                
                # Move o valor complementar para o conjunto de duplicados no caminho complementar
                node.v2_set.discard(v2_complementar)
                node.v2_duplicado_set.add(v2_complementar)
                
                # Move o valor original para o conjunto de duplicados no caminho original
                original_node = trie.root
                for value in path:
                    original_node = original_node.children[value]
                original_node.v2_set.discard(v2)
                original_node.v2_duplicado_set.add(v2)
            
            # Caso contrário, adiciona o valor complementar ao conjunto complementar
            else:
                node.v2_complementar_set.add(v2_complementar)
    if imp:
        l_dfs = l - 1
    else:
        l_dfs = l
    dfs_iterativo_limited(trie, l_dfs, processar_caminho)

def criar_trie_original(genoma, k, l, base_values, complemento_base):
    matches_o = Trie()
    base_to_int = base_values
    comp_table = str.maketrans('ATCG', 'TAGC')
    insert = matches_o.insert

    for kmer in sliding_window_kmers_fasta(genoma, k):
        if any(base not in base_to_int for base in kmer):
            continue

        v1 = tuple(base_to_int[b] for b in kmer[:l])
        v2_idx = 0
        for b in kmer[l:]:
            v2_idx = (v2_idx << 2) | base_to_int[b]
        insert(v1, v2_idx)

        # Complementar reverso
        kmer_comp = kmer.translate(comp_table)[::-1]
        v1c = tuple(base_to_int[b] for b in kmer_comp[:l])
        v2c_idx = 0
        for b in kmer_comp[l:]:
            v2c_idx = (v2c_idx << 2) | base_to_int[b]
        insert(v1c, v2c_idx)

    return matches_o

def criar_trie_complementar(trie_original, l, m):
    """
    Cria uma Trie complementar a partir de uma Trie original, reorganizando os conjuntos:
    - Todo conjunto complementar (v2_complementar_set) é escrito como conjunto original (v2_set).
    - Todo conjunto original (v2_set) é escrito como conjunto complementar (v2_complementar_set).
    - Os índices que não aparecem em nenhum dos conjuntos são escritos como duplicados (v2_duplicado_set).

    :param trie_original: Instância da Trie original.
    :param l: Comprimento dos caminhos (tamanho de v1).
    :param m: Valor máximo de v2 (4**l - 1).
    :return: Nova Trie complementar.
    """
    # Criar uma nova Trie complementar
    trie_complementar = Trie()

    def processar_caminho(path, node):
        """
        Processa cada caminho da Trie original e reorganiza os conjuntos na Trie complementar.
        """
        # Calcular todos os índices possíveis de v2
        todos_os_indices = set(range(m))

        # Reorganizar os conjuntos
        novo_v2_set = node.v2_complementar_set  # Complementar vira Original
        novo_v2_complementar_set = node.v2_set  # Original vira Complementar
        novo_v2_duplicado_set = todos_os_indices - (node.v2_set | node.v2_complementar_set | node.v2_duplicado_set)

        # Navegar até o nó correspondente na Trie complementar
        node_complementar = trie_complementar.root
        for value in path:
            if node_complementar.children[value] is None:
                node_complementar.children[value] = TrieNode()
            node_complementar = node_complementar.children[value]

        # Adicionar os conjuntos reorganizados ao nó da Trie complementar
        node_complementar.v2_set.update(novo_v2_set)
        node_complementar.v2_complementar_set.update(novo_v2_complementar_set)
        node_complementar.v2_duplicado_set.update(novo_v2_duplicado_set)

    # Percorre todos os caminhos da Trie original e processa os valores
    def dfs(node, path):
        if node.v2_set or node.v2_complementar_set or node.v2_duplicado_set:
            processar_caminho(path, node)
        for child_value, child_node in enumerate(node.children):  # Corrigido para usar enumerate()
            if child_node is not None:
                dfs(child_node, path + [child_value])

    # Iniciar a busca a partir da raiz da Trie original
    dfs(trie_original.root, [])

    return trie_complementar

def escrever_trie_em_txt(trie, arquivo_saida):
    """
    Escreve os caminhos e valores de uma Trie em um arquivo de texto no formato:
    caminho separado por "-" (exemplo: 0-0-0-0):
    - Valores originais de v2 separados por "," (exemplo: 23,26,27)
    - Valores complementares de v2 separados por "," (exemplo: 30,31)
    - Valores duplicados de v2 separados por "," (exemplo: 15,16)
    """
    with open(arquivo_saida, 'w') as f:
        def dfs(node, path):
            if node.v2_set or node.v2_complementar_set or node.v2_duplicado_set:
                caminho_formatado = "-".join(map(str, path))
                f.write(f"Caminho: {caminho_formatado}\n")
                if node.v2_set:
                    valores_originais = ",".join(map(str, sorted(node.v2_set)))
                    f.write(f"  Valores originais de v2: {valores_originais}\n")
                if node.v2_complementar_set:
                    valores_complementares = ",".join(map(str, sorted(node.v2_complementar_set)))
                    f.write(f"  Valores complementares de v2: {valores_complementares}\n")
                if node.v2_duplicado_set:
                    valores_duplicados = ",".join(map(str, sorted(node.v2_duplicado_set)))
                    f.write(f"  Valores duplicados de v2: {valores_duplicados}\n")
            for child_value, child_node in enumerate(node.children):
                if child_node is not None:
                    dfs(child_node, path + [child_value])

        dfs(trie.root, [])

def contar_caminhos(trie):
    """
    Conta o número de caminhos terminais na Trie.
    Um caminho terminal é aquele que possui valores associados no conjunto v2_set.

    :param trie: Instância da Trie.
    :return: Número total de caminhos terminais.
    """
    def dfs(node):
        count = 1 if node.v2_set else 0  # Conta o nó atual se tiver valores associados
        for child_node in node.children.values():
            count += dfs(child_node)
        return count

    return dfs(trie.root)

def print_todos_os_sets(trie):
    """
    Imprime todos os caminhos da Trie e os valores associados a cada conjunto:
    - v2_set (valores originais)
    - v2_complementar_set (valores complementares)
    - v2_duplicado_set (valores duplicados)
    """
    def dfs(node, path):
        if node.v2_set or node.v2_complementar_set or node.v2_duplicado_set:
            # Formatar o caminho
            caminho_formatado = "-".join(map(str, path))
            
            # Formatar os conjuntos
            valores_originais = ",".join(map(str, sorted(node.v2_set))) if node.v2_set else "[]"
            valores_complementares = ",".join(map(str, sorted(node.v2_complementar_set))) if node.v2_complementar_set else "[]"
            valores_duplicados = ",".join(map(str, sorted(node.v2_duplicado_set))) if node.v2_duplicado_set else "[]"
            
            # Imprimir os dados
            print(f"Caminho: {caminho_formatado}")
            print(f"  Valores originais de v2: {valores_originais}")
            print(f"  Valores complementares de v2: {valores_complementares}")
            print(f"  Valores duplicados de v2: {valores_duplicados}")
            print("-" * 50)  # Separador para melhor visualização
        
        # Percorrer os filhos do nó atual
        for child_value, child_node in node.children.items():
            dfs(child_node, path + [child_value])

    # Iniciar a busca a partir da raiz da Trie
    dfs(trie.root, [])

# ---------------------------------------------------------------------------------------------------------------#
# Abordagem via Tries + bitarray
# ---------------------------------------------------------------------------------------------------------------#

class TrieNodeBit:
    def __init__(self, m):
        self.children = [None] * 4  # Lista fixa para as bases A, T, C, G
        self.v2_set = bitarray(m)  # Array de bits para representar v2_set
        self.v2_set.setall(0)      # Inicializa todos os bits como 0
        self.v2_complementar_set = bitarray(m)
        self.v2_complementar_set.setall(0)
        self.v2_duplicado_set = bitarray(m)
        self.v2_duplicado_set.setall(0)

class TrieBit:
    def __init__(self, m):
        self.root = TrieNodeBit(m)  # Passa `m` para o nó raiz
    def iterate(self, path=None):
        if path is None:
            path = []
        # Retorna o caminho atual e o nó
        yield path, self
        # Percorre os filhos
        for i, child in enumerate(self.children):
            if child is not None:
                yield from child.iterate(path + [i])

    def insert(self, v1, v2):
        """
        Percorre os nós definidos por v1 e seta o bit v2 no nó terminal.
        """
        node = self.root
        for value in v1:
            if node.children[value] is None:
                node.children[value] = TrieNodeBit(len(node.v2_set))  # Usa o mesmo m
            node = node.children[value]
        node.v2_set[v2] = 1  # Seta o bit correspondente

def inicializar_triebit_completa(l, m):
    trie = TrieBit(m)
    def construir(node, depth):
        if depth == l:
            return
        for i in range(4):
            if node.children[i] is None:
                node.children[i] = TrieNodeBit(m)
            construir(node.children[i], depth + 1)
    construir(trie.root, 0)
    return trie

def nulomeros_trie_bit(trie_original, matches_n, m):
    """
    Cria uma Trie complementar a partir de uma Trie original, reorganizando os conjuntos:
    - Todo conjunto complementar (v2_complementar_set) é escrito como conjunto original (v2_set).
    - Todo conjunto original (v2_set) é escrito como conjunto complementar (v2_complementar_set).
    - Os índices que não aparecem em nenhum dos conjuntos são escritos como duplicados (v2_duplicado_set).

    :param trie_original: Instância da Trie original.
    :param m: Valor máximo de v2 (4**l - 1).
    :return: Nova Trie complementar, total de nulômeros.
    """
    trie_null = matches_n
    total_nulomeros = 0  # contador

    def processar_caminho(path, node):
        nonlocal total_nulomeros
        todos_os_indices = set(range(m))

        novo_v2_set = node.v2_complementar_set
        novo_v2_complementar_set = node.v2_set
        novo_v2_duplicado_set = todos_os_indices - (node.v2_set | node.v2_complementar_set | node.v2_duplicado_set)

        node_complementar = trie_null.root
        for value in path:
            if node_complementar.children[value] is None:
                node_complementar.children[value] = TrieNodeBit(m)
            node_complementar = node_complementar.children[value]

        # Adicionar os conjuntos reorganizados ao nó da Trie complementar e contar
        for index in novo_v2_set:
            idx = int(index)
            if not node_complementar.v2_set[idx]:
                node_complementar.v2_set[idx] = 1
                total_nulomeros += 1
        
        for index in novo_v2_complementar_set:
            idx = int(index)
            if not node_complementar.v2_complementar_set[idx]:
                node_complementar.v2_complementar_set[idx] = 1
                total_nulomeros += 1
        
        for index in novo_v2_duplicado_set:
            idx = int(index)
            if not node_complementar.v2_duplicado_set[idx]:
                node_complementar.v2_duplicado_set[idx] = 1
                total_nulomeros += 1

    def dfs(node, path):
        if node.v2_set or node.v2_complementar_set or node.v2_duplicado_set:
            processar_caminho(path, node)
        for child_value, child_node in enumerate(node.children):
            if child_node is not None:
                dfs(child_node, path + [child_value])

    dfs(trie_original.root, [])
    print(f"Total de nulômeros encontrados: {total_nulomeros}")
    return trie_null, total_nulomeros

def preencher_triebit_a_partir_trie(trie_original, trie_bit, l):
    """
    Marca na trie_bit todos os índices de v2 presentes nos sets da trie original.
    """
    def dfs(node_orig, node_bit, depth):
        if depth == l:
            # Marca os bits nos três sets
            for v2 in node_orig.v2_set:
                node_bit.v2_set[v2] = 1
            for v2 in node_orig.v2_complementar_set:
                node_bit.v2_complementar_set[v2] = 1
            for v2 in node_orig.v2_duplicado_set:
                node_bit.v2_duplicado_set[v2] = 1
            return
        for i in range(4):
            child_orig = node_orig.children[i] if node_orig.children[i] else None
            child_bit = node_bit.children[i] if node_bit.children[i] else None
            if child_orig is not None and child_bit is not None:
                dfs(child_orig, child_bit, depth + 1)
    dfs(trie_original.root, trie_bit.root, 0)

def contar_nulomeros_trie_bit_novo(trie_bit, l):
    """
    Conta o total de nulômeros em uma TrieBit.
    Um nulômero é cada (v1, v2) tal que v2 NÃO está setado em NENHUM dos três sets do nó terminal v1.
    """
    def dfs(node, depth):
        total = 0
        if depth == l:
            m = len(node.v2_set)
            for v2 in range(m):
                if not (node.v2_set[v2] or node.v2_complementar_set[v2] or node.v2_duplicado_set[v2]):
                    total += 1
            return total
        for child in node.children:
            if child is not None:
                total += dfs(child, depth + 1)
        return total
    return dfs(trie_bit.root, 0)

def print_trie_bit_nova(trie, l):
    """
    Imprime os caminhos da TrieBit e os índices dos bits setados em cada conjunto.
    :param trie: Instância de TrieBit.
    :param l: Comprimento do caminho (v1).
    """
    def dfs(node, path):
        if node.v2_set.any() or node.v2_complementar_set.any() or node.v2_duplicado_set.any():
            caminho_formatado = "-".join(map(str, path))
            print(f"Caminho: {caminho_formatado}")
            if node.v2_set.any():
                print("  v2_set:", [i for i, bit in enumerate(node.v2_set) if bit])
            if node.v2_complementar_set.any():
                print("  v2_complementar_set:", [i for i, bit in enumerate(node.v2_complementar_set) if bit])
            if node.v2_duplicado_set.any():
                print("  v2_duplicado_set:", [i for i, bit in enumerate(node.v2_duplicado_set) if bit])
            print("-" * 40)
        for child_value, child_node in enumerate(node.children):
            if child_node is not None:
                dfs(child_node, path + [child_value])
    dfs(trie.root, [])

def escrever_trie_em_txt_bitarray(trie, arquivo_saida, l):
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

def escrever_trie_em_csv_bitarray(trie, arquivo_saida, l):
    """
    Escreve os caminhos e valores de uma Trie baseada em bitarray em um arquivo CSV no formato:
    - v1: Índice lexicográfico de v1.
    - v2: Índices dos valores de v2.
    - Identificador: "+" para valores originais, "-" para complementares, "=" para duplicados.
    
    :param trie: Instância da Trie baseada em bitarray.
    :param arquivo_saida: Caminho para o arquivo CSV de saída.
    :param l: Comprimento de v1 (necessário para calcular o índice lexicográfico).
    """
    with open(arquivo_saida, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Escreve o cabeçalho
        writer.writerow(["v1", "v2", "Identificador"])

        def dfs(node, path):
            # Calcula o índice lexicográfico de v1
            indice_lexicografico = sum(base * (4 ** (l - i - 1)) for i, base in enumerate(path))
            
            # Processa os valores originais de v2
            for v2 in [i for i, bit in enumerate(node.v2_set) if bit]:
                writer.writerow([indice_lexicografico, v2, "+"])
            
            # Processa os valores complementares de v2
            for v2 in [i for i, bit in enumerate(node.v2_complementar_set) if bit]:
                writer.writerow([indice_lexicografico, v2, "-"])
            
            # Processa os valores duplicados de v2
            for v2 in [i for i, bit in enumerate(node.v2_duplicado_set) if bit]:
                writer.writerow([indice_lexicografico, v2, "="])

            # Percorre os filhos do nó atual
            for child_value, child_node in enumerate(node.children):
                if child_node is not None:
                    dfs(child_node, path + [child_value])

        # Inicia a busca a partir da raiz da Trie
        dfs(trie.root, [])

def escrever_dados_em_csv(dados, arquivo_saida):
    """
    Escreve o dicionário `dados` em um arquivo CSV, onde cada linha contém os valores de GC%, CpG e palíndromo.
    :param dados: Lista de tuplas ou estrutura similar contendo os valores (GC%, CpG, palíndromo).
    :param arquivo_saida: Caminho do arquivo CSV de saída.
    """
    with open(arquivo_saida, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Escreve o cabeçalho
        writer.writerow(["GC%", "CpG", "Palíndromo"])
        # Escreve os dados
        for gc, cpg, palindromo in dados:
            writer.writerow([gc, cpg, palindromo])

def calcular_dados_trie(matches_n, l, k):
    """
    Percorre a Trie `matches_n` e calcula GC%, CpG e palíndromo apenas para os últimos nós de cada caminho.
    
    :param matches_n: Instância da Trie.
    :param l: Comprimento dos caminhos (tamanho de v1).
    :param k: Comprimento total da sequência (v1 + v2).
    :return: Lista de tuplas contendo (GC%, CpG, palíndromo).
    """
    base_map = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    dados = []

    def converter_para_bases(indice, tamanho):
        """
        Converte um índice inteiro para uma sequência de bases no formato de string.
        """
        bases = []
        for _ in range(tamanho):
            bases.append(base_map[indice % 4])
            indice //= 4
        return ''.join(reversed(bases))

    def dfs(node, path):
        """
        Percorre a Trie e processa apenas os últimos nós de cada caminho.
        """
        # Verifica se o nó atual é um nó terminal (último nó do caminho)
        if all(child is None for child in node.children):  # Nó terminal
            if node.v2_set.any() or node.v2_complementar_set.any() or node.v2_duplicado_set.any():
                
                caminho_bases = ''.join(base_map[base] for base in path)  # Converte o caminho para bases
                
                # Processa os índices em `v2_set`
                for v2 in (i for i, bit in enumerate(node.v2_set) if bit):
                    v2_bases = converter_para_bases(v2, l)  # Converte v2 para bases
                    vector = caminho_bases + v2_bases  # Combina o caminho e v2
                    cpg = sum(1 for i in range(len(vector) - 1) if vector[i:i + 2] == 'CG')
                    gc = (sum(1 for char in vector if char in 'CG') * 100) / k
                    palindromo = eh_palindromo(vector)
                    dados.append((gc, cpg, palindromo))
                   

                # Processa os índices em `v2_complementar_set`
                for v2 in (i for i, bit in enumerate(node.v2_complementar_set) if bit):
                    v2_bases = converter_para_bases(v2, l)  # Converte v2 para bases
                    vector = caminho_bases + v2_bases  # Combina o caminho e v2
                    cpg = sum(1 for i in range(len(vector) - 1) if vector[i:i + 2] == 'CG')
                    gc = (sum(1 for char in vector if char in 'CG') * 100) / k
                    palindromo = eh_palindromo(vector)
                    dados.append((gc, cpg, palindromo))
                   

                # Processa os índices em `v2_duplicado_set`
                for v2 in (i for i, bit in enumerate(node.v2_duplicado_set) if bit):
                    v2_bases = converter_para_bases(v2, l)  # Converte v2 para bases
                    vector = caminho_bases + v2_bases  # Combina o caminho e v2
                    cpg = sum(1 for i in range(len(vector) - 1) if vector[i:i + 2] == 'CG')
                    gc = (sum(1 for char in vector if char in 'CG') * 100) / k
                    palindromo = eh_palindromo(vector)
                    dados.append((gc, cpg, palindromo))
                    
        # Continua a busca nos filhos
        for child_value, child_node in enumerate(node.children):
            if child_node is not None:
                dfs(child_node, path + [child_value])

    # Inicia a busca a partir da raiz da Trie
    dfs(matches_n.root, [])
    return dados

def print_trie_bit(trie, l):
    """
    Imprime os caminhos da TrieBit e os índices dos bits setados em cada conjunto.
    :param trie: Instância de TrieBit.
    :param l: Comprimento do caminho (v1).
    """
    def dfs(node, path):
        if node.v2_set.any() or node.v2_complementar_set.any() or node.v2_duplicado_set.any():
            caminho_formatado = "-".join(map(str, path))
            print(f"Caminho: {caminho_formatado}")
            if node.v2_set.any():
                print("  v2_set:", [i for i, bit in enumerate(node.v2_set) if bit])
            if node.v2_complementar_set.any():
                print("  v2_complementar_set:", [i for i, bit in enumerate(node.v2_complementar_set) if bit])
            if node.v2_duplicado_set.any():
                print("  v2_duplicado_set:", [i for i, bit in enumerate(node.v2_duplicado_set) if bit])
            print("-" * 40)
        for child_value, child_node in enumerate(node.children):
            if child_node is not None:
                dfs(child_node, path + [child_value])
    dfs(trie.root, [])

def precomputar_gc_cpg(l):

    base_map = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    dicionario = {}

    for indice in range(4**l):
        # Converte o índice para a sequência correspondente
        sequencia = []
        temp = indice
        for _ in range(l):
            sequencia.append(base_map[temp % 4])
            temp //= 4
        sequencia = ''.join(reversed(sequencia))

        # Calcula o total de G e C e as ocorrências de CG
        total_gc = sum(1 for char in sequencia if char in 'GC')
        ocorrencias_cg = sum(1 for i in range(len(sequencia) - 1) if sequencia[i:i+2] == 'CG')

        # Armazena os valores no dicionário
        dicionario[indice] = (total_gc, ocorrencias_cg)

    return dicionario

def calcular_parametros(matches_n, l, k):
    base_map = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    resultados = []

    def converter_para_bases(indice, tamanho):
        bases = []
        for _ in range(tamanho):
            bases.append(base_map[indice % 4])
            indice //= 4
        return ''.join(reversed(bases))

    def calcular_gc_percent(sequencia):
        return (sum(1 for base in sequencia if base in 'GC') * 100) / len(sequencia)

    def calcular_cpg(sequencia):
        return sum(1 for i in range(len(sequencia) - 1) if sequencia[i:i + 2] == 'CG')

    def eh_palindromo(sequencia):
        return sequencia == sequencia[::-1]

    # Iterar sobre os caminhos e valores da Trie
    for caminho, node in matches_n.iterate():  # Supondo que `iterate` percorre os nós da Trie
        caminho_v1 = ''.join(base_map[num] for num in caminho)  # Converter caminho v1 para bases
        for v2 in node.v2_set:  # Supondo que `v2_set` contém os valores de v2
            v2_bases = converter_para_bases(v2, l)
            sequencia_completa = caminho_v1 + v2_bases

            # Calcular os parâmetros
            gc_percent = calcular_gc_percent(sequencia_completa)
            cpg = calcular_cpg(sequencia_completa)
            palindromo = eh_palindromo(sequencia_completa)

            # Adicionar os resultados à lista
            resultados.append({
                "sequencia": sequencia_completa,
                "gc_percent": gc_percent,
                "cpg": cpg,
                "palindromo": int(palindromo)
            })

    return resultados

def importar_trie_em_txt_bitarray(arquivo_entrada, l, m):
    """
    Lê um arquivo no formato exportado por escrever_trie_em_txt_bitarray e reconstrói uma TrieBit.
    :param arquivo_entrada: Caminho do arquivo .txt.
    :param l: Comprimento de v1.
    :param m: Comprimento de v2 (tamanho do bitarray).
    :return: Instância de TrieBit reconstruída.
    """
    trie = TrieBit(m)
    with open(arquivo_entrada, 'r') as f:
        v1_lexico = None
        sets = {'+': [], '-': [], '=': []}
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Se já havia um v1 anterior, insere os valores na Trie
                if v1_lexico is not None:
                    v1 = []
                    temp = v1_lexico
                    for _ in range(l):
                        v1.append(temp // (4 ** (l - len(v1) - 1)))
                        temp = temp % (4 ** (l - len(v1)))
                    node = trie.root
                    for value in v1:
                        if node.children[value] is None:
                            node.children[value] = TrieNodeBit(m)
                        node = node.children[value]
                    # Preenche os bitarrays
                    for v2 in sets['+']:
                        node.v2_set[v2] = 1
                    for v2 in sets['-']:
                        node.v2_complementar_set[v2] = 1
                    for v2 in sets['=']:
                        node.v2_duplicado_set[v2] = 1
                # Novo v1
                v1_lexico = int(line[1:])
                sets = {'+': [], '-': [], '=': []}
            elif line in ['+', '-', '=']:
                tipo = line
                valores = next(f).strip()
                if valores:
                    sets[tipo] = [int(x) for x in valores.split(',') if x]
        # Último bloco
        if v1_lexico is not None:
            v1 = []
            temp = v1_lexico
            for _ in range(l):
                v1.append(temp // (4 ** (l - len(v1) - 1)))
                temp = temp % (4 ** (l - len(v1)))
            node = trie.root
            for value in v1:
                if node.children[value] is None:
                    node.children[value] = TrieNodeBit(m)
                node = node.children[value]
            for v2 in sets['+']:
                node.v2_set[v2] = 1
            for v2 in sets['-']:
                node.v2_complementar_set[v2] = 1
            for v2 in sets['=']:
                node.v2_duplicado_set[v2] = 1
    return trie

def importar_trie_em_csv_bitarray(arquivo_csv, l, m):
    """
    Reconstrói uma TrieBit a partir de um arquivo CSV gerado por escrever_trie_em_csv_bitarray.
    :param arquivo_csv: Caminho do arquivo CSV.
    :param l: Comprimento de v1.
    :param m: Tamanho do bitarray (4**l).
    :return: Instância de TrieBit reconstruída.
    """
    trie = TrieBit(m)
    import csv

    with open(arquivo_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            v1 = int(row['v1'])
            v2 = int(row['v2'])
            identificador = row['Identificador']
            # Converte v1 (índice lexicográfico) para caminho de bases
            path = []
            temp = v1
            for _ in range(l):
                path.append(temp // (4 ** (l - len(path) - 1)))
                temp = temp % (4 ** (l - len(path)))
            node = trie.root
            for value in path:
                if node.children[value] is None:
                    node.children[value] = TrieNodeBit(m)
                node = node.children[value]
            # Marca o bit correto conforme o identificador
            if identificador == '+':
                node.v2_set[v2] = 1
            elif identificador == '-':
                node.v2_complementar_set[v2] = 1
            elif identificador == '=':
                node.v2_duplicado_set[v2] = 1
    return trie

def contar_nulomeros_trie_bit(trie):
    """
    Conta o total de nulômeros em uma TrieBit.
    Um nulômero é cada bit setado em v2_set, v2_complementar_set ou v2_duplicado_set de qualquer nó.
    """
    def dfs(node):
        total = node.v2_set.count() + node.v2_complementar_set.count() + node.v2_duplicado_set.count()
        for child in node.children:
            if child is not None:
                total += dfs(child)
        return total
    return dfs(trie.root)

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

# ---------------------------------------------------------------------------------------------------------------#
# Abordagem via matches
# ---------------------------------------------------------------------------------------------------------------#

def checar_palíndromos(nulomeros,k, organismo, caminho_out): # busca entre os nulômeros os que são palíndromos

    palindromos = []

    with open(nulomeros) as f:

        while True:
        
            atual = f.read(k)

            if not atual:
                break
            prox = f.read(1)
            if eh_palindromo(atual):
                if prox == '+':
                    palindromos.append(atual + '+')
                    f.read(1)
                elif prox == '-':
                    palindromos.append(atual + '-')
                    f.read(1)
                else:
                    palindromos.append(atual)

    with open(caminho_out + 'palindromos' + organismo + '.txt', 'w') as f_out:
        f_out.write("Sequências palíndromas:\n")
        for palindromo in palindromos:
            f_out.write(palindromo + '\n')

    return len(palindromos)

def gerar_dicionario_combinado(matches, matchescomp, n1_out, k, organismo): # separa matches em únicos por fita e comuns
    """
    Gera o dicionário combinado de nulômeros e escreve os resultados em um arquivo.
    """

    l = int(k / 2)
    kmat = int(pow(4, l))
    kmer_total = gerar_kmers(l)

    # Calcula as diferenças e a interseção
    exclusivas_original = set(matches.keys()) - set(matchescomp.keys())
    exclusivas_complementar = set(matchescomp.keys()) - set(matches.keys())
    comuns = set(matches.keys()) & set(matchescomp.keys())
    # Cria o novo dicionário combinado
    dicionario_combinado = {}
    # Adiciona as chaves exclusivas do original com '+'
    for chave in exclusivas_original:
        dicionario_combinado[chave] = "-"
    # Adiciona as chaves exclusivas do complementar com '-'
    for chave in exclusivas_complementar:
        dicionario_combinado[chave] = "+"
    # Adiciona as chaves comuns com '='
    for chave in comuns:
        dicionario_combinado[chave] = ""
    # Escreve os resultados no arquivo
    with open(n1_out, 'w') as f:
        cont_null = 0
        null_original = 0
        null_complementar = 0
        for i in range(kmat):
            for j in range(kmat):
                chave = (i, j)
                if chave in exclusivas_original:
                    # Se a chave estiver em matches ou matchescomp, escreve o valor correspondente
                    f.write(f"{kmer_total[i]}{kmer_total[j]}{dicionario_combinado[chave]}\n")
                    null_complementar += 1
                elif chave in exclusivas_complementar:
                    # Se a chave estiver em matchescomp, escreve o valor correspondente
                    f.write(f"{kmer_total[i]}{kmer_total[j]}{dicionario_combinado[chave]}\n")
                    null_original += 1
                elif chave not in comuns:
                    # Apenas pares ausentes em ambos os dicionários são considerados nulômeros
                    f.write(f"{kmer_total[i]}{kmer_total[j]}\n")
                    cont_null += 1
    print(f"Gerando arquivo em: {n1_out}")
    print(f"Total de nulômeros encontrados em ambas fitas: {cont_null}.\nNulômeros na fita original: {null_original}.\nNulômeros na fita complementar: {null_complementar}.")
    return cont_null, null_original, null_complementar

def transformar_txt_para_csv(arquivo_txt, arquivo_csv):
    """
    Transforma o arquivo .txt gerado por gerar_dicionario_combinado em um arquivo .csv.
    O CSV conterá duas colunas: Nulômero e Identificador (+, -, ou vazio).

    :param arquivo_txt: Caminho para o arquivo .txt de entrada.
    :param arquivo_csv: Caminho para o arquivo .csv de saída.
    """
    with open(arquivo_txt, 'r') as txt_file, open(arquivo_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        # Escreve o cabeçalho do CSV
        writer.writerow(["Nulômero", "Identificador"])

        for line in txt_file:
            # Remove espaços e quebras de linha
            line = line.strip()
            if not line:
                continue  # Ignora linhas vazias

            # Separa o nulômero do identificador
            nulomero = line[:-1] if line[-1] in ['+', '-'] else line
            identificador = line[-1] if line[-1] in ['+', '-'] else ""

            # Escreve no CSV
            writer.writerow([nulomero, identificador])

def porc_gc(entrada, k, saida): # calcula a porcentagem de GC de cada nulômero
    """
    Calcula o conteúdo de GC nos primeiros k caracteres de cada linha do arquivo de entrada
    e escreve os resultados no arquivo de saída.
    """
    with open(entrada, 'r') as f, open(saida, 'w') as f_out:
        buffer = []  # Buffer para acumular linhas antes de escrever
        for line in f:
            # Pega apenas os primeiros k caracteres da linha
            vector = line[:k]
            if len(vector) < k:
                continue  # Ignora linhas menores que k caracteres

            # Calcula o conteúdo de GC
            cont = sum(1 for char in vector if char in 'CG')
            porc_gc = (cont * 100) / k

            # Adiciona o resultado ao buffer
            buffer.append(f"{vector},{porc_gc:.4f}\n")

            # Escreve o buffer no arquivo de saída em blocos
            if len(buffer) >= 10000:  # Ajuste o tamanho do buffer conforme necessário
                f_out.writelines(buffer)
                buffer = []

        # Escreve o restante do buffer
        if buffer:
            f_out.writelines(buffer)
                
def null_cpg(entrada, k, saida): # calcula o número de dinucleotídeos CG em cada nulômero
    """
    Conta o número de dinucleotídeos CG nos primeiros k caracteres de cada linha do arquivo de entrada
    e escreve no arquivo de saída apenas as sequências que possuem pelo menos um CG,
    junto com o número total de CGs e o identificador (+, -, ou vazio).
    """
    total_cpg = 0  # Contador total de CGs

    with open(entrada, 'r') as f, open(saida, 'w') as f_out:
        for line in f:
            # Pega os primeiros k caracteres da linha
            vector = line[:k]
            if len(vector) < k:
                continue  # Ignora linhas menores que k caracteres

            # Verifica o identificador (+, -, ou vazio)
            identificador = line[k:k+1] if len(line) > k else ''

            # Conta os dinucleotídeos CG
            num_cpg = sum(1 for i in range(len(vector) - 1) if vector[i] == 'C' and vector[i+1] == 'G')

            # Escreve no arquivo de saída apenas se houver pelo menos um CG
            if num_cpg > 0:
                f_out.write(f"{vector}{identificador}""{num_cpg}\n")
                total_cpg += 1

    return total_cpg

def matches(genoma, l, k): # marca quais sequencias existem no genoma e em sua fita complementar
    """
    Lê o genoma e calcula os pares de nulômeros (i1, i2) e seus complementares (i1c, i2c).
    Retorna dois dicionários: matches (original) e matchescomp (complementar).
    """
    current_time = dt.datetime.now()
    inicio = current_time
    print("Iniciando leitura do arquivo.")
    print(inicio)
    
    cromo_ciclo = 0  # Contador de ciclo para cromossomos
    contar_ciclo = 0  # Contador para total de kmers analisados
    matches = {}
    matchescomp = {}
    cromos = {}

    with open(genoma, 'r') as f:
        cromo = f.readline()
        print("Cromossomo:", cromo)
 
        vector = f.read(k)  # Lê os primeiros k caracteres do genoma

        while len(vector) == k:
            if vector[k-1] == '\n':
                vector = vector[:k-1]
                next_char = f.read(1)
                if next_char == '>':
                    cromo = f.readline()
                    cromo_ciclo += 1
                    print("Ciclos do cromossomo:", contar_ciclo)
                    print("Cromossomo:", cromo)
                    vector = f.read(k)
                    cromos[cromo_ciclo] = contar_ciclo
                else:
                    vector = vector + next_char
                    if len(vector) != k:
                        continue

            v1 = vector[0:l]
            v2 = vector[l:k]

            # Calcula o índice de cada kmer de v1 e v2
            i1 = int(calc_ind_lexicografico(v1, l))
            i2 = int(calc_ind_lexicografico(v2, l))

            # Calcula os índices complementares
            i1c = int(calc_ind_lexicografico_complementar(i1, l))
            i2c = int(calc_ind_lexicografico_complementar(i2, l))

            # Adiciona os pares (i1, i2) ao dicionário matches
            if (i1, i2) not in matches:
                matches[(i1, i2)] = []
            matches[(i1, i2)].append(contar_ciclo)

            # Adiciona os pares complementares (i1c, i2c) ao dicionário matchescomp
            if (i1c, i2c) not in matchescomp:
                matchescomp[(i1c, i2c)] = []
            matchescomp[(i1c, i2c)].append(contar_ciclo)

            # Atualiza o vetor para o próximo k-mer
            vector = vector + f.read(1)
            vector = vector[1:k+1]
            contar_ciclo += 1
    return matches, matchescomp, contar_ciclo

# ----------------------------------------------------------------------------------------------------------------#
# Funções em comum
# ----------------------------------------------------------------------------------------------------------------#

def calc_ind_lexicografico(bases,l): # calcula o ind tendo como base a ordem lexicográfica
    base_valores = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    indice = 0
    for i, base in enumerate(bases):
        peso = 4 ** (l - i - 1)  # Peso da posição
        indice += base_valores[base] * peso

    return indice

def indice_para_seq(indice, l): # converte o ind lexicográfico em bases
    """
    Converte um índice lexicográfico em uma sequência de bases.

    Args:
        indice (int): Índice lexicográfico.
        l (int): Comprimento do k-mer.

    Returns:
        str: Sequência de bases correspondente ao índice.
    """
    base_map = {0: 'A', 1: 'T', 2: 'C', 3: 'G'}
    sequencia = []
    for _ in range(l):
        sequencia.append(base_map[indice % 4])
        indice //= 4
    return ''.join(reversed(sequencia))

def calc_ind_lexicografico_complementar(i1, l): # calcula o ind tendo como base os matches da fita original
    """
    Calcula o índice lexicográfico complementar ao índice i1.

    Args:
        i1 (int): Índice lexicográfico da sequência original.
        l (int): Comprimento do k-mer.

    Returns:
        int: Índice lexicográfico complementar.
    """
    # Mapeamento de valores complementares
    complementar = {0: 1, 1: 0, 2: 3, 3: 2}

    # Desmontar o índice em base 4 com preenchimento de zeros
    bases = []
    for _ in range(l):
        bases.append(i1 % 4)  # Extrai o valor da base (mod 4)
        i1 //= 4              # Move para a próxima posição
    bases = bases[::-1]  # Inverte a ordem para que fique no formato correto

    # Aplicar a complementaridade
    bases_complementares = [complementar[base] for base in bases]

    # Recalcular o índice lexicográfico para a sequência complementar
    indice_complementar = 0
    for i, base in enumerate(bases_complementares):
        indice_complementar += base * (4 ** (l - i - 1))

    return indice_complementar

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

def create_matrix(rows, cols): # deprecated, cria matriz de zeros para comparações
    matrix = []
    for _ in range(rows):
        row = [0] * cols
        matrix.append(row)
    return matrix

def checar_ult_base(base,x,l): # deprecated

    base_map1 = {

        'A' : x,
        'T' : x * 2,
        'C' : x * 3,
        'G' : x * 4

    }

    base_map2 = {

        'A' : 1,
        'T' : x + 1,
        'C' : 2 * x + 1,
        'G' : 3 * x + 1
    }

    if base[l-1] == 'G':

        return base_map1.get(base[0], 0)

    elif base[l-1] == 'A':

        return base_map2.get(base[0],0)
    
def calcular_valores_x(l): # deprecated, calcula os valores de x para calculo de indice pg
    
    xdic = {}
    kaux = 0
    while l >= 0:

        x = int(1*pow(4,l-1))
        xdic[kaux] = x
        l -= 1
        kaux += 1

    return xdic

def remover_quebras_linha(caminho_arquivo): # deprecated, remove quebras de linha do arquivo com nulômeros
    with open(caminho_arquivo, 'r') as file:
        conteudo = file.read().replace('\n', '')  # Remove todas as quebras de linha
    with open(caminho_arquivo, 'w') as file:
        file.write(conteudo)
    return caminho_arquivo

def fita_complementar(arquivo_entrada, arquivo_saida): # deprecated, cria fita complementar

# Abrir o arquivo de entrada e processar
    with open(arquivo_entrada, "r") as f_in, open(arquivo_saida, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):  # Cabeçalho do FASTA
                f_out.write(line)  # Escreve o cabeçalho no arquivo de saída
            else:
                dna_seq = Seq(line.strip())
                comp_seq = dna_seq.complement()
                f_out.write(str(comp_seq) + "\n")  # Escreve a sequência complementar no arquivo de saída

def ler_nulomeros_array(genoma,l,k,xdic,organismo,n1_out): # maior demanda de tempo, substituido por matches

    current_time = dt.datetime.now()
    inicio = current_time
    print("Iniciando leitura do arquivo.")
    print(inicio)
    
    cromo_ciclo = 0 # Contador de ciclo para cromossomos
    contar_ciclo = 0 # Contador para total de kmers analisados
    kmat = int(pow(4,l))
    kmer_total = gerar_kmers(l)
    matches = {}
    cromos = {}

    with open(genoma, 'r') as f:
        cromo = f.readline()
        print("Cromossomo:", cromo)
 
        vector = f.read(k) # Lê os primeiros k caracteres do genoma

        while len(vector) == k:
    
            if vector[k-1] == '\n':
    
                vector = vector[:k-1]
    
                next_char = f.read(1)
    
                if next_char == '>':

                    cromo = f.readline()

                    cromo_ciclo += 1

                    print("Ciclos do cromossomo:", contar_ciclo)

                    print("Cromossomo:", cromo)

                    vector = f.read(k)

                    cromos[cromo_ciclo] = contar_ciclo

                else:

                    vector = vector + next_char

                    if len(vector) != k:

                        continue

            v1 = vector[0:l]
            v2 = vector[l:k]
            
            # Calcula o índice de cada kmer de v1 segundo PG

            i1 = int(calc_ind_lexicografico(v1,l))

            i2 = int(calc_ind_lexicografico(v2,l))

            if (i1,i2) not in matches:
                matches[(i1,i2)] = []
                
            matches[(i1,i2)].append(contar_ciclo)

            vector = vector + f.read(1)
    
            vector = vector[1:k+1]

            contar_ciclo += 1


    print("Iniciando anotação")
    final = dt.datetime.now() - inicio

    with open(n1_out + 'nulomeros_' + organismo + '_' + str(k) + '.txt', 'w') as f:
        f.write(f"Organismo: {organismo}\n")
        f.write(f"Tempo total de análise: {final}\n")
        cont_null = 0
        for i in range(kmat):
            for j in range(kmat):
                if (i, j) not in matches:
                    f.write(kmer_total[i] + kmer_total[j] + '\n')
                    cont_null += 1
        

    with open(n1_out + 'pares_ocorridos_' + organismo + '_' + str(k) + '.txt', 'w') as f:
        f.write(f"Organismo: {organismo}\n")
        f.write(f"Tempo total de análise: {final}\n")
        for (i1, i2), ciclos in matches.items():
            # Escreve o par de nulômeros e os ciclos em que ocorreram
            f.write(f"Par {kmer_total[i1]}{kmer_total[i2]} (match: ({i1}, {i2})) ocorreu nos ciclos: {ciclos}\n")
    return cont_null

def calcular_indice_pg(base,l,xdic): # substituido por calc_ind_lexicografico
        
        index = np.zeros(l)
        for ind in range(l):
            
                x = xdic[ind]
                if base[ind] == 'A':
                    index[ind] = 1
                elif base[ind] == 'T':
                    index[ind] = x + 1
                elif base[ind] == 'C': 
                    index[ind] = 2 * x + 1
                elif base[ind] == 'G':
                    index[ind] = 3 * x + 1

        for ind in range(l - 1, 0, -1):
                index[ind] -= 1
                index[ind - 1] += index[ind]
        
        return index[0]

def calc_ind_lexicografico_trie(seq1,seq2,trie,l): # calcula o ind tendo como base a ordem lexicográfica
    base_valores = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    indice = 0
    for i, seq in enumerate(seq1):
        peso = 4 ** (l - i - 1)  # Peso da posição
        indice = base_valores[seq] * peso
        trie.insert(indice)  # Adiciona o índice à Trie

    for i, seq2 in enumerate(seq2):
        peso = 4 ** (l - i - 1)  # Peso da posição
        indice += base_valores[seq] * peso

    return indice

def remover_null_duplos(original, complementar, organismo, k, caminho_out): # deprecated, substituido por matches
    """
    Remove nulômeros duplicados entre dois arquivos e escreve as sequências exclusivas
    e comuns em um arquivo de saída, com marcações '+' ou '-'.
    """
    # Lê os arquivos e filtra sequências de tamanho k
    with open(original, 'r') as f1, open(complementar, 'r') as f2:
        seq1 = {line.strip() for line in f1 if len(line.strip()) == k}
        seq2 = {line.strip() for line in f2 if len(line.strip()) == k}

    # Calcula as diferenças e interseção
    apenas_seq1 = seq1 - seq2
    apenas_seq2 = seq2 - seq1
    intersecao = seq1 & seq2
    total = len(apenas_seq1) + len(apenas_seq2) + len(intersecao)  

    # Escreve todas as sequências no arquivo de saída
    caminho_saida = f"{caminho_out}sem_duplos_{organismo}.txt"
    with open(caminho_saida, 'w') as f_out:
        f_out.write(f"Total de nulômeros: {total}\n")
        # Escreve as sequências exclusivas de seq1 com '+'
        f_out.writelines(f"{seq}+\n" for seq in sorted(apenas_seq1))
        # Escreve as sequências exclusivas de seq2 com '-'
        f_out.writelines(f"{seq}-\n" for seq in sorted(apenas_seq2))
        # Escreve as sequências comuns sem marcação
        f_out.writelines(f"{seq}\n" for seq in sorted(intersecao))

    return caminho_saida

def gerar_trie(v1,v2,l): # gera trie para os kmers, deprecated
    base_values = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    v1 = [base_values[base] for base in v1]
    v2 = calc_ind_lexicografico(v2, l)  # `v2` é um número inteiro
    v1 = v1 + [v2]
    # Insere v1 na Trie
    trie.insert(v1)
    trie.print_trie()
    return trie

def checar_inicio_fim(sequencia):
    """
    Verifica se a sequência começa com G e se ela termina com C.
    :param sequencia: Sequência de nucleotídeos.
    :return: Valor binário correspondente a cada caso, seqg e seqc.
    """
    seqg = 1 if sequencia.startswith('G') else 0
    seqc = 1 if sequencia.endswith('C') else 0

    return seqg, seqc