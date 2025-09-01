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

def ajustar_nome_orgs(nome_org):
    nome_org = nome_org.replace(" ", "_")  # Substitui espaços por underscores
    return re.sub(r'[^a-zA-Z0-9_]', '_', nome_org)  # Remove caracteres especiais, mantendo apenas letras, números e underscores

#----------------------------------------------------------------------#
# TrieBit approach
#----------------------------------------------------------------------#

# Base definitions

def revcomp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))

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
