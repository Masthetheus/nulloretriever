from funcoes_mestrado import *
import os
import csv
import pandas as pd


def final_analysis_testing(
        csv_saida,
        null_file,
        k,
        org,
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
    resultados = []
    l = int(k // 2)
    parametros = parametros_nulomeros(null_file, l)
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

config_file = "../../config/config.yaml"
k = input("Which k value was used?")
if not k:
    k = 10
org = input("Which organism to perform the final analysis?")
if not org:
    null_file = "../data/teste_test_run"
    org = "teste"
else:
    null_file = "../data/" + org + "_test_run"
csv_out = "../data/" + org + "_data.csv"
final_analysis_testing(csv_out,null_file,k,org,config_file)

