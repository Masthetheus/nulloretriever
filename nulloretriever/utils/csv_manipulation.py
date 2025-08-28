"""General utils aimed to csv analysis files manipulation and data retrieving
UNDER DEVELOPMENT:
    testing needed
    compatibility needed
    integration with main pipeline needed
"""

import csv
import pandas as pd

def add_genome_size(csv_nullomers, csv_organisms):
    """Add genome size information to a csv containing genomes accession codes
    Args:
        csv_nullomers: csv file containing organism nullomer analysis data, generated from the current pipeline
        csv_organisms: csv file containing organism global data, including genome NCBI accession codes
    Returns:
        output: adds genome size to each row on csv_nullomers according to it's relative organism
    """
    df_null = pd.read_csv(csv_nullomers, sep=",")
    df_orgs = pd.read_csv(csv_organisms, sep="\t")

    df_orgs['Organism Name Normalized'] = df_orgs['Organism Name'].apply(normalizar_nome)
    df_null['Organism Normalized'] = df_null['Organism'].apply(normalizar_nome)

    df_null = pd.merge(
        df_null,
        df_orgs[['Organism Name Normalizaded', 'Genome Size']],
        left_on='Organismo Normalized',
        right_on='Organism Name Normalized',
        how='left'
    )

    df_null = df_null.drop(columns=['Organism Name Normalized', 'Organismo Normalized'])
    df_null.to_csv(csv_nullomers, sep=",", index=False)
    if 'Genome Size' in df_null.columns:
        df_null = df_null.drop(columns=['Genome Size'])

def add_nullomer_relative_percentage(file):
    """Retrieves nullomer percentage relative to genome size
    Args:
        file(str): path to main nullomer csv file
    Returns:
        output: adds nullomer percentage column
    Observations:
        file(str): must have genome size column
    """
    df = pd.read_csv(file, sep=",")
    df['Nullomer percent'] = (df['Total nullomer'] / df['Genome Size']) * 100
    df['Nullomer percent'] = df['Nullomer percent'].round(2)
    df.to_csv(file, sep=",", index=False)