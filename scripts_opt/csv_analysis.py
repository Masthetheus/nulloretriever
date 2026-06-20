"""Analysis of gathered data csv from the snakemake pipeline."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

columns_erase = {'k', 'organism', 'organism_name',
                 'tax_id', 'motifs_homopolymers', 'class_id'}
df = pd.read_csv('workflow/results/nullomer_statistics_summarized_8_to_15.csv')
# for column in columns:
#     sns.lineplot(data=df, x="k",
#                  y=column, hue=' .')
#     plt.show()

dfj = pd.read_json('workflow/data/organisms.json').transpose()
mapping = dfj.set_index('SpeciesName')['phylum_id']
df['phylum_id'] = df['organism_name'].map(mapping)
df_class = df.groupby(['phylum_id'])
df_k = df.groupby(['k'])
columns = list(df)
columns_set = set(columns)
columns = list(columns_set-columns_erase)
# for column in columns:
#     sns.histplot(data=df, x="k",
#                  y=column, hue='class_id')
#     plt.show()
soma = df_class.describe()
# for name, group in df_class:
#     for column in columns:
#         sns.boxplot(data=group, x='k', y=column, hue='class_id')
#         plt.show()
#         sns.violinplot(data=group, x='k', y=column, hue='class_id')
#         plt.show()
# count = (group[1]['counter'] != 0).sum()
# print(count)
print(columns)
for name, group in df_k:
    print(group.describe())
for name, group in df_k:
    for column in columns:
        plt.figure(figsize=(12, 6))
        current = sns.boxplot(data=group, x='k', y=column,
                              hue='phylum_id', palette='Paired')
        sns.move_legend(current, "upper right", bbox_to_anchor=(1.05, 1))
        plt.savefig(f'graphs/{name[0]}/{column}', dpi=200)
        plt.close()
