"""Analysis of gathered data csv from the snakemake pipeline."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv('workflow/results/nullomer_statistics_summarized_9_to_14.csv')
organismos = df.groupby(['organism', 'k'])
print(list(df))
sns.lineplot(data=df, x="k",
             y='motifs_cpg_global_mean', hue='organism')
plt.show()
