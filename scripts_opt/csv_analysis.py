"""Analysis of gathered data csv from the snakemake pipeline."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv('workflow/results/nullomer_statistics_summarized_9_to_14.csv')
organismos = df.groupby(['organism', 'k'])
columns = list(df)
for column in columns:
    sns.lineplot(data=df, x="k",
                 y=column, hue=' .')
    plt.show()
