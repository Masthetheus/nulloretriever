"""Calculate statistical measures from nullomer csv data."""

import pandas as pd
from sklearn.metrics import r2_score

K_FOCUS = 12
ORDER_COL = "order_id"

file = "workflow/results/filtered_nullomer_statistics_summarized.csv"

df = pd.read_csv(file)
df = df[df["k"] == K_FOCUS].copy()
df["maws"] = df["counter"] - df["trivial"]
df["observed_kmers"] = (4**(df["k"])) - df["counter"]
#df["log_genome"] = np.log10(df["Genome length"])
gc_frac = df["composition"] / 100
df["pal_density"] = df["motifs_palindromy_count"] / (df["Genome length"] / 1e6)

print(df.columns)
y = df["observed_kmers"].dropna()
f = df["counter"].dropna()
from sklearn.linear_model import LinearRegression
X = df[["observed_kmers"]]
y = df["counter"]
model = LinearRegression().fit(X, y)
r2 = model.score(X, y)
print(r2)
