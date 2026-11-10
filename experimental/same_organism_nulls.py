import re
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

DIR = Path("workflow/results_scerevisiae_datasets_sequence/")
K = 11
MIN_COUNT, MAX_COUNT = 2, 20        # nullomer presente em 2..20 genomas
DIST_CUT = 0.7                      # corte do dendrograma (1 - Jaccard); ajuste olhando o dendrograma
METADATA = "metadados_leveduras_gca.csv"
GCA_RE = re.compile(r"GCA_\d{9}\.\d+")

LUT = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
POW = 4 ** np.arange(K - 1, -1, -1, dtype=np.int64)


def ler_codigos(path):
    """Lê um arquivo de nullomers e devolve um array de inteiros únicos."""
    raw = Path(path).read_bytes().upper().replace(b"\r", b"").replace(b"\n", b"")
    arr = LUT[np.frombuffer(raw, dtype=np.uint8)].reshape(-1, K)
    if (arr == 255).any():
        raise ValueError(f"Caracteres fora de ACGT em {path}")
    return np.unique(arr.astype(np.int64) @ POW)


def decodificar(code):
    return "".join("ACGT"[(code // 4 ** (K - 1 - i)) % 4] for i in range(K))


files = sorted(DIR.glob(f"*_{K}_result"))
nomes = [GCA_RE.search(f.name).group(0) if GCA_RE.search(f.name) else f.name for f in files]

# ---- passo 1: contagem global (um bincount de 4^12 posições) ----
contagem = np.zeros(4 ** K, dtype=np.int32)
for f in files:
    contagem[ler_codigos(f)] += 1
    print("contado:", f.name)

# histograma (substitui o plt.hist(d))
presentes = contagem[contagem > 0]
plt.hist(presentes, bins=np.arange(1, presentes.max() + 2) - 0.5)
plt.axvline(MAX_COUNT + 0.5, color="r", ls="--")
plt.xlabel("nº de genomas em que o nullomer ocorre"); plt.ylabel("nº de nullomers")
plt.yscale("log"); plt.savefig("hist_nullomers.png", dpi=150); plt.close()

# ---- passo 2: matriz esparsa organismo x nullomer raro ----
raros = np.flatnonzero((contagem >= MIN_COUNT) & (contagem <= MAX_COUNT))
print(f"{len(raros)} nullomers raros")
col_de = {c: i for i, c in enumerate(raros)}   # código -> coluna
lut_col = np.full(4 ** K, -1, dtype=np.int32)
lut_col[raros] = np.arange(len(raros))

linhas, colunas = [], []
for r, f in enumerate(files):
    cols = lut_col[ler_codigos(f)]
    cols = cols[cols >= 0]
    linhas.append(np.full(len(cols), r)); colunas.append(cols)
M = sparse.csr_matrix(
    (np.ones(sum(map(len, colunas)), dtype=np.int32),
     (np.concatenate(linhas), np.concatenate(colunas))),
    shape=(len(files), len(raros)))

# ---- passo 3: Jaccard + clustering hierárquico ----
inter = (M @ M.T).toarray().astype(float)
tam = np.asarray(M.sum(axis=1)).ravel()
uniao = tam[:, None] + tam[None, :] - inter
jacc = np.divide(inter, uniao, out=np.zeros_like(inter), where=uniao > 0)
dist = 1 - jacc
np.fill_diagonal(dist, 0)

Z = linkage(squareform(dist, checks=False), method="average")
grupos = fcluster(Z, t=DIST_CUT, criterion="distance")

df = pd.DataFrame({"GCA": nomes, "Cluster": grupos, "N_nullomers_raros": tam})
df["Tamanho_cluster"] = df.groupby("Cluster")["GCA"].transform("size")
df.to_csv("grupos_nullomers_raros.csv", index=False)

# ---- passo 4: cruzar com os metadados ----
meta = pd.read_csv(METADATA, dtype=str)
df["GCA"] = df["GCA"].str.replace(".fasta_\d\d_result","", regex = True)
meta["GCA"] = meta["GCA"].str.replace(".","_")
print(meta.columns.tolist())
print(df["GCA"].head().tolist())
print(meta["GCA"].head().tolist() if "GCA" in meta else "sem coluna GCA")
print("casam:", df["GCA"].isin(meta["GCA"]).sum(), "de", len(df))
cruz = df.merge(meta, on="GCA", how="left").sort_values(["Tamanho_cluster", "Cluster"], ascending=[False, True])
cruz.to_csv("grupos_com_metadados.csv", index=False, encoding="utf-8-sig")

# resumo por cluster: valor mais frequente de cada campo de origem
def top(s):
    s = s.dropna().replace("", np.nan).dropna()
    return "; ".join(f"{v} ({n})" for v, n in s.value_counts().head(3).items())

campos = ["Isolation_source", "Host", "Geo_loc", "BioProject", "Submitter"]
resumo = (cruz.groupby("Cluster")
          .agg(N=("GCA", "size"), **{c: (c, top) for c in campos})
          .sort_values("N", ascending=False))
resumo.to_csv("resumo_clusters.csv", encoding="utf-8-sig")
print(resumo.head(20))
