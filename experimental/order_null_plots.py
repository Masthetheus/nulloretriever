import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
from scipy import stats
import os

# ── configuração ─────────────────────────────────────────────────────────────
GRAPH_OUTPUT = "workflow/results/graphs/by_order"
K_FOCUS      = 11
MIN_N        = 3       # ordens com n < MIN_N descartadas
GROUP_COL    = "order_id"
PHYLUM_IDS   = [4890, 5204]   # int — ascomycota, basidiomycota

# nomes legíveis por order_id (NCBI taxid → nome da ordem)
ORDER_NAMES = {
    5125:    "Hypocreales",
    4892:    "Saccharomycetales",
    5042:    "Eurotiales",
    5234:    "Tremellales",
    162474:  "Malasseziales",
    2916678: "Metschnikowiales",
    92860:   "Pleosporales",
    3243772: "Trichomonascales",
    3243779: "Dipodascales",
    3243775: "Phaffomycetales",
    5338:    "Agaricales",
    5139:    "Sordariales",
    33183:   "Onygenales",
    2726947: "Myriangiales",
    1028384: "Glomerellales",
    5178:    "Helotiales",
    5114:    "Diaporthales",
    162475:  "Microstromatales",
    1851469: "Cystofilobasidiales",
    231213:  "Sporidiobolales",
    2926619: "Saccharomycopsidales",
    3402561: "Amphisphaeriales",
    5267:    "Ustilaginales",
    5303:    "Polyporales",
    37989:   "Xylariales",
    451869:  "Botryosphaeriales",
    1484953: "Phacidiales",
    5404:    "Exobasidiales",
    36064:   "Cantharellales",
    34346:   "Schizosaccharomycetales",
}

# ── carrega e filtra ──────────────────────────────────────────────────────────
df = pd.read_csv("workflow/results/filtered_nullomer_statistics_summarized.csv")

# phylum_id pode vir como int ou float dependendo do pandas
df["phylum_id"] = pd.to_numeric(df["phylum_id"], errors="coerce")
df["order_id"]  = pd.to_numeric(df["order_id"],  errors="coerce")

df_fungi = df[df["phylum_id"].isin(PHYLUM_IDS)].copy()
df_k     = df_fungi[df_fungi["k"] == K_FOCUS].copy()

# descarta ordens com n < MIN_N (no k foco)
valid_orders = (
    df_k[GROUP_COL].value_counts()
    [lambda s: s >= MIN_N].index
)
df_k_f   = df_k[df_k[GROUP_COL].isin(valid_orders)].copy()
df_fun_f = df_fungi[df_fungi[GROUP_COL].isin(valid_orders)].copy()

# label legível: usa ORDER_NAMES se disponível, senão "ord_<id>"
def label(oid):
    return ORDER_NAMES.get(int(oid), f"ord_{int(oid)}")

orders     = sorted(df_k_f[GROUP_COL].dropna().unique())
order_lbls = [label(o) for o in orders]
n_ord      = len(orders)

# paleta automática
cmap = mpl.colormaps['tab20']
PALETTE = {o: cmap(i) for i, o in enumerate(orders)}
PALETTE_LBL = {label(o): cmap(i) for i, o in enumerate(orders)}

# filo de cada ordem (pra colorir borda se quiser)
order_phylum = df_k_f.groupby(GROUP_COL)["phylum_id"].agg(
    lambda x: x.mode()[0]
).to_dict()

print(f"Ordens válidas (n≥{MIN_N}): {n_ord}")
for o in orders:
    n = (df_k_f[GROUP_COL] == o).sum()
    ph = "Asco" if order_phylum.get(o) == 4890 else "Basidio"
    print(f"  {label(o):30s} n={n:3d}  [{ph}]")

# ── helpers ───────────────────────────────────────────────────────────────────
def spearman_label(x, y):
    mask = x.notna() & y.notna()
    r, p = stats.spearmanr(x[mask], y[mask])
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    return f"ρ = {r:.2f}, p {sig}"

def save(fig, name):
    os.makedirs(GRAPH_OUTPUT, exist_ok=True)
    fig.savefig(f"{GRAPH_OUTPUT}/{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

def legend_outside(ax, n):
    ncol = max(1, n // 12)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left",
              frameon=False, fontsize=8, ncol=ncol)

# ── fig 1 — genome length vs trivial % ───────────────────────────────────────
def fig_genome_vs_trivial():
    fig, ax = plt.subplots(figsize=(9, 6))
    for o in orders:
        sub = df_k_f[df_k_f[GROUP_COL] == o]
        # borda: tracejada pra basidio, sólida pra asco
        ec = "#5204AA" if order_phylum.get(o) == 5204 else "none"
        ax.scatter(sub["Genome length"], sub["trivial_porc"],
                   color=PALETTE[o], edgecolors=ec, linewidths=0.5,
                   alpha=0.75, s=50, label=label(o))
    ax.set_xscale("log")
    ax.set_xlabel("Genome length (bp, log scale)", fontsize=11)
    ax.set_ylabel(f"Trivial nullomers (%, k={K_FOCUS})", fontsize=11)
    ax.set_title(spearman_label(df_k_f["Genome length"], df_k_f["trivial_porc"]),
                 fontsize=10, color="gray")
    legend_outside(ax, n_ord)
    sns.despine(ax=ax)
    save(fig, f"ord_fig1_genome_trivial_k{K_FOCUS}")

# ── fig 2 — trivial % por k, mediana por ordem ───────────────────────────────
def fig_trivial_by_k():
    fig, ax = plt.subplots(figsize=(9, 6))
    for o in orders:
        sub = df_fun_f[df_fun_f[GROUP_COL] == o]
        for _, org_df in sub.groupby("organism"):
            ax.plot(org_df["k"], org_df["trivial_porc"],
                    color=PALETTE[o], alpha=0.06, linewidth=0.7)
        med = sub.groupby("k")["trivial_porc"].median()
        ax.plot(med.index, med.values,
                color=PALETTE[o], linewidth=2.2, label=label(o))
    ax.set_xlabel("k", fontsize=11)
    ax.set_ylabel("Trivial nullomers (%)", fontsize=11)
    legend_outside(ax, n_ord)
    sns.despine(ax=ax)
    save(fig, "ord_fig2_trivial_by_k")

# ── fig 3 — resíduo após controle por tamanho ────────────────────────────────
def fig_residuals():
    sub = df_k_f.dropna(subset=["Genome length", "trivial_porc"]).copy()
    log_g = np.log10(sub["Genome length"])
    slope, intercept, *_ = stats.linregress(log_g, sub["trivial_porc"])
    sub["residual"] = sub["trivial_porc"] - (slope * log_g + intercept)

    fig, ax = plt.subplots(figsize=(9, 6))
    for o in orders:
        s = sub[sub[GROUP_COL] == o]
        ax.scatter(s["Genome length"], s["residual"],
                   color=PALETTE[o], alpha=0.75, edgecolors="none",
                   s=50, label=label(o))
    ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel("Genome length (bp, log scale)", fontsize=11)
    ax.set_ylabel(f"Residual trivial % (k={K_FOCUS})", fontsize=11)
    ax.set_title("Desvio do esperado por tamanho de genoma", fontsize=10, color="gray")
    legend_outside(ax, n_ord)
    sns.despine(ax=ax)
    save(fig, f"ord_fig3_residuals_k{K_FOCUS}")

# ── fig 4 — boxplot ordens ordenadas por mediana ─────────────────────────────
def fig_boxplot_orders():
    df_plot = df_k_f.copy()
    df_plot["order_lbl"] = df_plot[GROUP_COL].map(label)
    order_med = df_plot.groupby("order_lbl")["trivial_porc"].median().sort_values()
    ordered_lbls = list(order_med.index)

    fig, ax = plt.subplots(figsize=(max(10, n_ord), 5))
    sns.boxplot(data=df_plot, x="order_lbl", y="trivial_porc",
                order=ordered_lbls, palette=PALETTE_LBL,
                flierprops={"marker": ".", "alpha": 0.4},
                linewidth=0.8, ax=ax)
    ax.set_xlabel("Order", fontsize=11)
    ax.set_ylabel(f"Trivial nullomers (%, k={K_FOCUS})", fontsize=11)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    sns.despine(ax=ax)
    save(fig, f"ord_fig4_boxplot_k{K_FOCUS}")

# ── fig 5 — heatmap mediana por ordem × k ────────────────────────────────────
def fig_heatmap():
    df_h = df_fun_f.copy()
    df_h["order_lbl"] = df_h[GROUP_COL].map(label)
    pivot = df_h.groupby(["order_lbl", "k"])["trivial_porc"].median().unstack()
    # ordena linhas pela mediana em k foco
    pivot = pivot.loc[pivot[K_FOCUS].sort_values().index]

    fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 1.2),
                                    max(6, len(pivot) * 0.45)))
    sns.heatmap(pivot, cmap="YlOrRd_r", annot=True, fmt=".1f",
                linewidths=0.3, ax=ax,
                cbar_kws={"label": "Trivial %", "shrink": 0.6})
    ax.set_xlabel("k", fontsize=11)
    ax.set_ylabel("")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8)
    save(fig, "ord_fig5_heatmap_order_k")

# ── kruskal-wallis ────────────────────────────────────────────────────────────
def stats_kruskal():
    groups = [df_k_f[df_k_f[GROUP_COL] == o]["trivial_porc"].dropna()
              for o in orders]
    H, p = stats.kruskal(*groups)
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"\nKruskal-Wallis entre ordens (k={K_FOCUS}): H={H:.3f}, p={p:.2e} {sig}")

    # pairwise asco vs basidio por ordem
    asco_orders   = [o for o in orders if order_phylum.get(o) == 4890]
    basidio_orders = [o for o in orders if order_phylum.get(o) == 5204]
    g1 = df_k_f[df_k_f[GROUP_COL].isin(asco_orders)]["trivial_porc"].dropna()
    g2 = df_k_f[df_k_f[GROUP_COL].isin(basidio_orders)]["trivial_porc"].dropna()
    u, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
    print(f"Mann-Whitney Ascomycota vs Basidiomycota: U={u:.0f}, p={pval:.4f}")

# ── tabela descritiva ─────────────────────────────────────────────────────────
def descriptive_table():
    df_t = df_k_f.copy()
    df_t["order_lbl"] = df_t[GROUP_COL].map(label)
    df_t["phylum"]    = df_t["phylum_id"].map({4890: "Ascomycota", 5204: "Basidiomycota"})
    tbl = df_t.groupby(["phylum", "order_lbl"])["trivial_porc"].agg(
        n="count", median="median",
        q25=lambda x: x.quantile(0.25),
        q75=lambda x: x.quantile(0.75)
    ).round(2).sort_values("median")
    print(f"\nDescritiva por ordem (k={K_FOCUS}):")
    print(tbl.to_string())
    tbl.to_csv(f"{GRAPH_OUTPUT}/descriptive_by_order_k{K_FOCUS}.csv")

# ── executa ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    os.makedirs(GRAPH_OUTPUT, exist_ok=True)
    fig_genome_vs_trivial()
    fig_trivial_by_k()
    #fig_residuals()
    fig_boxplot_orders()
    fig_heatmap()
    stats_kruskal()
    descriptive_table()
    print("\nPronto.")
