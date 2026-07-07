import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats

# ── configuração ────────────────────────────────────────────────────────────
GRAPH_OUTPUT = "workflow/results/graphs"
K_FOCUS = 12          # k principal pra análise
CONSERV_THRESHOLD = 0.80  # 80% dos organismos do grupo

PHYLUM_FILTER = {
    "Ascomycota":    [4890],   # adaptar ao nome exato no dataframe
    "Basidiomycota": [5204],
}
PALETTE = {
    4890:    "#378ADD",
    5204: "#1D9E75",
}

# ── carrega e filtra ─────────────────────────────────────────────────────────
df = pd.read_csv("treated_final_results.csv")

# filtra só asco e basidio — adapta o nome da coluna de filo se necessário
print(df["phylum_id"])
df_fungi = df[df["phylum_id"].isin(
    PHYLUM_FILTER["Ascomycota"] + PHYLUM_FILTER["Basidiomycota"]
)].copy()
print(df_fungi)

df_k = df_fungi[df_fungi["k"] == K_FOCUS].copy()

# ── helpers ──────────────────────────────────────────────────────────────────
def spearman_label(x, y):
    mask = x.notna() & y.notna()
    r, p = stats.spearmanr(x[mask], y[mask])
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    return f"ρ = {r:.2f}, p {sig}"

def save(fig, name):
    fig.savefig(f"{GRAPH_OUTPUT}/{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

# ── figura 1 — tamanho de genoma vs % trivial (scatter central) ──────────────
def fig_genome_vs_trivial():
    fig, ax = plt.subplots(figsize=(7, 5))
    for phylum, color in PALETTE.items():
        sub = df_k[df_k["phylum_id"].isin([phylum])].copy()
        print(f"Sub = {sub}")
        ax.scatter(
            sub["genome_length"], sub["trivial_porc"],
            color=color, alpha=0.6, edgecolors="none", s=40, label=phylum
        )
    ax.set_xscale("log")
    ax.set_xlabel("Genome length (bp, log scale)", fontsize=11)
    ax.set_ylabel(f"Trivial nullomers (%, k={K_FOCUS})", fontsize=11)
    ax.set_title(spearman_label(df_k["genome_length"], df_k["trivial_porc"]),
                 fontsize=10, color="gray")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    save(fig, f"fig1_genome_vs_trivial_k{K_FOCUS}")

# ── figura 2 — proporção trivial através de k (mediana por filo) ─────────────
def fig_trivial_by_k():
    fig, ax = plt.subplots(figsize=(7, 5))
    for phylum, color in PALETTE.items():
        sub = df_fungi[df_fungi["phylum_id"] == phylum]
        # linhas individuais — baixa opacidade
        for _, org in sub.groupby("organism"):
            ax.plot(org["k"], org["trivial_porc"],
                    color=color, alpha=0.08, linewidth=0.8)
        # mediana grossa
        med = sub.groupby("k")["trivial_porc"].median()
        ax.plot(med.index, med.values,
                color=color, linewidth=2.5, label=f"{phylum} (median)")
    ax.set_xlabel("k", fontsize=11)
    ax.set_ylabel("Trivial nullomers (%)", fontsize=11)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    save(fig, "fig2_trivial_by_k")

# ── figura 3 — resíduo após controle por tamanho ─────────────────────────────
def fig_residuals():
    sub = df_k.dropna(subset=["genome_length", "trivial_porc"]).copy()
    log_genome = np.log10(sub["genome_length"])
    slope, intercept, *_ = stats.linregress(log_genome, sub["trivial_porc"])
    sub["residual"] = sub["trivial_porc"] - (slope * log_genome + intercept)

    fig, ax = plt.subplots(figsize=(7, 5))
    for phylum, color in PALETTE.items():
        s = sub[sub["phylum_id"] == phylum]
        ax.scatter(s["genome_length"], s["residual"],
                   color=color, alpha=0.6, edgecolors="none", s=40, label=phylum)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
    ax.set_xscale("log")
    ax.set_xlabel("Genome length (bp, log scale)", fontsize=11)
    ax.set_ylabel(f"Residual trivial % (k={K_FOCUS})", fontsize=11)
    ax.set_title("Deviation from genome-size expectation", fontsize=10, color="gray")
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    save(fig, f"fig3_residuals_k{K_FOCUS}")

# ── figura 4 — composição GC% nullomers vs GC% genoma ───────────────────────
def fig_gc_comparison():
    if "composition" not in df_k.columns or "gc_perc" not in df_k.columns:
        print("colunas de GC não encontradas, pulando fig4")
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    for phylum, color in PALETTE.items():
        s = df_k[df_k["phylum_id"] == phylum]
        ax.scatter(s["gc_perc"], s["composition"],
                   color=color, alpha=0.6, edgecolors="none", s=40, label=phylum)
    lims = [
        min(df_k["gc_perc"].min(), df_k["composition"].min()),
        max(df_k["gc_perc"].max(), df_k["composition"].max()),
    ]
    ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="y = x")
    ax.set_xlabel("Genome GC%", fontsize=11)
    ax.set_ylabel("Nullomer GC%", fontsize=11)
    ax.legend(frameon=False)
    sns.despine(ax=ax)
    save(fig, f"fig4_gc_nullomers_k{K_FOCUS}")

def fig_boxplot_by_k():
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(
        data=df_fungi, x="k", y="trivial_porc", hue="phylum_id",
        palette=PALETTE, flierprops={"marker": ".", "alpha": 0.4},
        linewidth=0.8, ax=ax
    )
    ax.set_xlabel("k", fontsize=11)
    ax.set_ylabel("Trivial nullomers (%)", fontsize=11)
    ax.legend(title="", frameon=False)
    sns.despine(ax=ax)
    save(fig, "fig5_boxplot_trivial_by_k")

def stats_kruskal():
    groups = [
        df_k[df_k["phylum_id"] == p]["trivial_porc"].dropna()
        for p in PALETTE
    ]
    stat, p = stats.kruskal(*groups)
    print(f"\nKruskal-Wallis trivial_porc (k={K_FOCUS}): H={stat:.3f}, p={p:.4f}")
    for p1, p2 in [(4890, 5204)]:
        g1 = df_k[df_k["phylum_id"] == p1]["trivial_porc"].dropna()
        g2 = df_k[df_k["phylum_id"] == p2]["trivial_porc"].dropna()
        u, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")
        print(f"  Mann-Whitney {p1} vs {p2}: U={u:.0f}, p={pval:.4f}")
# ── tabela descritiva ────────────────────────────────────────────────────────
def descriptive_table():
    tbl = df_fungi.groupby(["phylum_id", "k"])["trivial_porc"].agg(
        n="count", median="median",
        q25=lambda x: x.quantile(0.25),
        q75=lambda x: x.quantile(0.75)
    ).round(2)
    print("\nDescritiva trivial_porc por filo e k:")
    print(tbl.to_string())
    tbl.to_csv(f"{GRAPH_OUTPUT}/descriptive_trivial.csv")

# ── executa tudo ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import os; os.makedirs(GRAPH_OUTPUT, exist_ok=True)
    fig_genome_vs_trivial()
    fig_trivial_by_k()
    #fig_residuals()
    fig_gc_comparison()
    fig_boxplot_by_k()
    stats_kruskal()
    descriptive_table()
    print("\nPronto.")
