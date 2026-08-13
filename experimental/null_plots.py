"""
Script completo para análise de nullômeros
- Derivação de colunas
- Regressão para counter, trivial_count e maws
- Cálculo de resíduos e análise por ordem
- Geração de gráficos e tabelas
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import statsmodels.formula.api as smf
import os
import warnings
warnings.filterwarnings('ignore')

# ---- CONFIGURAÇÃO ----
K_FOCUS = 12
ORDER_COL = "order_id"
OUTPUT_DIR = "workflow/results"
GRAPH_DIR = f"{OUTPUT_DIR}/graphs"
TABLE_DIR = OUTPUT_DIR
os.makedirs(GRAPH_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

print("="*60)
print("ANÁLISE DE NULLÔMEROS - SCRIPT COMPLETO")
print("="*60)

# ---- 1. CARREGAR DADOS ----
print("\n[1] Carregando dados...")
df = pd.read_csv("workflow/shuffled_results/filtered_nullomer_statistics_summarized.csv")
print(f"Total de linhas no CSV: {len(df)}")

# Verifica colunas obrigatórias
required_cols = ["k", "counter", "trivial", "Genome length", "composition", 
                 "motifs_cpg_total", "motifs_palindromy_count", 
                 "motifs_cpg_nullomers_with_cpg", "motifs_palindromy_relative_fraction",
                 ORDER_COL]
for col in required_cols:
    if col not in df.columns:
        raise KeyError(f"Coluna '{col}' não encontrada. Colunas disponíveis: {df.columns.tolist()}")

# ---- 2. DERIVAÇÃO DE COLUNAS ----
print("\n[2] Derivando colunas...")
df["trivial_count"] = df["trivial"]
df["maws"] = df["counter"] - df["trivial"]
df["log_genome"] = np.log10(df["Genome length"])
df["observed_kmers"] = 4**df["k"] - df["counter"]

# CpG O/E (razão Observado/Esperado)
gc_frac = df["composition"] / 100
c_freq = gc_frac / 2
g_freq = gc_frac / 2
expected_cpg = c_freq * g_freq * df["Genome length"]
df["cpg_oe"] = np.where(expected_cpg > 0, df["motifs_cpg_total"] / expected_cpg, np.nan)

# Densidade de palindromia (por Mb)
df["pal_density"] = df["motifs_palindromy_count"] / (df["Genome length"] / 1e6)

print(f"Colunas derivadas adicionadas. Total de colunas: {len(df.columns)}")

# ---- 3. FILTRAGEM E LIMPEZA ----
print(f"\n[3] Filtrando k={K_FOCUS} e removendo NaNs...")
df_k12 = df[df["k"] == K_FOCUS].copy()

# Verifica NaNs nas colunas de interesse
cols_for_model = [
    "counter", "trivial_count", "maws", 
    "log_genome", "composition", "cpg_oe", "pal_density",
    "motifs_cpg_nullomers_with_cpg", "motifs_palindromy_relative_fraction",
    ORDER_COL
]

# Remove linhas com NaN em qualquer uma dessas colunas
df_k12_clean = df_k12[cols_for_model].dropna()
print(f"Linhas antes da limpeza: {len(df_k12)}")
print(f"Linhas após remover NaNs: {len(df_k12_clean)}")

# Se não houver dados suficientes, interrompe
if len(df_k12_clean) < 10:
    raise ValueError("Poucos dados após remoção de NaNs. Verifique seus dados.")

# ---- 4. PREPARAÇÃO PARA MODELAGEM ----
print("\n[4] Preparando modelos...")
df_model = df_k12_clean.copy()
print(f"Usando {len(df_model)} organismos (k={K_FOCUS})")

# Variáveis preditoras contínuas (excluindo ordem)
predictor_cols = [
    "log_genome", "composition", "cpg_oe", "pal_density",
    "motifs_cpg_nullomers_with_cpg", "motifs_palindromy_relative_fraction"
]
X = df_model[predictor_cols]

# ---- 5. MODELO PARA COUNTER ----
print("\n[5.1] Modelo para counter (número total de nullômeros)")
y_counter = df_model["counter"]
model_counter = LinearRegression().fit(X, y_counter)
r2_counter = model_counter.score(X, y_counter)
df_model["residual_counter"] = y_counter - model_counter.predict(X)

# ---- 6. MODELO PARA TRIVIAL_COUNT ----
print("\n[5.2] Modelo para trivial_count")
y_trivial = df_model["trivial_count"]
model_trivial = LinearRegression().fit(X, y_trivial)
r2_trivial = model_trivial.score(X, y_trivial)
df_model["residual_trivial"] = y_trivial - model_trivial.predict(X)

# ---- 7. MODELO PARA MAWS ----
print("\n[5.3] Modelo para maws (nullômeros não-triviais)")
y_maws = df_model["maws"]
model_maws = LinearRegression().fit(X, y_maws)
r2_maws = model_maws.score(X, y_maws)
df_model["residual_maws"] = y_maws - model_maws.predict(X)

# ---- 8. RESULTADOS DOS MODELOS ----
print("\n[6] Resultados dos modelos numéricos (sem ordem):")
print(f"  R² counter       : {r2_counter:.4f}")
print(f"  R² trivial_count : {r2_trivial:.4f}")
print(f"  R² maws          : {r2_maws:.4f}")

# ---- 9. ANÁLISE POR ORDEM (resíduos médios) ----
print("\n[7] Resíduos médios por ordem...")

# Agrupa por ordem e calcula a média dos resíduos para cada variável
order_residuals = df_model.groupby(ORDER_COL).agg({
    "residual_counter": "mean",
    "residual_trivial": "mean",
    "residual_maws": "mean",
    "counter": "mean",
    "trivial_count": "mean",
    "maws": "mean",
    "log_genome": "mean",
    "composition": "mean"
}).reset_index()

# Ordena por resíduo de maws (mais positivo = mais MAWs que o esperado)
order_residuals_sorted = order_residuals.sort_values("residual_maws", ascending=False)

# Salva tabela completa
order_residuals_sorted.to_csv(f"{TABLE_DIR}/order_residuals_k{K_FOCUS}.csv", index=False)
print(f"Tabela de resíduos por ordem salva em: {TABLE_DIR}/order_residuals_k{K_FOCUS}.csv")

# ---- 10. MODELO COM ORDEM (para comparar ganho) ----
print("\n[8] Modelo com ordem (efeito fixo)...")

# Adiciona ordem como variável categórica
df_model["order_cat"] = df_model[ORDER_COL].astype("category")

# Modelo sem ordem (já temos, mas vamos refazer com statsmodels para teste LR)
X_const = sm.add_constant(X)
model_no_order = sm.OLS(y_maws, X_const).fit()

# Modelo com ordem (usando formula API para facilitar)
model_with_order = smf.ols("maws ~ " + " + ".join(predictor_cols) + " + C(order_cat)", data=df_model).fit()

# R² de ambos
r2_no_order = model_no_order.rsquared
r2_with_order = model_with_order.rsquared

# Teste de razão de verossimilhança
lr_stat = 2 * (model_with_order.llf - model_no_order.llf)
df_lr = model_with_order.df_model - model_no_order.df_model
p_lr = 1 - stats.chi2.cdf(lr_stat, df_lr)

print(f"R² sem ordem: {r2_no_order:.4f}")
print(f"R² com ordem: {r2_with_order:.4f}")
print(f"Ganho:        {r2_with_order - r2_no_order:.4f}")
print(f"LR test:      LR = {lr_stat:.2f}, p = {p_lr:.4e}")

# ---- 11. GRÁFICOS ----
print("\n[9] Gerando gráficos...")

# 11.1 Top 20 ordens com maior resíduo de MAWs
fig, ax = plt.subplots(figsize=(12, 6))
top20 = order_residuals_sorted.head(20)
ax.barh(top20[ORDER_COL].astype(str), top20["residual_maws"], color="steelblue")
ax.axvline(0, color='red', linestyle='--', linewidth=1)
ax.set_xlabel("Resíduo médio de MAWs (observado - predito)")
ax.set_title(f"Top 20 ordens com MAIS MAWs que o esperado (k={K_FOCUS})")
ax.invert_yaxis()  # Para a maior barra ficar no topo
fig.tight_layout()
fig.savefig(f"{GRAPH_DIR}/top20_order_residuals_maws_k{K_FOCUS}.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 11.2 Comparação de R² dos modelos
fig, ax = plt.subplots(figsize=(8, 5))
vars_names = ["Counter", "Trivial", "MAWs"]
r2_values = [r2_counter, r2_trivial, r2_maws]
bars = ax.bar(vars_names, r2_values, color=["#2c7bb6", "#d7191c", "#fdae61"])
ax.axhline(r2_with_order, color='green', linestyle='--', label=f"Com ordem (MAWs) = {r2_with_order:.3f}")
ax.set_ylabel("R²")
ax.set_title(f"R² dos modelos numéricos (sem ordem, k={K_FOCUS})")
ax.legend()
for bar, val in zip(bars, r2_values):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.02, f"{val:.3f}", ha='center', va='bottom')
fig.tight_layout()
fig.savefig(f"{GRAPH_DIR}/r2_comparison_k{K_FOCUS}.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 11.3 Scatter: resíduo de MAWs vs tamanho do genoma (por ordem)
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(order_residuals["log_genome"], order_residuals["residual_maws"], alpha=0.6)
ax.axhline(0, color='red', linestyle='--', linewidth=0.8)
ax.set_xlabel("log10(Genome length) médio por ordem")
ax.set_ylabel("Resíduo médio de MAWs")
ax.set_title(f"Resíduo de MAWs vs tamanho médio do genoma por ordem (k={K_FOCUS})")
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(f"{GRAPH_DIR}/residual_vs_genome_order_k{K_FOCUS}.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# 11.4 Scatter: resíduo de MAWs vs GC% médio por ordem
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(order_residuals["composition"], order_residuals["residual_maws"], alpha=0.6)
ax.axhline(0, color='red', linestyle='--', linewidth=0.8)
ax.set_xlabel("GC% médio por ordem")
ax.set_ylabel("Resíduo médio de MAWs")
ax.set_title(f"Resíduo de MAWs vs GC% médio por ordem (k={K_FOCUS})")
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(f"{GRAPH_DIR}/residual_vs_gc_order_k{K_FOCUS}.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ---- 12. TABELA RESUMO FINAL ----
summary_df = pd.DataFrame({
    "Variável": ["counter", "trivial_count", "maws"],
    "R² (numérico)": [r2_counter, r2_trivial, r2_maws],
    "R² (com ordem)": [r2_with_order]*3,  # mesmo valor para todos, apenas para referência
    "Ganho da ordem": [r2_with_order - r2_counter, r2_with_order - r2_trivial, r2_with_order - r2_maws]
})
summary_df.to_csv(f"{TABLE_DIR}/model_summary_k{K_FOCUS}.csv", index=False)
print(f"\nTabela resumo salva em: {TABLE_DIR}/model_summary_k{K_FOCUS}.csv")

# ---- 13. FINALIZAÇÃO ----
print("\n" + "="*60)
print("ANÁLISE CONCLUÍDA COM SUCESSO")
print("="*60)
print(f"Gráficos salvos em: {GRAPH_DIR}")
print(f"Tabelas salvas em: {TABLE_DIR}")
print("\nArquivos gerados:")
print(f"  - {GRAPH_DIR}/top20_order_residuals_maws_k{K_FOCUS}.png")
print(f"  - {GRAPH_DIR}/r2_comparison_k{K_FOCUS}.png")
print(f"  - {GRAPH_DIR}/residual_vs_genome_order_k{K_FOCUS}.png")
print(f"  - {GRAPH_DIR}/residual_vs_gc_order_k{K_FOCUS}.png")
print(f"  - {TABLE_DIR}/order_residuals_k{K_FOCUS}.csv")
print(f"  - {TABLE_DIR}/model_summary_k{K_FOCUS}.csv")
