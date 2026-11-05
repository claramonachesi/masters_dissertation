# Estatística para metaciclogênese
# Global: qui-quadrado de independência; Comparações dirigidas: teste exato de Fisher (MCM vs controles); Correção múltipla: Holm
# Script criado com auxílio do ChatGPT
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

sns.set_context("paper")
sns.set_style("whitegrid", {"grid.color": "#D3D3D3"})

df = pd.read_csv("metaciclogenese_2026_UERJ.tsv", sep="\t")
df.columns = [c.strip() for c in df.columns]

df["Epimastigotas"] = pd.to_numeric(df["Epimastigotas"], errors="raise").astype(int)
df["Diferenciadas"] = pd.to_numeric(df["Diferenciadas"], errors="raise").astype(int)

# Total contado por tratamento
df["Total"] = df["Epimastigotas"] + df["Diferenciadas"]
if not (df["Total"].isin([0, 300])).all():
    print("\n[AVISO] Há tratamentos com total diferente de 300 (ou 0).\n")

# Remover tratamentos sem células contadas (Total = 0)
df_zero = df[df["Total"] == 0].copy()
if not df_zero.empty:
    print("\nTratamentos removidos por total = 0:")
    print(df_zero[["Tratamentos", "Epimastigotas", "Diferenciadas", "Total"]])
df = df[df["Total"] > 0].copy().reset_index(drop=True)

## Derived variables
# Proporção de diferenciadas: desfecho principal
df["prop_diferenciadas"] = df["Diferenciadas"] / df["Total"]

# Proporção de epimastigotas: complementar
df["prop_epimastigotas"] = df["Epimastigotas"] / df["Total"]

# Probabilidade de permanecer como epimastigota versus diferenciar
df["odds_E_vs_D"] = np.where(df["Diferenciadas"] == 0, np.nan, df["Epimastigotas"] / df["Diferenciadas"])

print(df[["Tratamentos", "Epimastigotas", "Diferenciadas", "prop_epimastigotas", "prop_diferenciadas", "odds_E_vs_D"]])

## Global chi-square
# Tabela global: linhas = tratamentos; colunas = E e D
contingency_global = df[["Epimastigotas", "Diferenciadas"]].to_numpy()
chi2, p_chi2, dof, expected = chi2_contingency(contingency_global)
chi2_results = pd.DataFrame({"chi2": [chi2], "p_value": [p_chi2], "dof": [dof]})
chi2_results.to_csv("META_qui_quadrado_global.tsv", sep="\t", index=False)

print(f"\nQui-quadrado global: {chi2_results}")

# Salvar frequências esperadas
expected_df = pd.DataFrame(expected, columns=["Expected_Epimastigotas", "Expected_Diferenciadas"])
expected_df.insert(0, "Tratamentos", df["Tratamentos"])
expected_df.to_csv("META_frequencias_esperadas_qui_quadrado.tsv", sep="\t", index=False)

## Pairwise Fisher: cada MCM vs cada controle
controls = ["Controle", "Controle+DMSO", "LPC 10 uM"]
control_symbol_map = {"Controle": "*", "Controle+DMSO": "**", "LPC 10 uM": "***"}

# Selecionar MCMs
mcm_mask = df["Tratamentos"].str.upper().str.startswith("MCM")
mcm_list = df.loc[mcm_mask, "Tratamentos"].tolist()

results = []
for mcm in mcm_list:
    row_mcm = df[df["Tratamentos"] == mcm].iloc[0]
    if row_mcm["Total"] == 0:
        continue
    for ctrl in controls:
        row_ctrl = df[df["Tratamentos"] == ctrl]
        if row_ctrl.empty:
            print(f"[AVISO] Controle não encontrado: {ctrl}")
            continue
        row_ctrl = row_ctrl.iloc[0]
        if row_ctrl["Total"] == 0:
            continue

        # Tabela 2x2 – linhas: MCM, controle; colunas: Epimastigotas, Diferenciadas
        table_2x2 = np.array([
            [row_mcm["Epimastigotas"], row_mcm["Diferenciadas"]],
            [row_ctrl["Epimastigotas"], row_ctrl["Diferenciadas"]]])

        # Fisher
        odds_ratio, p_value = fisher_exact(table_2x2, alternative="two-sided")

        # Diferença de proporção de diferenciadas
        prop_diff_mcm = row_mcm["Diferenciadas"] / row_mcm["Total"]
        prop_diff_ctrl = row_ctrl["Diferenciadas"] / row_ctrl["Total"]
        delta_prop_diff = prop_diff_mcm - prop_diff_ctrl

        results.append({"MCM": mcm, "Controle": ctrl, "simbolo_controle": control_symbol_map[ctrl], "Epimastigotas_MCM": row_mcm["Epimastigotas"], "Diferenciadas_MCM": row_mcm["Diferenciadas"], "Epimastigotas_Controle": row_ctrl["Epimastigotas"], "Diferenciadas_Controle": row_ctrl["Diferenciadas"], "prop_diferenciadas_MCM": prop_diff_mcm, "prop_diferenciadas_Controle": prop_diff_ctrl, "delta_prop_diferenciadas": delta_prop_diff, "odds_ratio_E_vs_D": odds_ratio, "p_value_raw": p_value })

## Correção de múltiplas comparações de Holm
fisher_df = pd.DataFrame(results)
reject, pvals_holm, _, _ = multipletests(fisher_df["p_value_raw"], alpha=0.05, method="holm")

fisher_df["p_value_holm"] = pvals_holm
fisher_df["significativo_holm"] = reject

# Interpretação simples
fisher_df["interpretacao"] = np.where(fisher_df["delta_prop_diferenciadas"] < 0, "MCM menor diferenciacao que o controle", np.where( fisher_df["delta_prop_diferenciadas"] > 0, "MCM maior diferenciacao que o controle", "Mesma proporcao observada"))

# Ordenação: 1) menor proporção de diferenciadas na MCM; 2) maior OR E vs D; 3) menor p corrigido
fisher_df = fisher_df.sort_values(
    by=["prop_diferenciadas_MCM", "odds_ratio_E_vs_D", "p_value_holm"],
    ascending=[True, False, True]).reset_index(drop=True)

fisher_df.to_csv("META_fisher_MCM_vs_controles.tsv", sep="\t", index=False)

print("\nFisher: MCM vs controles")
print(fisher_df[["MCM", "Controle", "prop_diferenciadas_MCM", "prop_diferenciadas_Controle", "delta_prop_diferenciadas", "odds_ratio_E_vs_D", "p_value_raw", "p_value_holm", "significativo_holm"]])

## Ranking das MCMs por inibição da diferenciação
# Ranking descritivo simples: menor proporção de diferenciadas
ranking_mcm = df[df["Tratamentos"].isin(mcm_list)].copy()
ranking_mcm = ranking_mcm.sort_values(by=["prop_diferenciadas", "odds_E_vs_D"], ascending=[True, False]).reset_index(drop=True)

ranking_mcm["ranking_menor_diferenciacao"] = np.arange(1, len(ranking_mcm) + 1)

ranking_mcm.to_csv("META_ranking_MCM_menor_diferenciacao.tsv", sep="\t", index=False)

print("\nRanking descritivo das MCMs (menor diferenciação):")
print(ranking_mcm[["Tratamentos", "Epimastigotas", "Diferenciadas", "prop_diferenciadas", "odds_E_vs_D"]])

# Tabela final resumida por MCM
summary_rows = []
for mcm in ranking_mcm["Tratamentos"]:
    sub = fisher_df[fisher_df["MCM"] == mcm].copy()

    # Quais controles foram significativos após Holm
    sub_sig = sub[sub["significativo_holm"]].copy()

    controles_sig = sub_sig["Controle"].tolist()
    simbolos_sig = sub_sig["simbolo_controle"].tolist()

    # texto resumido dos controles significativos
    if len(controles_sig) == 0:
        controles_sig_txt = "Nenhum"
        simbolos_sig_txt = ""
    else:
        controles_sig_txt = "; ".join(controles_sig)
        simbolos_sig_txt = "\n".join(simbolos_sig)  

    # Melhor comparação em termos de menor p ajustado
    best_row = sub.sort_values("p_value_holm", ascending=True).iloc[0]

    # Quantos controles a MCM superou no sentido biológico desejado
    n_menor_que_controles = (sub["delta_prop_diferenciadas"] < 0).sum()

    row_main = ranking_mcm[ranking_mcm["Tratamentos"] == mcm].iloc[0]

    summary_rows.append({"MCM": mcm, "ranking_menor_diferenciacao": int(row_main["ranking_menor_diferenciacao"]), "Epimastigotas": int(row_main["Epimastigotas"]), "Diferenciadas": int(row_main["Diferenciadas"]), "prop_diferenciadas": row_main["prop_diferenciadas"], "prop_epimastigotas": row_main["prop_epimastigotas"], "odds_E_vs_D": row_main["odds_E_vs_D"], "n_controles_com_menor_diferenciacao": int(n_menor_que_controles), "controles_significativos_holm": controles_sig_txt, "simbolos_para_grafico": simbolos_sig_txt, "melhor_controle_comparado": best_row["Controle"], "melhor_OR_E_vs_D": best_row["odds_ratio_E_vs_D"], "melhor_p_holm": best_row["p_value_holm"], "melhor_delta_prop_diferenciadas": best_row["delta_prop_diferenciadas"]})

summary_mcm_df = pd.DataFrame(summary_rows)

## Plot: proporção de diferenciadas
# Cor diferente para controles
colors = []
for trt in df["Tratamentos"]:
    if trt in controls:
        colors.append("#7a7a7a")
    else:
        colors.append("#ef80ac")

fig, ax = plt.subplots(figsize=(14,7))
ctrl_mean = df[df["Tratamentos"].isin(controls)]["prop_diferenciadas"].mean()
sns.barplot(data=df, x="Tratamentos", y="prop_diferenciadas", dodge=False, palette=colors, edgecolor="black", linewidth=1, ax=ax)

# Linha média dos controles
ax.axhline(ctrl_mean, color="black", linestyle="--", linewidth=2, label="Média dos controles")

for i, row in df.iterrows():
    if row["significativo"]:
        ax.text(i, row["prop_diferenciadas"] + 0.03, "***", ha="center", va="bottom", fontsize=16, fontweight="bold")
texto = (f"Qui-quadrado global\n"
    f"χ² = {chi2:.2f}\n"
    f"p = {p_chi2:.3e}")
ax.text(0.98, 0.98, texto, transform=ax.transAxes, ha="right", va="top", fontsize=11, bbox=dict(facecolor="white", edgecolor="black", boxstyle="round"))
ax.set_ylabel("Proporção de diferenciadas")
ax.set_xlabel("")
ax.set_ylim(0, 1)
plt.xticks(rotation=90)
sns.despine()
plt.tight_layout()
plt.savefig("META_grafico_proporcao_diferenciadas.png")

