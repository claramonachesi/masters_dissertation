# Análise de viabilidade celular com correção por branco, normalização por Controle+DMSO, média das triplicatas técnicas, Friedman global e Wilcoxon com correção de Holm pareado vs Controle+DMSO.
# Script criado com auxílio do ChatGPT
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy.stats import friedmanchisquare, wilcoxon
from statsmodels.stats.multitest import multipletests

# teste com apenas as duas replicatas boas
files = {"BioRep2": "viabilidade_2026_UERJ_raw_96h_2.tsv", "BioRep3": "viabilidade_2026_UERJ_raw_96h_3.tsv"}
normalization_control = "Controle+DMSO" #Nome do grupo usado para normalização
blank_name = "Branco(Meio)" #Nome do branco usado para correção de fundo

# Ordem desejada no gráfico/tabelas
group_order = ["Branco(Meio)", "Controle", "Controle+DMSO", "Triton 1%", "LPC", "MCM1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MCM8", "MCM9", "MCM10", "MCM11", "MCM12", "MCM13", "MCM14", "MCM15", "MCM16", "MCM17", "MCM18"]

# Comparações planejadas vs Controle+DMSO
planned_comparisons = ["Controle", "Triton 1%", "LPC", "MCM1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "MCM8", "MCM9", "MCM10", "MCM11", "MCM12", "MCM13", "MCM14", "MCM15", "MCM16", "MCM17", "MCM18"]

# Apelido para ajudar o script a reconhecer nomes, caso variem um pouco
aliases = {"branco(meio)": "Branco(Meio)", "branco": "Branco(Meio)", "controle+dmso": "Controle+DMSO", "controle + dmso": "Controle+DMSO", "controle dmso": "Controle+DMSO", "controle": "Controle", "triton 1%": "Triton 1%", "lpc": "LPC",}

## FUNÇÕES AUXILIARES
def canonicalize_group_name(name: str) -> str:
    """Padroniza o nome da condição."""
    s = str(name).strip()
    s_low = s.lower().replace("µ", "u")
    return aliases.get(s_low, s)

def convert_numeric_columns(df: pd.DataFrame, columns: list) -> pd.DataFrame:
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors="raise")
    return df

def load_one_biorep(file_path: str, biorep_name: str) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep="\t")

    # cada coluna = um grupo; cada linha = uma réplica técnica
    df = convert_numeric_columns(df, df.columns.tolist()).copy()

    # transpõe: índice vira nome do grupo
    df_t = df.T.reset_index()
    df_t.columns = ["Group", "tech_rep1", "tech_rep2", "tech_rep3"]

    df_t["Group"] = df_t["Group"].apply(canonicalize_group_name)
    df_t["BioRep"] = biorep_name
    df_t["tech_mean_raw"] = df_t[["tech_rep1", "tech_rep2", "tech_rep3"]].mean(axis=1)
    df_t["tech_sd_raw"] = df_t[["tech_rep1", "tech_rep2", "tech_rep3"]].std(axis=1, ddof=1)

    return df_t[["BioRep", "Group", "tech_rep1", "tech_rep2", "tech_rep3", "tech_mean_raw", "tech_sd_raw"]]

def process_one_biorep(df_bio: pd.DataFrame, blank_name: str = blank_name, norm_control: str = normalization_control) -> pd.DataFrame:
    """ corrige pelo branco e normaliza por Controle+DMSO """
    df = df_bio.copy()

    # Média técnica do branco
    blank_rows = df[df["Group"] == blank_name]
    if blank_rows.empty:
        raise ValueError(f"Não encontrei o grupo '{blank_name}' neste experimento.")
    blank_mean = blank_rows["tech_mean_raw"].mean()

    # Correção pelo branco
    df["blank_mean_raw"] = blank_mean
    df["tech_mean_blankcorr"] = df["tech_mean_raw"] - blank_mean

    # Controle de normalização
    ctrl_rows = df[df["Group"] == norm_control]
    if ctrl_rows.empty:
        raise ValueError(f"Não encontrei o grupo '{norm_control}' neste experimento.")
    ctrl_mean_blankcorr = ctrl_rows["tech_mean_blankcorr"].mean()
    if np.isclose(ctrl_mean_blankcorr, 0):
        raise ValueError(f"O valor corrigido do grupo '{norm_control}' ficou zero; não é possível normalizar.")

    # Viabilidade relativa ao Controle+DMSO
    df["norm_control_blankcorr"] = ctrl_mean_blankcorr
    df["viability_pct"] = (df["tech_mean_blankcorr"] / ctrl_mean_blankcorr) * 100.0

    return df

def add_significance_annotations(fig, stats_df, x_order, y_base, y_step):
    """Adiciona barras/anotações de significância ao gráfico."""
    sig_df = stats_df[stats_df["p_holm"] < 0.05].copy()
    current_y = y_base

    for _, row in sig_df.iterrows():
        g1 = row["Group_A"]
        g2 = row["Group_B"]
        if g1 not in x_order or g2 not in x_order:
            continue
        x1 = x_order.index(g1)
        x2 = x_order.index(g2)
        if x1 > x2:
            x1, x2 = x2, x1
        p = row["p_holm"]
        stars = "***" if p < 0.001 else "**" if p < 0.01 else "*"

        fig.add_shape(type="line", x0=x1, x1=x2, y0=current_y, y1=current_y, line=dict(color="black", width=1.5))
        fig.add_shape(type="line", x0=x1, x1=x1, y0=current_y, y1=current_y - y_step * 0.25, line=dict(color="black", width=1.5))
        fig.add_shape(type="line", x0=x2, x1=x2, y0=current_y, y1=current_y - y_step * 0.25, line=dict(color="black", width=1.5))
        fig.add_annotation(x=(x1 + x2) / 2, y=current_y + y_step * 0.08, text=f"{stars}<br>adj p = {p:.3f}", showarrow=False, font=dict(size=18, color="black"),align="center")
        current_y += y_step

    return fig

## LEITURA E PROCESSAMENTO DOS 3 EXPERIMENTOS
all_bio = []
for biorep_name, file_path in files.items():
    df_bio = load_one_biorep(file_path, biorep_name)
    df_bio = process_one_biorep(df_bio)
    all_bio.append(df_bio)

data = pd.concat(all_bio, ignore_index=True)
data.to_csv("VIAB2_raw_dados_processados_por_bioreplicata.tsv", sep="\t", index=False)

## RESUMO FINAL PARA PLOTAGEM (MÉDIA E DP BIOLÓGICOS)
summary = (data.groupby("Group", observed=False)["viability_pct"].agg(["mean", "std", "median", "count"]).reset_index().rename(columns={"mean": "bio_mean_pct", "std": "bio_sd_pct", "median": "bio_median_pct", "count": "n_biological"}))
summary.to_csv("VIAB2_raw_resumo_viabilidade_normalizada.tsv", sep="\t", index=False)

## FRIEDMAN GLOBAL 
# Para o teste global, remove o Branco por padrão.
groups_for_friedman = [g for g in group_order if g in data["Group"].astype(str).unique() and g != blank_name]

pivot_global = (data[data["Group"].isin(groups_for_friedman)].pivot_table(index="BioRep", columns="Group", values="viability_pct", aggfunc="mean"))

# Garante que todas as colunas existam e na ordem desejada
pivot_global = pivot_global[[g for g in groups_for_friedman if g in pivot_global.columns]]

# Remove colunas com algum NA, caso exista problema de ausência
pivot_global = pivot_global.dropna(axis=1, how="any")

if pivot_global.shape[1] < 3:
    raise ValueError("O teste de Friedman precisa de pelo menos 3 grupos completos. Verifique se todas as condições existem em todos os 3 experimentos biológicos.")

friedman_stat, friedman_p = friedmanchisquare(*[pivot_global[col].values for col in pivot_global.columns])
friedman_df = pd.DataFrame({"test": ["Friedman"], "n_biological_replicates": [pivot_global.shape[0]], "n_groups": [pivot_global.shape[1]], "statistic": [friedman_stat], "p_value": [friedman_p]})
friedman_df.to_csv("VIAB2_raw_friedman_global.tsv", sep="\t", index=False)

## WILCOXON PAREADO VS CONTROLE+DMSO
# Comparações planejadas
available_groups = set(data["Group"].astype(str).unique())
planned = [g for g in planned_comparisons if g in available_groups]

pivot_posthoc = (data[data["Group"].isin(planned + [normalization_control])].pivot_table(index="BioRep", columns="Group", values="viability_pct", aggfunc="mean"))

results = []
for group in planned:
    if group == normalization_control:
        continue
    if group not in pivot_posthoc.columns or normalization_control not in pivot_posthoc.columns:
        continue
    sub = pivot_posthoc[[normalization_control, group]].dropna()

    # Com n pequeno, zero differences podem atrapalhar. Guardamos as diferenças.
    diffs = sub[group] - sub[normalization_control]
    nonzero_diffs = diffs[~np.isclose(diffs, 0)]

    # Se todas as diferenças forem zero, Wilcoxon não é informativo
    if len(nonzero_diffs) == 0:
        stat = np.nan
        p = 1.0
        note = "Todas as diferenças pareadas foram zero."
    else:
        try:
            # zero_method='wilcox' ignora diferenças exatamente zero
            stat, p = wilcoxon(sub[group], sub[normalization_control], zero_method="wilcox", alternative="two-sided", mode="auto")
            note = ""
        except ValueError as e:
            stat = np.nan
            p = np.nan
            note = str(e)

    results.append({"Group_A": normalization_control, "Group_B": group, "n_pairs": len(sub), "median_A": sub[normalization_control].median(), "median_B": sub[group].median(), "mean_A": sub[normalization_control].mean(), "mean_B": sub[group].mean(), "sd_A": sub[normalization_control].std(ddof=1), "sd_B": sub[group].std(ddof=1), "mean_difference_B_minus_A": diffs.mean(), "median_difference_B_minus_A": diffs.median(), "wilcoxon_stat": stat, "p_raw": p, "note": note})

wilcox_df = pd.DataFrame(results)

# Correção de múltiplas comparações (Holm)
mask_valid = wilcox_df["p_raw"].notna()
rej, p_holm, _, _ = multipletests(wilcox_df.loc[mask_valid, "p_raw"], alpha=0.05, method="holm")

wilcox_df["p_holm"] = np.nan
wilcox_df["reject_holm_0.05"] = False
wilcox_df.loc[mask_valid, "p_holm"] = p_holm
wilcox_df.loc[mask_valid, "reject_holm_0.05"] = rej
wilcox_df.to_csv("VIAB2_raw_wilcoxon_vs_controle_dmso.tsv", sep="\t", index=False)

## GRÁFICO
plot_df = summary.copy()
plot_df = plot_df[plot_df["Group"].isin([g for g in group_order if g in plot_df["Group"].astype(str).tolist()])]
plot_df["Group"] = pd.Categorical(plot_df["Group"], categories=group_order, ordered=True)
plot_df = plot_df.sort_values("Group").reset_index(drop=True)
fig = go.Figure()

# Barras com média ± DP biológico
fig.add_trace(go.Bar(x=plot_df["Group"], y=plot_df["bio_mean_pct"], error_y=dict(type="data", array=plot_df["bio_sd_pct"].fillna(0), visible=True), marker_color="#ef80ac", showlegend=False, name="Média ± DP biológico"))

# Pontos individuais das réplicas biológicas
bio_colors = {"BioRep1": "black", "BioRep2": "black", "BioRep3": "black",}

for biorep in data["BioRep"].unique():
    sub = data[data["BioRep"] == biorep].copy()
    sub = sub.sort_values("Group")
    fig.add_trace(go.Scatter(x=sub["Group"], y=sub["viability_pct"], mode="markers", name=biorep, marker=dict(size=9, symbol="circle", line=dict(width=1, color="black")), opacity=0.9))

# Anotação do Friedman
fig.add_annotation(text=f"Friedman: χ² = {friedman_stat:.3f}, p = {friedman_p:.4f}", xref="paper", yref="paper", x=0.5, y=1.08, showarrow=False, font=dict(size=22, color="black"))

# Legenda de significância
fig.add_annotation(text="Wilcoxon pareado vs Controle+DMSO<br>* adj p < 0.05<br>** adj p < 0.01<br>*** adj p < 0.001", xref="paper", yref="paper", x=1.14, y=0.02, showarrow=False, align="left", font=dict(size=16, color="black"), bordercolor="black", borderwidth=1, borderpad=6, bgcolor="white")

# Barras de significância
y_max = np.nanmax(plot_df["bio_mean_pct"] + plot_df["bio_sd_pct"].fillna(0))
if np.isnan(y_max):
    y_max = np.nanmax(plot_df["bio_mean_pct"])

y_base = y_max * 1.10
y_step = max(5, y_max * 0.08)

fig = add_significance_annotations(fig=fig, stats_df=wilcox_df, x_order=[g for g in plot_df["Group"].astype(str).tolist()], y_base=y_base, y_step=y_step)
fig.update_layout(title="Ensaio de viabilidade celular", title_x=0.5, font_family="Manrope", title_font_size=34, xaxis_title="Condição", yaxis_title="Viabilidade (% do Controle+DMSO, após correção pelo branco)", template="simple_white", font=dict(size=20), legend_title="Réplicas biológicas", margin=dict(b=120, t=130),)
fig.update_xaxes(tickangle=45)
fig.write_image("VIAB2_raw_grafico_viabilidade_friedman_wilcoxon.png", scale=3, width=1600, height=1100)
fig.write_html("VIAB2_raw_grafico_viabilidade_friedman_wilcoxon.html")
