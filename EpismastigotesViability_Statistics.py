# Estatística para triagem das moléculas em epimastigotas
# Script criado com o auxílio do ChatGPT
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("paper")
sns.set_style("whitegrid", {"grid.color": "#F0F0F0"})

df = pd.read_csv("triagem_mols_todas_proporcional.tsv", sep="\t")

# Identify replicate columns
replicate_cols = [c for c in df.columns if "poço" in c]

# Convert comma decimals to floats
for col in replicate_cols + ["media"]:
    df[col] = df[col].astype(str).str.replace(",", ".").astype(float)

# Prepare long-format DataFrame
data_long = []
for i, mol in enumerate(df["Molecula"]):
    for val in df.loc[i, replicate_cols]:
        data_long.append((mol, val))
data_long = pd.DataFrame(data_long, columns=["Molecule", "Value"])

df["mean"] = df[replicate_cols].mean(axis=1)
df["std"] = df[replicate_cols].std(axis=1)

# Build figure
fig = sns.barplot(data=df, x="Molecula", y="mean", color="#ef80ac")
fig.tick_params(axis="x", rotation=45)

for i, row in df.iterrows():
    if row["std"] <= 0.05:
        texto = f"{row['std']:.2f}"
        x_pos = i
        y_pos = row["mean"] + row["std"] + 0.2
        fig.text( x=x_pos, y=y_pos, s=texto, ha="center", va="bottom", color="black", fontsize=10)

name_map = {
    "poço1 (Mi cél/mL)": "Triplicata 1", "poço2 (Mi cél/mL)": "Triplicata 2", "poço3 (Mi cél/mL)": "Triplicata 3"}

# Scatter points (replicates)
for col in replicate_cols:
    sns.scatterplot(data=df, x="Molecula", y=col, label=name_map[col], s=20,                 marker="o", edgecolor="black", linewidth=1 alpha=0.9, ax=fig)            

fig.errorbar(x=range(len(df)), y=df["mean"], yerr=df["std"], fmt="none", c="black", capsize=5)
fig.set_title("Ensaio de viabilidade celular de epimastigotas", pad=20, loc="center")  
fig.set_xlabel("Molécula", labelpad=15)
fig.set_ylabel("Contagem de epimastigotas (células * 10⁶/mL)", labelpad=15)
fig.tick_params(axis="both", labelsize=10)
fig.legend(title="Triplicata técnica", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, facecolor="white")

sns.despine(ax=fig)
fig.set_facecolor("white")

fig.figure.tight_layout()

plt.savefig("dd_wdesvio_padrao_epimastigota_todas_normalizado.png", dpi=300, bbox_inches="tight")
