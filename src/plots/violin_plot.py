import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.makedirs("figures", exist_ok=True)

# Cargar tabla
df = pd.read_csv("data/conteocell_todos.tsv", sep="\t")

# Todos los subtipos
subtipos = [
    "Astrocyte", "Chandelier", "Endothelial",
    "L2-3-IT", "L4-IT", "L5-6-NP", "L5-ET", "L5-IT",
    "L6b", "L6-CT", "L6-IT-Car3", "L6-IT",
    "Lamp5-Lhx6", "Lamp5", "Microglia-PVM",
    "Oligodendrocyte", "OPC", "Pax6", "Pvalb",
    "Sncg", "Sst-Chodl", "Sst", "Vip", "VLMC"
]

# Filtrar donantes con todos ceros
df = df[df[subtipos].sum(axis=1) > 0]

# Ordenar por condición AD
orden_ad = ["Not AD", "Low", "Intermediate", "High"]
df["AD_State"] = pd.Categorical(df["AD_State"], categories=orden_ad, ordered=True)

# Long format
df_long = df.melt(
    id_vars=["ID", "AD_State", "Sex"],
    value_vars=subtipos,
    var_name="Subtipo",
    value_name="Num_celulas"
)

df_long = df_long[df_long["Num_celulas"] > 0]

# Colores por condición
colores_ad = {"Not AD": "#2ecc71", "Low": "#f1c40f", "Intermediate": "#e67e22", "High": "#e74c3c"}

# Figura horizontal
fig, ax = plt.subplots(figsize=(18, 8))  # ancho > alto para horizontal

sns.violinplot(
    data=df_long,
    x="Subtipo", y="Num_celulas",    # x y y volteados
    hue="AD_State", hue_order=orden_ad,
    palette=colores_ad,
    order=subtipos,
    inner="quartile",
    density_norm="width",
    cut=0,
    linewidth=0.8,
    ax=ax
)

ax.set_ylabel("Número de células", fontsize=12)
ax.set_xlabel("Subtipo celular", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)  # rotar etiquetas para que no se empalmen

ax.axvline(x=0, color="gray", linewidth=0.5, linestyle="--")  # línea en 0 de referencia
ax.legend(title="AD State", bbox_to_anchor=(1.01, 1), loc="upper left")

plt.tight_layout()
plt.savefig("figures/violin_subtipos.png", bbox_inches="tight", dpi=150)
print("Listo")