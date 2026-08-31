import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.makedirs("figures", exist_ok=True)

# Cargar tabla
df = pd.read_csv("data/conteocell_todos.tsv", sep="\t")

# Solo estos subtipos
subtipos = ["Lamp5", "Sst-Chodl", "Vip", "Sncg", "Lamp5-Lhx6"]
df = df[df[subtipos].sum(axis=1) > 0]

# Ordenar por condición AD
orden_ad = ["Not AD", "Low", "Intermediate", "High"]
df["AD_State"] = pd.Categorical(df["AD_State"], categories=orden_ad, ordered=True)



# Pasar a formato largo (long format) para seaborn
df_long = df.melt(
    id_vars=["ID", "AD_State", "Sex"],
    value_vars=subtipos,
    var_name="Subtipo",
    value_name="Num_celulas"
)



# Colores por condición
colores_ad = {"Not AD": "#2ecc71", "Low": "#f1c40f", "Intermediate": "#e67e22", "High": "#e74c3c"}

# Figura
fig, ax = plt.subplots(figsize=(14, 7))

# Boxplot + stripplot encima
sns.boxplot(
    data=df_long,
    x="Subtipo", y="Num_celulas",
    hue="AD_State", hue_order=orden_ad,
    palette=colores_ad,
    width=0.6,
    fliersize=0,   # ocultamos outliers del box porque los muestra el strip
    ax=ax
)

sns.stripplot(
    data=df_long,
    x="Subtipo", y="Num_celulas",
    hue="AD_State", hue_order=orden_ad,
    palette=colores_ad,
    dodge=True,
    size=4, alpha=0.6, jitter=True,
    legend=False,
    ax=ax
)

ax.set_xlabel("Subtipo celular", fontsize=12)
ax.set_ylabel("Número de células", fontsize=12)
ax.set_title("Distribución de células por subtipo y condición AD", fontsize=13)
ax.legend(title="AD State", bbox_to_anchor=(1.01, 1), loc="upper left")

plt.tight_layout()
plt.savefig("figures/boxplot_subtipos.png", bbox_inches="tight", dpi=150)
print("Listo")