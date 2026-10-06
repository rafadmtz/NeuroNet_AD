import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = Path(f"{BASE}/output/protein_coding/mi_matrices")
OUTFIGDIR = Path(f"{BASE}/output/protein_coding/figures")
OUTDATADIR = Path(f"{BASE}/output/protein_coding/data")

CHANGE= (OUTDATADIR / "gene_degrees_all_groups_protein_coding.csv")

def overall_change():
    df_degree = pd.read_csv(CHANGE)

    overall_change_in_degree={}
    
    overall_change_in_degree["degree0_across"]= len(df_degree[df_degree["max_range"] == 0])
    print(f"genes 0 across {overall_change_in_degree['degree0_across']}")

    same_nonzero_count = len(df_degree[(df_degree["max_range"]==0) & (df_degree["Not_AD"]!=0)])
    if same_nonzero_count:
        overall_change_in_degree["same_degree_across"] = same_nonzero_count
        print(df_degree[(df_degree["max_range"]==0) & (df_degree["Not_AD"]!=0)][["gene"] + GROUPS])

    diffs = df_degree[GROUPS].diff(axis=1).iloc[:, 1:]

    overall_change_in_degree["increase"]= len(df_degree[(diffs > 0).any(axis=1)])
    print(f"genes increase {overall_change_in_degree['increase']}")
    
    overall_change_in_degree["decrease"] = len(df_degree[(diffs < 0).any(axis=1)])
    print(f"genes decrease {overall_change_in_degree['decrease']}")

    plt.figure(figsize=(6, 4))
    categories = list(overall_change_in_degree.keys())
    values = list(overall_change_in_degree.values())

    plt.bar(categories, values, color=["gray", "steelblue", "firebrick", "seagreen"])
    plt.ylabel("Número de genes")
    plt.title("Cambio global de grado a través de los estadios")
    plt.xticks(rotation=15)

    plt.tight_layout()
    plt.savefig(OUTFIGDIR / "overall_degree_change_counts.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()

    print(f"Guardado en {OUTFIGDIR / 'overall_degree_change_counts.png'}")
    
def plot_change():
    
    df_degree = pd.read_csv(CHANGE)

    df_degree = df_degree[df_degree["max_range"] > 0]

    max_range_genes= 40
    
    # x genes with more change
    df_degree = df_degree.nlargest(max_range_genes, "max_range")

    heatmap_data = df_degree.set_index("gene")[
    ["Not_AD", "Low", "Intermediate", "High"]
    ]

    heatmap_z = heatmap_data.apply(
        lambda x: (x - x.mean()) / x.std(),
        axis=1
    )

    # Tamaño de los números dentro del heatmap
    ANNOT_SIZE = 6

    plt.figure(figsize=(16, 4))

    sns.heatmap(
        heatmap_z.T,
        cmap="coolwarm",
        center=0,
        annot=True,
        fmt=".2f",
        annot_kws={"size": ANNOT_SIZE}
    )

    plt.xlabel(f"Top Genes by max_range {max_range_genes}")
    plt.ylabel("AD severity")

    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig(OUTFIGDIR / "heatmap_plt_change_z.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()
    print(f"Guardado en {OUTFIGDIR / 'heatmap_plt_change_z.png'}")
    
def monotonic_change_in_degree():
    
    df_degree = pd.read_csv(CHANGE)
    
    diffs = df_degree[GROUPS].diff(axis=1).iloc[:, 1:]

    is_increasing = (diffs > 0).all(axis=1)
    is_decreasing = (diffs < 0).all(axis=1)

    df_degree["monotonic"] = "no"
    df_degree.loc[is_increasing, "monotonic"] = "increasing"
    df_degree.loc[is_decreasing, "monotonic"] = "decreasing"

    df_monotonic = df_degree[df_degree["monotonic"] != "no"].copy()

    df_monotonic = df_monotonic.sort_values("max_range", ascending=False)

    outpath = OUTDATADIR / "monotonic_degree_changes_protein_coding.csv"
    df_monotonic.to_csv(outpath, index=False)

    print(f"Genes con cambio monotónico: {len(df_monotonic)} "
          f"(increasing: {is_increasing.sum()}, decreasing: {is_decreasing.sum()})")
    print(f"Guardado en {outpath}")

    return df_monotonic
        
def plot_monotonic_changes(df_monotonic):
    
    type_of_change= ["increasing", "decreasing"]
    
    max_range_genes= 20

    for change in type_of_change:
        df_change= df_monotonic[df_monotonic["monotonic"] == change]
        
        df_change = df_change.nlargest(max_range_genes, "max_range")
        
        heatmap_data = df_change.set_index("gene")[
            ["Not_AD", "Low", "Intermediate", "High"]
            ]
        
        heatmap_z = heatmap_data.apply(
            lambda x: (x - x.mean()) / x.std(),
            axis=1
        )
    
        # Tamaño de los números dentro del heatmap
        ANNOT_SIZE = 6
    
        plt.figure(figsize=(16, 4))
    
        sns.heatmap(
            heatmap_z.T,
            cmap="coolwarm",
            center=0,
            annot=True,
            fmt=".2f",
            annot_kws={"size": ANNOT_SIZE}
        )
    
        plt.xlabel(f"Top Genes by max_range {max_range_genes}, monotonic {change}")
        plt.ylabel("AD severity")
    
        plt.xticks(rotation=90)
    
        plt.tight_layout()
        plt.savefig(OUTFIGDIR / f"heatmap_plt_change_z_monotonic_{change}.png", dpi=300, bbox_inches="tight")
        plt.show()
        plt.close()
        print(f"Guardado en {OUTFIGDIR / f'heatmap_plt_change_z_monotonic_{change}.png'}")

    
if __name__ == "__main__":
    plot_change()
    overall_change()
    plot_monotonic_changes(monotonic_change_in_degree())