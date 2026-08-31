import polars as pl
import networkx as nx
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import gseapy as gp

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = BASE / "output" / "red_mi_vip_completo" / "edgelist_vip_completo"
GROUPS = ["Not_AD", "Low", "Intermediate", "High"]

def enrichment_analysis():
    for group in GROUPS:
        edgelist01 = pd.read_csv(INDIR / f"infomap_{group}.csv")
        all_genes = edgelist01["gene"].to_list()

        genes_por_comunidad = edgelist01.groupby("community")["gene"].apply(list).to_dict()

        df_results = []

        for comunidad, genes in genes_por_comunidad.items():
            try:
                enr = gp.enrichr(
                    gene_list  = genes,
                    gene_sets  = "GO_Biological_Process_2023",
                    organism   = "human",
                    cutoff     = 0.05,
                    background = all_genes
                )
                df_results.append(
                    enr.results.assign(Comunidad=comunidad, tamaño_comunidad=len(genes))
                )
            except Exception as e:
                print(f"[{group}] comunidad {comunidad} (n={len(genes)} genes) falló: {e}")

        todos_resultados = pd.concat(df_results, ignore_index=True)
        todos_resultados = todos_resultados.sort_values("Adjusted P-value")
        todos_resultados.to_csv(INDIR / f"resultados_enrichment_vip_{group}_nl01_BP.csv", index=False)
        print(f"resultados_enrichment_vip_{group}_nl01_BP.csv")
            
import matplotlib.pyplot as plt

def distribuciones_comunidades():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.flatten()

    colores = ["#1f77b4", "#f1c40f", "#e67e22", "#e74c3c"]

    comunidades_por_grupo = {}

    for ax, group, color in zip(axes, GROUPS, colores):
        edgelist01 = pd.read_csv(INDIR / f"infomap_{group}.csv")

        # tamaño de cada comunidad
        genes_por_comunidad = edgelist01.groupby("community")["gene"].apply(list).to_dict()
        comunidades_por_grupo[group] = genes_por_comunidad

        tamanos = [len(genes) for genes in genes_por_comunidad.values()]

        ax.hist(tamanos, bins=30, color=color, edgecolor="black")
        ax.set_title(f"{group} (n={len(tamanos)} comunidades)")
        ax.set_xlabel("Community size (# genes)")
        ax.set_ylabel("Number of communities")

    fig.suptitle("Community size distribution by AD severity group")
    fig.tight_layout()
    fig.savefig(BASE / "src" / "figures" / "community_size_distribution_4groups_vip01.png", dpi=300)
    plt.close(fig)

    return comunidades_por_grupo

if __name__ == "__main__":
    distribuciones_comunidades()
    