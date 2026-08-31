import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
INDIR = Path("/export/space3/users/silvanac/NeuroNet_AD/output/red_mi_vip_completo")
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")

def read_graphs():
    edgelists = {
        group: pl.read_parquet(
            INDIR / "edgelist_vip_completo" / f"edgelist_{group}_MI_nodelink01.parquet"
        )
        for group in GROUPS
        }    
    graphs = {
        group: nx.from_pandas_edgelist(
            df.to_pandas(), source="source", target="target", edge_attr="MI"
        )
        for group, df in edgelists.items()
        }
    return graphs

def degree_variation_analysis():
    graphs = read_graphs() 
    
    for group in GROUPS:
        # Degree de cada gen core en cada grupo
        degree_data = []
        
        for gene in graphs[group].nodes():
            degree_data.append({
                "gene": gene,
                "degree": graphs[group].degree(gene)
            })

        df_degree = pd.DataFrame(degree_data)
        df_degree = df_degree.sort_values("degree", ascending=False)

        df_degree.to_csv(
            INDIR / "edgelist_vip_completo" / f"gene_degrees_{group}_01.csv",
            index=False
        )
        
        
genes_list = [
    "BIN1",
    "KCNQ2",
    "CHD5",
    "KDM4B",
    "SUGP2",
    "LUC7L",
    "ADGRL1",
    "SRCIN1",
    "STX16",
]

def top_place_genes(genes_list):
    
    degree_data = []
    
    for group in GROUPS:
        df_degree = pd.read_csv(
            INDIR / "edgelist_vip_completo" / f"gene_degrees_{group}_01.csv",
            index_col="gene"
        )
        
        for gene in genes_list:
            if gene in df_degree.index:
                place = df_degree.index.get_loc(gene) + 1
            else:
                place = None
            
            degree_data.append({
                            "gene": gene,
                            "group": group,
                            "place": place
                        })
    
    df_result = pd.DataFrame(degree_data)
    df_result = df_result.pivot(index="gene", columns="group", values="place")
    df_result = df_result.reindex(genes_list)  # conserva el orden original de genes_list
    df_result = df_result[GROUPS]  # conserva el orden original de GROUPS como columnas

    df_result.to_csv(
        INDIR / "edgelist_vip_completo" / "top_place_genes_vip01.csv"
    )
        
def shortest_path_analysis():
    graphs = read_graphs()
    results = []

    for group in GROUPS:
        G = graphs[group]

        # trabajar solo sobre el componente conexo más grande
        largest_cc = max(nx.connected_components(G), key=len)
        G_giant = G.subgraph(largest_cc)

        avg_path_length = nx.average_shortest_path_length(G_giant)

        results.append({
            "group": group,
            "giant_component_nodes": G_giant.number_of_nodes(),
            "avg_shortest_path_length": avg_path_length,
        })

    df_results = pd.DataFrame(results)
    df_results.to_csv(
        INDIR / "edgelist_vip_completo" / "avg_shortest_path_length.csv",
        index=False,
    )
    return df_results

def distribucion_grado_plot():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    colores = ["#1f77b4", "#f1c40f", "#e67e22", "#e21f0a"]

    for ax, group, color in zip(axes, GROUPS, colores):
        df = pd.read_csv(
            INDIR / "edgelist_vip_completo" / f"gene_degrees_{group}_01.csv"
        )

        ax.hist(df["degree"], bins=30, color=color, edgecolor="black")
        ax.set_title(group)
        ax.set_yscale("log")
        ax.set_xlabel("Degree")
        ax.set_ylabel("Number of genes")

    fig.suptitle("Degree distribution by AD severity group")
    fig.tight_layout()
    fig.savefig(BASE / "src" / "figures" /"degree_distribution_4groups_vip01.png", dpi=300)
    plt.close(fig)



def distribucion_grado_plot_log():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    colores = ["#1f77b4", "#f1c40f", "#e67e22", "#e21f0a"]

    for ax, group, color in zip(axes, GROUPS, colores):
        df = pd.read_csv(
            INDIR / "edgelist_vip_completo" / f"gene_degrees_{group}_01.csv"
        )

        degree_counts = df["degree"].value_counts().sort_index()

        ax.scatter(degree_counts.index, degree_counts.values, 
                   color=color, edgecolor="black", s=20)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(group)
        ax.set_xlabel("Degree (log)")
        ax.set_ylabel("Number of genes (log)")

    fig.suptitle("Degree distribution by AD severity group")
    fig.tight_layout()
    fig.savefig(BASE / "src" / "figures" /"degree_distribution_4groups_vip01_log.png", dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    distribucion_grado_plot()
    distribucion_grado_plot_log()
    