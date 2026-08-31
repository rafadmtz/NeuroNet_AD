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

MI_CUTOFFS = dict(zip(GROUPS, [0.209, 0.216, 0.24, 0.278]))

def degree_variation_analysis():
    outdir = INDIR / "edgelist_vip_completo" / "gene_degrees"
    outdir.mkdir(parents=True, exist_ok=True)

    for group in GROUPS:
        edgelist = pl.read_parquet(
            INDIR / "edgelist_vip_completo" / f"edgelist_{group}_MI_nodelink01.parquet"
        )

        G = nx.from_pandas_edgelist(
            edgelist.to_pandas(), source="source", target="target", edge_attr="MI"
        )

        # universo completo de genes evaluados en la matriz MI original
        mat = pl.read_parquet(INDIR / "mi_matrices" / f"Vip_{group}_MI.parquet")
        id_col = mat.columns[0]
        all_genes = mat[id_col].to_list()

        G.add_nodes_from(all_genes)  # los que no tenían arista quedan con degree 0

        degree_data = [
            {"gene": gene, "degree": G.degree(gene)}
            for gene in G.nodes()
        ]

        df_degree = pd.DataFrame(degree_data)
        df_degree = df_degree.sort_values("degree", ascending=False)

        df_degree.to_csv(
            outdir / f"gene_degrees_{group}_01.csv",
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
            INDIR / "edgelist_vip_completo"/ "gene_degrees" / f"gene_degrees_{group}_01.csv",
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
        INDIR / "edgelist_vip_completo" / "gene_degrees" / "top_place_genes_vip01.csv"
    )
    
if __name__ == "__main__":
    top_place_genes(genes_list)