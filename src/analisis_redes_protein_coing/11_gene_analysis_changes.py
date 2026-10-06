import pandas as pd
import numpy as np
import polars as pl
import networkx as nx
from pathlib import Path
from itertools import combinations

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = Path(f"{BASE}/output/protein_coding/mi_matrices")
COMUNIDADES= Path(f"{BASE}/output/protein_coding/data/comunidades_infomap")
OUTFIGDIR = Path(f"{BASE}/output/protein_coding/figures")
OUTDATADIR = Path(f"{BASE}/output/protein_coding/data")

CHANGE= (OUTDATADIR / "gene_degrees_all_groups_protein_coding.csv")

def read_filteres_graphs_pc():
    edgelists = {
        group: pl.read_parquet(
            INDIR / "filtered_edgelists" / f"Vip_{group}_protein_coding_mi_edgelist_filtered01.parquet"
        )
        for group in GROUPS
        }    
    graphs = {
        group: nx.from_pandas_edgelist(
            df.to_pandas(), source="source", target="target", edge_attr="weight"
        )
        for group, df in edgelists.items()
        }
    return graphs

def read_degree_change():
    df = pd.read_csv(CHANGE)
    return(df)

def read_comunidades(group):
    df_comunidades = pd.read_csv(f"{COMUNIDADES}/infomap_comunidades_{group}.csv")

    num_comunidades= df_comunidades["modulo"].unique()

    comunidades = {i: df_comunidades.loc[df_comunidades["modulo"] ==i, "gen"].tolist()
                   for i in num_comunidades}
    return(comunidades)

def grado_normalizado(df_degrees):
    df_degrees = df_degrees.drop(columns=["max_range"])
    nodes = {
        "Not_AD": 1005,
        "Low": 989,
        "Intermediate": 838,
        "High": 773
    }

    for group, n_nodes in nodes.items():
        df_degrees[group] = df_degrees[group] / n_nodes
        
    return(df_degrees)

def grado_intra_e_inter_comunidad(df_degrees, comunidades, G, group):
    
    genes_grado = df_degrees.loc[df_degrees[group] > 0, "gene"].unique().tolist()

    comunidad_por_gen = {
        gen: comunidad
        for comunidad, genes in comunidades.items()
        for gen in genes
    }

    resultados = {}

    for gen in genes_grado:
        intra = 0
        inter = 0
        conteo_comunidades = {}

        comunidad_gen = comunidad_por_gen[gen]
        grado_total = G.degree(gen)

        for vecino in G.neighbors(gen):
            comunidad_vecino = comunidad_por_gen[vecino]

            conteo_comunidades[comunidad_vecino] = (
                conteo_comunidades.get(comunidad_vecino, 0) + 1
            )

            if comunidad_vecino == comunidad_gen:
                intra += 1
            else:
                inter += 1

        if grado_total > 0:
            participacion = 1 - sum(
                (n / grado_total) ** 2
                for n in conteo_comunidades.values()
            )
        else:
            participacion = np.nan

        resultados[gen] = {
            "comunidad": comunidad_gen,
            "intra": intra,
            "inter": inter,
            "participacion": participacion
        }

    df_resultados = (
        pd.DataFrame.from_dict(resultados, orient="index")
        .rename_axis("gene")
        .reset_index()
    )

    media_intra = df_resultados.groupby("comunidad")["intra"].transform("mean")
    sd_intra = df_resultados.groupby("comunidad")["intra"].transform(
        lambda x: x.std(ddof=0)
    )

    df_resultados["z_intramodular"] = (
        (df_resultados["intra"] - media_intra)
        / sd_intra.replace(0, np.nan)
    )

    intra_por_gen = df_resultados.set_index("gene")["intra"]
    inter_por_gen = df_resultados.set_index("gene")["inter"]
    participacion_por_gen = df_resultados.set_index("gene")["participacion"]
    z_por_gen = df_resultados.set_index("gene")["z_intramodular"]

    return intra_por_gen, inter_por_gen, participacion_por_gen, z_por_gen
    
def gene_topological_profiles():
    degrees = read_degree_change()
    
    graphs= read_filteres_graphs_pc()
    
    df_degrees = grado_normalizado(degrees)
    
    for group in GROUPS:
        G = graphs[group]
        
        betweenness = nx.betweenness_centrality(G, normalized=True)    
        df_degrees[f"betweenness_{group}"] = df_degrees["gene"].map(betweenness)
        
        cc = nx.clustering(G)
        df_degrees[f"cc_{group}"] = df_degrees["gene"].map(cc)
        
        comunidades= read_comunidades(group)
        
        intra, inter, participacion, z = grado_intra_e_inter_comunidad(df_degrees, comunidades, G, group)

        df_degrees[f"intra_{group}"] = df_degrees["gene"].map(intra)
        df_degrees[f"inter_{group}"] = df_degrees["gene"].map(inter)
        df_degrees[f"participacion_{group}"] = df_degrees["gene"].map(participacion)
        df_degrees[f"z_intramodular_{group}"] = df_degrees["gene"].map(z)

        df_degrees[f"embeddedness_{group}"] = (df_degrees[f"intra_{group}"] / (df_degrees[f"intra_{group}"] + df_degrees[f"inter_{group}"]))
        
        df_degrees.to_csv(
            OUTDATADIR / "gene_topological_profiles_by_AD_stage_protein_coding.csv",
            index=False
        )
        
        return df_degrees
        
def gene_topological_analysis(df_degrees):
    
    metricas = [
        "grado_normalizado",
        "betweenness",
        "cc",
        "intra",
        "inter",
        "embeddedness",
        "participacion",
        "z_intramodular"
    ]
    
    for group1, group2 in combinations(GROUPS, 2):
        
        
        
def main():
    df= gene_topological_profiles()
    
    
        
if __name__ == "__main__":
    main()
    
    
        