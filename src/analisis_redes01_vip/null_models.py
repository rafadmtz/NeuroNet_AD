import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import numpy as np

from grado_nodos_vip01 import read_graphs


GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
INDIR = Path("/export/space3/users/silvanac/NeuroNet_AD/output/red_mi_vip_completo")
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")

def null_model_comparison():
    results = []
    graphs = read_graphs()
    for group in GROUPS:
        G = graphs[group]
        G_random = nx.gnm_random_graph(n=G.number_of_nodes(), m=G.number_of_edges())
        results.append({
            "group": group,
            "cc_real": nx.average_clustering(G),
            "cc_random": nx.average_clustering(G_random),
        })
    return pd.DataFrame(results)

def null_model_shortest_path():
    graphs = read_graphs()
    results = []

    for group in GROUPS:
        G = graphs[group]

        # giant component real, solo para sacar n y m de referencia
        largest_cc_real = max(nx.connected_components(G), key=len)
        G_giant_real = G.subgraph(largest_cc_real)
        n_giant = G_giant_real.number_of_nodes()
        m_giant = G_giant_real.number_of_edges()

        # --- red aleatoria con el mismo n y m del giant component real ---
        G_random = nx.gnm_random_graph(n=n_giant, m=m_giant)

        largest_cc_random = max(nx.connected_components(G_random), key=len)
        G_giant_random = G_random.subgraph(largest_cc_random)

        avg_path_random = nx.average_shortest_path_length(G_giant_random)
        cc_random = nx.average_clustering(G_giant_random)

        results.append({
            "group": group,
            "n_input": n_giant,
            "m_input": m_giant,
            "n_giant_random": G_giant_random.number_of_nodes(),
            "pct_nodes_retained_random": G_giant_random.number_of_nodes() / n_giant,
            "avg_path_random": avg_path_random,
            "cc_random": cc_random,
        })

    df_results = pd.DataFrame(results)
    df_results.to_csv(
        INDIR / "edgelist_vip_completo" / "null_model_shortest_path_vip01_maxcomponent.csv",
        index=False,
    )
    return df_results

if __name__ == "__main__":
    null_model_shortest_path()
