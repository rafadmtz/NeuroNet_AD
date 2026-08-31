import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import networkx as nx
import pandas as pd

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = BASE / "output" / "red_mi_vip_completo"
GROUPS = ["Not_AD", "Low", "Intermediate", "High"]


def edgelist_from_mi():
    for group in GROUPS:
        mat = pl.read_parquet(INDIR / "mi_matrices" / f"Vip_{group}_MI.parquet")
        id_col = mat.columns[0]

        edgelist = (
            mat.unpivot(index=id_col, variable_name="target", value_name="MI")
                .rename({id_col: "source"})
                .filter(pl.col("source") != pl.col("target"))
                .filter(pl.col("MI").is_not_null())
                .filter(pl.col("source") < pl.col("target"))
                .with_columns([
                    pl.col("source").cast(pl.Categorical),
                    pl.col("target").cast(pl.Categorical),
                    pl.col("MI").cast(pl.Float32),
            ])
        )
        edgelist.write_parquet(f"edgelist_{group}_MI_completo.parquet")
    
def node_to_link_01():
    for group in GROUPS:

        mi_cutoffs = []
        values_ratio = []
        values_nodes = []

        edgelist_full= pl.read_parquet(INDIR / "edgelist_vip_completo" / f"edgelist_{group}_MI_completo.parquet")
        
        for i in range(1000, 5000, 5):
            mi= i / 10000
            edgelist=edgelist_full.filter(pl.col("MI") > mi)

            nodes = pl.concat([edgelist["source"], edgelist["target"]]).n_unique()
            n_edges = len(edgelist)

            ratio= nodes / n_edges if n_edges > 0 else 0
            if 0.09 < ratio < 0.11:
                values_ratio.append(ratio)
                mi_cutoffs.append(str(mi))
                values_nodes.append(nodes)
                
        fig, ax1 = plt.subplots()

        ax1.bar(mi_cutoffs, values_ratio, color="steelblue", label="Node/Edge ratio")
        ax1.set_xlabel("MI cutoff")
        ax1.set_ylabel("Node/Edge ratio", color="steelblue")
        ax1.tick_params(axis="y", labelcolor="steelblue")
        ax1.set_xticklabels(mi_cutoffs, rotation=45, ha="right")

        ax2 = ax1.twinx()
        ax2.plot(mi_cutoffs, values_nodes, color="firebrick", marker="o", label="N° nodos")
        ax2.set_ylabel("Número de nodos", color="firebrick")
        ax2.tick_params(axis="y", labelcolor="firebrick")

        plt.title(f"Node-to-link ratio y N° de nodos por cutoff de MI — {group}")
        fig.tight_layout()
        plt.savefig(f"src/figures/Vip_{group}_close_to0.1.png", dpi=150, bbox_inches="tight")
        plt.close()
        
        
def node_to_link_01_write_parquet():
    mi_values01=[0.209, 0.216, 0.24, 0.278]
    for group, i in zip(GROUPS, mi_values01):
        edgelist_full= pl.read_parquet(INDIR / "edgelist_vip_completo" / f"edgelist_{group}_MI_completo.parquet")
        edgelist=edgelist_full.filter(pl.col("MI") > i)
        edgelist.write_parquet(INDIR / "edgelist_vip_completo" /f"edgelist_{group}_MI_nodelink01.parquet")


def network_metrics(G):
    components = list(nx.connected_components(G))
    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "connected_components": len(components),
        "clustering_coefficient": nx.average_clustering(G),
    }

def analyze_01_vip():
    results = {}
    for group in GROUPS:
        edgelist_full=pl.read_parquet(INDIR / "edgelist_vip_completo" /f"edgelist_{group}_MI_nodelink01.parquet")
        G = nx.from_pandas_edgelist(
            edgelist_full.to_pandas(),
            source="source",
            target="target",
            edge_attr="MI",   # conserva el peso MI como atributo del edge
        )
        results[group] = network_metrics(G)
    
    df_metrics = pd.DataFrame(results).T
    df_metrics.to_csv(INDIR / "edgelist_vip_completo" / "metrics_vip_01.csv")
    return df_metrics
        
def jaccard_nodes(G1, G2):
    n1, n2 = set(G1.nodes()), set(G2.nodes())
    return len(n1 & n2) / len(n1 | n2)

def jaccard_edges(G1, G2):
    e1 = set(frozenset(e) for e in G1.edges())
    e2 = set(frozenset(e) for e in G2.edges())
    return len(e1 & e2) / len(e1 | e2)

def jaccard_analysis():
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

    results = []
    seen = set()

    for group in GROUPS:
        for jac in GROUPS:
            if jac == group:
                continue
            pair = frozenset({group, jac})
            if pair in seen:
                continue
            seen.add(pair)

            results.append({
                "group_1": group,
                "group_2": jac,
                "jaccard_nodes": jaccard_nodes(graphs[group], graphs[jac]),
                "jaccard_edges": jaccard_edges(graphs[group], graphs[jac]),
            })

    df_jaccard = pd.DataFrame(results)
    df_jaccard.to_csv(INDIR / "edgelist_vip_completo" / "jaccard_analysis.csv", index=False)
    return df_jaccard
    
from venn import venn

def plot_venn_nodes():
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
    node_sets = {
        group: set(graphs[group].nodes())
        for group in GROUPS
    }

    fig, ax = plt.subplots(figsize=(10, 10))
    venn(node_sets, ax=ax)
    plt.title("Node composition overlap across groups")
    plt.savefig(BASE / "src" / "figures" / "venn_nodes.png", dpi=300, bbox_inches="tight")
    plt.close()
    
def core_genes():
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
    node_sets = {
        group: set(graphs[group].nodes())
        for group in GROUPS
    }
    
    cores=node_sets["Not_AD"]
    for group in GROUPS:
        cores &= node_sets[group]
        
    print(f"Core genes (presentes en los {len(GROUPS)} grupos): {len(cores)}")
    
    core_df = pd.DataFrame({"gene": sorted(cores)})
    core_df.to_csv(INDIR / "edgelist_vip_completo" / "core_genes.csv", index=False)
    return cores, graphs

def subgraph_coregenes():
    cores, graphs= core_genes()
    
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
    
    core_subgraphs = {
        group: graphs[group].subgraph(cores).copy()
        for group in GROUPS
    }
    
    summary = []
    for group, G_core in core_subgraphs.items():
        summary.append({
            "group": group,
            "core_nodes": G_core.number_of_nodes(),
            "core_edges": G_core.number_of_edges(),
            "core_connected_components": nx.number_connected_components(G_core),
            "core_clustering_coefficient": nx.average_clustering(G_core),
        })

    df_summary = pd.DataFrame(summary)
    df_summary.to_csv(INDIR / "edgelist_vip_completo" / "core_subgraphs_metrics.csv", index=False)
    print(df_summary)

    return core_subgraphs, df_summary

def degree_variation_analysis():
    cores, graphs = core_genes()

    # Degree de cada gen core en cada grupo (red completa, no el subgrafo)
    degree_data = []
    for gene in cores:
        row = {"gene": gene}
        for group in GROUPS:
            row[f"degree_{group}"] = graphs[group].degree(gene)
        degree_data.append(row)

    df_degree = pd.DataFrame(degree_data)

    # Variación: diferencia entre el degree máximo y mínimo de cada gen entre grupos
    degree_cols = [f"degree_{group}" for group in GROUPS]
    df_degree["max_degree"] = df_degree[degree_cols].max(axis=1)
    df_degree["min_degree"] = df_degree[degree_cols].min(axis=1)
    df_degree["degree_range"] = df_degree["max_degree"] - df_degree["min_degree"]

    # También el ratio, para capturar cambios relativos (ej. de 3 a 21 es más dramático que de 100 a 118)
    df_degree["degree_ratio"] = df_degree["max_degree"] / df_degree["min_degree"].replace(0, 1)

    df_degree = df_degree.sort_values("degree_range", ascending=False)

    df_degree.to_csv(INDIR / "edgelist_vip_completo" / "core_gene_degree_variation.csv", index=False)

    return df_degree

def plot_network_architecture():
    cores, graphs = core_genes()
    colors = {
        "Not_AD": "#8B5CF6",     # morado
        "Low": "#3B82F6",       # azul
        "Intermediate": "#10B981",  # verde
        "High": "#EF4444",      # rojo
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    axes = axes.flatten()

    for ax, group in zip(axes, GROUPS):
        G = graphs[group]  # tu grafo completo ya construido

        # Layout: spring_layout separa el componente gigante de las islas
        # k controla la separación entre nodos; con redes de ~1000 nodos, k pequeño ayuda
        pos = nx.spring_layout(G, k=0.15, iterations=50, seed=42)

        # Tamaño de nodo proporcional al degree, para resaltar hubs
        degrees = dict(G.degree())
        node_sizes = [10 + degrees[n] * 0.8 for n in G.nodes()]

        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_size=node_sizes,
            node_color=colors[group],
            alpha=0.7,
            linewidths=0.3,
            edgecolors="black",
        )
        nx.draw_networkx_edges(
            G, pos, ax=ax,
            alpha=0.3,
            width=2,
            edge_color="grey",
        )

        ax.set_title(
            f"{group}\n({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)",
            fontsize=14, fontweight="bold"
        )
        ax.axis("off")

    plt.suptitle("Network architecture across AD severity groups (VIP)", fontsize=18, y=0.98)
    plt.tight_layout()
    plt.savefig(BASE / "src" / "figures" / "network_architecture.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Guardado en src/figures/network_architecture.png")

if __name__ == "__main__":
    plot_network_architecture()