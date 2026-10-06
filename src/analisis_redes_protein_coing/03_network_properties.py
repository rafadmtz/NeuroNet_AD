import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations


GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = Path(f"{BASE}/output/protein_coding/mi_matrices")
OUTFIGDIR = Path(f"{BASE}/output/protein_coding/figures")
OUTDATADIR = Path(f"{BASE}/output/protein_coding/data")

OUTFIGDIR.mkdir(parents=True, exist_ok=True)
OUTDATADIR.mkdir(parents=True, exist_ok=True)

cutoffs = {
    "Not_AD": 0.2086,
    "Low": 0.2161,
    "Intermediate": 0.2401,
    "High": 0.2772
}

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


def mi_distribution():
    colores = ["#1f77b4", "#f1c40f", "#e67e22", "#e74c3c"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    for ax, group, color in zip(axes, GROUPS, colores):
        edgelist = pl.read_parquet(INDIR / f"Vip_{group}_protein_coding_mi_edgelist.parquet")
        ax.hist(edgelist["weight"].to_numpy(), bins=100, color=color, edgecolor="black", linewidth=0.3)
        ax.set_yscale("log")
        ax.set_title(group)
        ax.set_xlabel("MI")
        ax.set_ylabel("Frecuencia")

    fig.suptitle("Distribución de valores de MI por grupo de severidad")
    fig.tight_layout()
    fig.savefig(OUTFIGDIR / "mi_distribution_4groups_protein_coding.png", dpi=300)
    plt.close(fig)
        
        


def node_to_link_01():
    colores = ["#1f77b4", "#f1c40f", "#e67e22", "#e74c3c"]  # Not_AD, Low, Intermediate, High

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for ax1, group, color in zip(axes, GROUPS, colores):

        mi_cutoffs = []
        values_ratio = []
        values_nodes = []
        cutoffs_optimos = {}

        edgelist_full = pl.read_parquet(INDIR / f"Vip_{group}_protein_coding_mi_edgelist.parquet")

        for i in range(2000, 2800, 1):
            mi = i / 10000
            edgelist = edgelist_full.filter(pl.col("weight") > mi)

            nodes = pl.concat([edgelist["source"], edgelist["target"]]).n_unique()
            n_edges = len(edgelist)

            ratio = nodes / n_edges if n_edges > 0 else 0
            if 0.098 < ratio < 0.102:
                values_ratio.append(ratio)
                mi_cutoffs.append(str(mi))
                values_nodes.append(nodes)


        values_ratio_arr = np.array(values_ratio)
        distancias = np.abs(values_ratio_arr - 0.10)
        idx = np.argmin(distancias)

        cutoffs_optimos[group] = {
            "mi_cutoff": float(mi_cutoffs[idx]),
            "ratio": values_ratio[idx],
            "n_nodos": values_nodes[idx],
        }

        for group, info in cutoffs_optimos.items():
            print(f"{group}: MI={info['mi_cutoff']:.4f}  ratio={info['ratio']:.4f}  nodos={info['n_nodos']}")

        ax1.bar(mi_cutoffs, values_ratio, color=color, label="Node/Edge ratio")
        ax1.set_xlabel("MI cutoff")
        ax1.set_ylabel("Node/Edge ratio", color=color)
        ax1.tick_params(axis="y", labelcolor=color)
        ax1.set_xticks(range(len(mi_cutoffs)))
        ax1.set_xticklabels(mi_cutoffs, rotation=45, ha="right")

        ax2 = ax1.twinx()
        ax2.plot(mi_cutoffs, values_nodes, color="firebrick", marker="o", label="N° nodos")
        ax2.set_ylabel("Número de nodos", color="firebrick")
        ax2.tick_params(axis="y", labelcolor="firebrick")
        ax2.ticklabel_format(axis="y", useOffset=False, style="plain")

        ax1.set_title(group)

    fig.suptitle("Node-to-link ratio y N° de nodos por cutoff de MI")
    fig.tight_layout()
    fig.savefig(OUTFIGDIR / "node_to_link_01_4groups.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def degree_variation_analysis(cutoffs, graphs):

    for group in GROUPS:
        edgelist = pl.read_parquet(INDIR / f"Vip_{group}_protein_coding_mi_edgelist.parquet")

        todos = set(pl.concat([edgelist["source"], edgelist["target"]]).unique().to_list())

        sobreviven_df = edgelist.filter(pl.col("weight") > cutoffs[group])
        sobreviven = set(
            pl.concat([sobreviven_df["source"], sobreviven_df["target"]]).unique().to_list()
        )

        zero_degree = todos - sobreviven

        degree_data = [
            {"gene": gene, "degree": graphs[group].degree(gene)}
            for gene in graphs[group].nodes()
        ]
        df_degree = pd.DataFrame(degree_data)

        df_zero = pd.DataFrame({
            "gene": sorted(zero_degree),
            "degree": 0,
        })

        df_degree = pd.concat([df_degree, df_zero], ignore_index=True)
        df_degree = df_degree.sort_values("degree", ascending=False).reset_index(drop=True)

        nodos_grafo = set(graphs[group].nodes())
        print(f"{group}: {len(nodos_grafo)} nodos en la red + {len(zero_degree)} de grado 0 = {len(df_degree)} (universo: {len(todos)})")
        if nodos_grafo != sobreviven:
            print(f"  AVISO: el grafo no coincide con el cutoff {cutoffs[group]} "
                  f"({len(nodos_grafo - sobreviven)} nodos extra, {len(sobreviven - nodos_grafo)} faltantes)")
        if nodos_grafo & zero_degree:
            print(f"  AVISO: {len(nodos_grafo & zero_degree)} genes aparecen en ambos conjuntos")

        df_degree.to_csv(OUTDATADIR / f"gene_degrees_{group}_protein_coding.csv", index=False)
 
def distribucion_grado_plot():
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.flatten()

    bins = np.arange(0, 320, 5)

    colores = ["#1f77b4", "#f1c40f", "#e67e22", "#e21f0a"]

    for ax, group, color in zip(axes, GROUPS, colores):
        df = pd.read_csv(OUTDATADIR /  f"gene_degrees_{group}_protein_coding.csv")

        ax.hist(df["degree"], bins=bins, color=color, edgecolor="black")
        ax.set_title(group)
        ax.set_yscale("log")
        ax.set_xlabel("Degree")
        ax.set_ylabel("Number of genes")

    fig.suptitle("Degree distribution by AD severity group")
    fig.tight_layout()
    fig.savefig(OUTFIGDIR/"degree_distribution_4groups_vip_protein_coding_01.png", dpi=300)
    plt.close(fig)
        
def network_metrics(G):
    components = list(nx.connected_components(G))
    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "connected_components": len(components),
        "clustering_coefficient": nx.average_clustering(G),
    }
    
def jaccard_index(G1, G2):
    e1 = set(frozenset(e) for e in G1.edges())
    e2 = set(frozenset(e) for e in G2.edges())
    return len(e1 & e2) / len(e1 | e2)
        
def networks_properties(graphs):
    results = []

    for group in GROUPS:
        G = graphs[group]

        largest_cc = max(nx.connected_components(G), key=len)
        G_giant = G.subgraph(largest_cc)

        avg_path_length = nx.average_shortest_path_length(G_giant)

        fila = {
            "group": group,
            **network_metrics(G),
            "giant_component_nodes": G_giant.number_of_nodes(),
            "avg_shortest_path_length": avg_path_length,
        }
        results.append(fila)

    df_results = pd.DataFrame(results)
    df_results.to_csv(
        OUTDATADIR / "network_properties_protein_coding_vip01.csv",
        index=False,
    )
    
    jaccard_results= []

    for i, k in combinations(GROUPS, 2):
        jac = jaccard_index(graphs[i], graphs[k])
        jaccard_results.append({
            "group_1": i,
            "group_2": k,
            "jaccard_index": jac,
        })
        print(f"{i} vs {k}: jaccard = {jac:.4f}")

    df_jaccard = pd.DataFrame(jaccard_results)
    df_jaccard.to_csv(
        OUTDATADIR / "jaccard_index_protein_coding.csv",
        index=False,
    )
        

def plot_network_architecture(graphs):
    colors = {
        "Not_AD": "#8B5CF6",     # morado
        "Low": "#3B82F6",       # azul
        "Intermediate": "#10B981",  # verde
        "High": "#EF4444",      # rojo
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    axes = axes.flatten()

    for ax, group in zip(axes, GROUPS):
        G = graphs[group]  # grafo completo ya construido

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
    plt.savefig(OUTFIGDIR / "network_architecture.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Guardado en OUTFIGDIR / network_architecture.png")

def degree_changue():
    
    all_genes=set()
    for group in GROUPS:
        all_genes.update(pd.read_csv(OUTDATADIR /  f"gene_degrees_{group}_protein_coding.csv")["gene"])
    print(f"all genes len {len(all_genes)}")
    
    df_degree = pd.DataFrame({"gene": sorted(all_genes)})

    for group in GROUPS:
        df_group = pd.read_csv(OUTDATADIR / f"gene_degrees_{group}_protein_coding.csv")
        df_group = df_group.rename(columns={"degree": group}) 
        df_degree = df_degree.merge(df_group, on="gene", how="left")

    # genes eliminados por conteo 0 en un grupo -> grado 0 en ese grupo
    df_degree[GROUPS] = df_degree[GROUPS].fillna(0)

    df_degree["max_range"] = df_degree[GROUPS].max(axis=1) - df_degree[GROUPS].min(axis=1)
    df_degree = df_degree.sort_values("max_range", ascending=False).reset_index(drop=True)

    df_degree.to_csv(OUTDATADIR / "gene_degrees_all_groups_protein_coding.csv", index=False)
    
    
        

if __name__ == "__main__":
    graphs=read_filteres_graphs_pc()
    # degree_variation_analysis(cutoffs, graphs)
    #distribucion_grado_plot()
    #networks_properties(graphs)
    # plot_network_architecture(graphs)
    degree_changue()