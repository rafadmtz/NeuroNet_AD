import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import numpy as np
import gseapy as gp
import infomap
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
INDIR = Path("/export/space3/users/silvanac/NeuroNet_AD/output/red_mi_vip_completo")
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR_EGO = BASE / "output" / "red_mi_vip_completo" / "ego_groups"

def read_graphs(nodes):
    edgelists = {
            group: pl.read_parquet(
                INDIR_EGO / f"ego01_{nodes}_{group}.parquet"
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

def jaccard_nodes(G1, G2):
    n1, n2 = set(G1.nodes()), set(G2.nodes())
    return len(n1 & n2) / len(n1 | n2)

def jaccard_edges(G1, G2):
    e1 = set(frozenset(e) for e in G1.edges())
    e2 = set(frozenset(e) for e in G2.edges())
    return len(e1 & e2) / len(e1 | e2)

def enrichment_neighrborhood(nodes, communities=False):
    graphs = read_graphs(nodes) 
    df_results = []
    
    for group in GROUPS:
        
        G = graphs[group]
        
        try:
            enr = gp.enrichr(
                gene_list = list(G.nodes()),
                gene_sets = "GO_Biological_Process_2023",
                organism = "human",
                cutoff = 0.05,
            )
            
            df_group = enr.results.assign(
                Comunidad=group, 
                tamaño_comunidad=len(list(G.nodes()))
            )
            
            df_results.append(df_group)
            
            outdir = INDIR_EGO / f"enrichment_ego01_{nodes}"
            outdir.mkdir(parents=True, exist_ok=True)
            
            df_group = df_group.sort_values("Adjusted P-value")
            df_group.to_csv(
                outdir / f"enrichment_{nodes}_{group}_01.csv",
                index=False
            )
            
            print(f"resultados_enrichment_{nodes}_{group}_01.csv guardados")
            
        except Exception as e:
            print(
                f"[{group}] comunidad {nodes} "
                f"(n={len(list(G.nodes()))} genes) falló: {e}"
            )


def subgraph_neighbors(nodes, group):
    graphs = read_graphs(nodes)
    G = graphs[group]

    nodes_to_keep = set()

    for node_interest in nodes:
        if node_interest in G:
            neighbors = list(G.neighbors(node_interest))
            nodes_to_keep.update(neighbors)
            nodes_to_keep.add(node_interest)
        else:
            print(f"{node_interest} no se encuentra en el grupo {group}.")

    if not nodes_to_keep:
        return None

    subgraph = G.subgraph(nodes_to_keep)
    return subgraph

def infomap_communities(nodes,two_level= True):
    for group in GROUPS:
        df = pl.read_parquet(
            INDIR_EGO / f"ego01_{nodes}_{group}.parquet"
            )

        G_subgraph = nx.from_pandas_edgelist(
            df.to_pandas(),
            source="source",
            target="target",
            edge_attr="MI"
            )
        
        comm=infomap.run(G_subgraph, 
                        seed= 2,
                        num_trials=1000,
                        two_level=two_level,
                        directed=False
                        )

        df_comm=comm.to_dataframe()

        df_comm = df_comm.rename(columns={"module_id": "community", "name": "gene"})[["gene", "community"]]
        
        # Tamaño de cada comunidad
        community_sizes = (
            df_comm
            .groupby("community")
            .size()
            .sort_values(ascending=False)
        )

        print(group)
        print()
        
        for community, size in community_sizes.items():
            print(f"Comunidad {community}: {size} genes")
        
        print()
        
        outdir = INDIR_EGO / f"infomap_communities_ego01_{nodes}"
        outdir.mkdir(parents=True, exist_ok=True)

        # Guardarlo
        df_comm.to_csv(outdir /f"comunidades_ego01_{nodes}_{group}_infomap.csv", index=False)
        
def jaccard_analysis(nodes):
    
    edgelists = {
        group: pl.read_parquet(
            INDIR_EGO / f"ego01_{nodes}_{group}.parquet"
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
    df_jaccard.to_csv(INDIR_EGO / f"jaccard_analysis_{nodes}.csv", index=False)
    return df_jaccard

def heatmap_jaccard():
    
    df_jaccard = pd.read_csv(INDIR_EGO / f"jaccard_analysis_{nodes}.csv")
    
    df_nodes = df_jaccard.pivot(index="group_1", columns="group_2", values="jaccard_nodes")
    df_edges = df_jaccard.pivot(index="group_1", columns="group_2", values="jaccard_edges")
    
    groups = ["Not_AD", "Low", "Intermediate", "High"]
    
    df_nodes = df_nodes.reindex(index=groups, columns=groups)
    df_edges = df_edges.reindex(index=groups, columns=groups)
    
    df_nodes = df_nodes.combine_first(df_nodes.T)
    df_edges = df_edges.combine_first(df_edges.T)
    
    df_nodes.index.name = None
    df_nodes.columns.name = None

    df_edges.index.name = None
    df_edges.columns.name = None
    
    mask = np.triu(np.ones_like(df_nodes, dtype=bool))
    
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        df_nodes,
        mask=mask,
        annot=True,
        fmt=".3f",
        square=True,
        cmap="viridis",
        vmin=0,
        vmax=0.4,
        cbar_kws={"label": "Jaccard index"}
    )
    
    plt.title("Jaccard index between disease stages")
    plt.tight_layout()
    plt.savefig(INDIR_EGO / f"jaccard_nodes_heatmap_{nodes}.png", dpi=300)
    plt.show()
    plt.close()
   
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        df_edges,
        mask=mask,
        annot=True,
        fmt=".3f",
        square=True,
        cmap="viridis",
        vmin=0,
        vmax=0.08,
        cbar_kws={"label": "Jaccard index"}
    )
    
    plt.title("Jaccard index between disease stages")
    plt.tight_layout()
    plt.savefig(INDIR_EGO / f"jaccard_edges_heatmap_{nodes}.png", dpi=300)
    plt.show()
    plt.close()
    
    

def overlap_nodes(G1, G2):
    
    n1, n2 = set(G1.nodes()), set(G2.nodes())
    
    return len(n1 & n2) / min(len(n1), len(n2))


def overlap_edges(G1, G2):
    
    e1 = set(frozenset(e) for e in G1.edges())
    e2 = set(frozenset(e) for e in G2.edges())
    
    return len(e1 & e2) / min(len(e1), len(e2))


def overlap_analysis(nodes):
    
    edgelists = {
        group: pl.read_parquet(
            INDIR_EGO / f"ego01_{nodes}_{group}.parquet"
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
        for ov in GROUPS:
            if ov == group:
                continue
            pair = frozenset({group, ov})
            if pair in seen:
                continue
            seen.add(pair)

            results.append({
                "group_1": group,
                "group_2": ov,
                "overlap_nodes": overlap_nodes(graphs[group], graphs[ov]),
                "overlap_edges": overlap_edges(graphs[group], graphs[ov]),
            })

    df_overlap = pd.DataFrame(results)
    df_overlap.to_csv(
        INDIR_EGO / f"overlap_analysis_{nodes}.csv",
        index=False
    )
    
    df_nodes = df_overlap.pivot(
        index="group_1",
        columns="group_2",
        values="overlap_nodes"
    )
    
    df_edges = df_overlap.pivot(
        index="group_1",
        columns="group_2",
        values="overlap_edges"
    )
    
    groups = ["Not_AD", "Low", "Intermediate", "High"]
    
    df_nodes = df_nodes.reindex(index=groups, columns=groups)
    df_edges = df_edges.reindex(index=groups, columns=groups)
    
    df_nodes = df_nodes.combine_first(df_nodes.T)
    df_edges = df_edges.combine_first(df_edges.T)
    
    df_nodes.index.name = None
    df_nodes.columns.name = None

    df_edges.index.name = None
    df_edges.columns.name = None
    
    mask = np.triu(np.ones_like(df_nodes, dtype=bool))
    
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        df_nodes,
        mask=mask,
        annot=True,
        fmt=".3f",
        square=True,
        cmap="viridis",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Overlap coefficient"}
    )
    
    plt.title("Node overlap between disease stages")
    plt.tight_layout()
    plt.savefig(
        INDIR_EGO / f"overlap_nodes_heatmap_{nodes}.png",
        dpi=300
    )
    plt.show()
    plt.close()
    
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        df_edges,
        mask=mask,
        annot=True,
        fmt=".3f",
        square=True,
        cmap="viridis",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Overlap coefficient"}
    )
    
    plt.title("Edge overlap between disease stages")
    plt.tight_layout()
    plt.savefig(
        INDIR_EGO / f"overlap_edges_heatmap_{nodes}.png",
        dpi=300
    )
    plt.show()
    plt.close()
    
    return df_overlap

def basic_network_metrics(nodes):
    
    graphs = read_graphs(nodes)
    
    df_results = []
    
    for group, G in graphs.items():
        print(f"Group: {group}")
        print(f"Number of nodes: {G.number_of_nodes()}")
        print(f"Number of edges: {G.number_of_edges()}")
        
        G_core = nx.k_core(G)
        G_giant = max(nx.connected_components(G), key=len)
        G_giant = G.subgraph(G_giant).copy()
        G_connnected_components = len(list(nx.connected_components(G)))
        
        print(f"K-core number of nodes: {G_core.number_of_nodes()}")
        print(f"K-core number of edges: {G_core.number_of_edges()}")
        
        print(f"Giant component number of nodes: {G_giant.number_of_nodes()}")
        print(f"Giant component number of edges: {G_giant.number_of_edges()}")
        
        clustering_coefficient = nx.average_clustering(G_core)
        avg_path_length = nx.average_shortest_path_length(G_giant)
        
        print(f"Average clustering coefficient (k-core): {clustering_coefficient:.4f}")
        print(f"Average shortest path length (giant component): {avg_path_length:.4f}")
        print("-" * 40)
    
        df_results.append({
            "group": group,
            "nodes": G.number_of_nodes(),
            "edges": G.number_of_edges(),
            "k_core_nodes": G_core.number_of_nodes(),
            "k_core_edges": G_core.number_of_edges(),
            "giant_component_nodes": G_giant.number_of_nodes(),
            "giant_component_edges": G_giant.number_of_edges(),
            "avg_clustering_coefficient": clustering_coefficient,
            "avg_shortest_path_length": avg_path_length,
            "connected_components": G_connnected_components
        })
        
    df_results = pd.DataFrame(df_results)
    df_results.to_csv(
        INDIR_EGO / f"basic_network_metrics_{nodes}.csv",
        index=False
    )
        
from venn import venn

def plot_venn_nodes(nodes):
    
    graphs = read_graphs(nodes)
    node_sets = {
        group: set(graphs[group].nodes())
        for group in GROUPS
    }

    fig, ax = plt.subplots(figsize=(10, 10))
    venn(node_sets, ax=ax)
    plt.title(f"Node composition overlap across groups in {nodes}")
    plt.savefig(INDIR_EGO / f"venn_{nodes}_nodes.png", dpi=300, bbox_inches="tight")
    plt.close()
    
if __name__ == "__main__":
    
    
    
    construir_subgrafo=0
    if construir_subgrafo:
        nodes = ["SRCIN1", "SUGP2", "LUC7L"]
        groups = ["Not_AD", "Low", "Intermediate", "High"]

        OUTPUT_DIR = BASE / "output" / "red_mi_vip_completo" / "ego_groups"

        for group in groups:
            sg = subgraph_neighbors(nodes, group)
            if sg is None:
                print(f"Sin nodos de interes encontrados en {group}, se omite.")
                continue

            print(f"{group}: {sg.number_of_nodes()} nodos, {sg.number_of_edges()} aristas")

            edgelist = nx.to_pandas_edgelist(sg, source="source", target="target")
            df = pl.from_pandas(edgelist)
            df.write_parquet(f"{OUTPUT_DIR}/ego01_SR_SUG_LUC_{group}.parquet")
            
    nodes="SR_SUG_LUC"
    
    infomap_com=0
    if infomap_com:
        infomap_communities(nodes)
        
    enrich=1
    if enrich:
        enrichment_neighrborhood(nodes)
        
    #jaccard_analysis(nodes)
    
    jaccard_heatmap=0
    if jaccard_heatmap:
        heatmap_jaccard()
    
    #overlap_analysis(nodes)
    # plot_venn_nodes(nodes)
    # basic_network_metrics(nodes)