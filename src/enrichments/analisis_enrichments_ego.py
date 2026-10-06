import pandas as pd
from pathlib import Path
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import polars as pl

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR_EGO = BASE / "output" / "red_mi_vip_completo" / "ego_groups"

GROUPS=["Not_AD", "Low", "Intermediate", "High"]
GROUPS_x_notad=["Low", "Intermediate", "High"]

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

def jaccard_ids(id1, id2):
    union= id1 | id1
    if not union:
        return 0
    return len(id1 & id2) / len(union)

def jaccard_BP_enrich(nodes):
    
    df_results=[]
    
    for group, jacc in combinations(GROUPS, 2):
        
            df_group=pd.read_csv(INDIR_EGO / "enrichment_GO_BP" / f"GO_BP_ego01_{nodes}_{group}.csv")
            df_jacc=pd.read_csv(INDIR_EGO / "enrichment_GO_BP" / f"GO_BP_ego01_{nodes}_{jacc}.csv")
            
            ids_group = set(df_group["ID"].dropna())
            ids_jacc = set(df_jacc["ID"].dropna())
            
            jaccard = jaccard_ids(ids_jacc, ids_group)
            
            df_results.append({
                "group_1": group,
                "group_2": jacc,
                "jaccard": jaccard
            })
    
    df_results = pd.DataFrame(df_results)
    
    print(df_results)
    
    return df_results

def plot_jaccard_BP(df_results, nodes):

    matrix = df_results.pivot(
        index="group_1",
        columns="group_2",
        values="jaccard"
    )

    matrix = matrix.combine_first(matrix.T)
    matrix = matrix.reindex(index=GROUPS, columns=GROUPS)

    mask = np.triu(np.ones_like(matrix, dtype=bool), k=0)

    plt.figure(figsize=(6, 5))

    sns.heatmap(
        matrix,
        mask=mask,
        annot=True,
        fmt=".3f",
        square=True,
        cmap="viridis",
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Jaccard coefficient"}
    )

    plt.title("Jaccard index of GO BP terms between disease stages")
    plt.tight_layout()
    plt.savefig(
        INDIR_EGO / "figures" / f"jaccard_BP_heatmap_{nodes}.png",
        dpi=300
    )
    plt.show()
    plt.close()

def plot_hits_enrichment_allnet(nodes):

    graphs=read_graphs(nodes)
        
    colors = {
            "No_hit": "#005DF3",      
            "Hit": "#EA1D1D",  
    }  

    graph_ids = {}

    for group in GROUPS:
                df_group=pd.read_csv(INDIR_EGO / "enrichment_GO_BP" / f"GO_BP_ego01_{nodes}_{group}.csv")
                graph_ids[group] = set(df_group["ID"].dropna())

    common_ids = set.intersection(*graph_ids.values())

    print(common_ids)

    hits_by_state= {}


    for group in GROUPS:
        df_group=pd.read_csv(INDIR_EGO / "enrichment_GO_BP" / f"GO_BP_ego01_{nodes}_{group}.csv")
        df_group=df_group[df_group["ID"]=='GO:2000779']
        genes = df_group["geneID"].iloc[0].split("/")
        hits_by_state[group]= genes
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    axes = axes.flatten()
    
    
    for ax, group in zip(axes, GROUPS):
        G = graphs[group] 
        
        pos = nx.spring_layout(G, k=0.5, iterations=100, seed=42)


        degrees = dict(G.degree())
        node_sizes = [20 + degrees[n] * 0.0001 for n in G.nodes()]

        nx.draw_networkx_nodes(
            G, pos, ax=ax,
            node_size=node_sizes,
            node_color=[
            colors["Hit"] if n in hits_by_state[group]
            else colors["No_hit"]
            for n in G.nodes()
            ],
            alpha=0.85,
            linewidths=0.2,
            edgecolors="black",
        )

        nx.draw_networkx_edges(
            G, pos, ax=ax,
            alpha=0.10,
            width=0.8,
            edge_color="black",
        )

        ax.set_title(
            f"{group}\n(hits del termino 'GO:2000779')",
            fontsize=14, fontweight="bold"
        )
        ax.axis("off")

    plt.suptitle(f"'GO:2000779' hits across AD severity groups {nodes}", fontsize=18, y=0.98)
    plt.tight_layout()
    plt.savefig(INDIR_EGO / "figures" / f"enrichment_hits_GO:2000779_{nodes}.png", dpi=300, bbox_inches="tight")
    plt.close()

def plot_hits_new_genes(nodes):
    graphs=read_graphs(nodes)
     
    new_genes_dict = {}
         
    new_genes= pd.read_csv(INDIR_EGO / f"new_not_reappeared_genes_{nodes}.csv")
             
     
    for _, row in new_genes.iterrows():
        genes = row["new_genes"].strip("[]").replace("'", "").split(", ")
        new_genes_dict[row["group"]] = genes   
     
    subgraphs = {group: graphs[group].subgraph(new_genes_dict[group]) for group in GROUPS_x_notad}
        
    colors = {  
            "No_hit": "#005DF3",      
            "Hit": "#EA1D1D",  
    }  

    hits_by_state = {}
    id_by_state = {}

    for group in GROUPS_x_notad:
        df_group = pd.read_csv(
            INDIR_EGO / "enrichment_GO_BP_new_genes" / f"GO_BP_new_genes_{nodes}_{group}.csv"
        )

        first_id = df_group["ID"].iloc[0]

        df_group = df_group[df_group["ID"] == first_id]

        genes = df_group["geneID"].iloc[0].split("/")
        hits_by_state[group] = genes

        id_by_state[group] = {
            "ID": first_id,
            "Description": df_group["Description"].iloc[0]
        }
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    axes = axes.flatten()

    for ax, group in zip(axes, GROUPS_x_notad):
        G = subgraphs[group]

        pos = nx.spring_layout(
            G,
            k=0.5,
            iterations=100,
            seed=42
        )

        degrees = dict(G.degree())
        node_sizes = [20 + degrees[n] * 0.6 for n in G.nodes()]

        nx.draw_networkx_nodes(
            G,
            pos,
            ax=ax,
            node_size=node_sizes,
            node_color=[
                colors["Hit"] if n in hits_by_state[group]
                else colors["No_hit"]
                for n in G.nodes()
            ],
            alpha=0.85,
            linewidths=0.2,
            edgecolors="black",
        )

        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            alpha=0.10,
            width=0.8,
            edge_color="black",
        )

        ax.set_title(
            f"{group}\n"
            f"{id_by_state[group]['ID']}\n"
            f"{id_by_state[group]['Description']}",
            fontsize=14,
            fontweight="bold"
        )

        ax.axis("off")

    plt.suptitle(
        f"Subgraphs of Top GO BP enrichment hits in new genes across AD severity groups ({nodes})",
        fontsize=18,
        y=0.98
    )

    plt.tight_layout()

    plt.savefig(
        INDIR_EGO / "figures" / f"enrichment_hits_topGO_{nodes}_new_genes_subgraphs.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
   
if __name__ == "__main__":
    nodes="SR_SUG_LUC"
    plot_hits_enrichment_allnet(nodes)