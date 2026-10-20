import pandas as pd
import numpy as np
import polars as pl
import networkx as nx
from pathlib import Path
from scipy.stats import hypergeom

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
CONTRASTS= [("Not_AD", "Low"), ("Low", "Intermediate"), ("Intermediate", "High"), ("Not_AD", "High")]

MIN_MODULE_SIZE= 10

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = Path(f"{BASE}/output/protein_coding/mi_matrices")
COMUNIDADES= Path(f"{BASE}/output/protein_coding/data/comunidades_infomap")
OUTFIGDIR = Path(f"{BASE}/output/protein_coding/figures")
OUTDATADIR = Path(f"{BASE}/output/protein_coding/data")
OUTDATADIR_GENE= Path(f"{OUTDATADIR}/gene_analysis_changes")

OUTDATADIR_GENE.mkdir(parents=True, exist_ok=True)

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

def read_comunidades(group):
    df_comunidades = pd.read_csv(f"{COMUNIDADES}/infomap_comunidades_{group}.csv")

    num_comunidades= df_comunidades["modulo"].unique()

    comunidades = {i: df_comunidades.loc[df_comunidades["modulo"] ==i, "gen"].tolist()
                   for i in num_comunidades}

    return comunidades

def neighbors_by_state(gene, G):
    if G.has_node(gene):
        return set(G.neighbors(gene))
    return set()

def gene_neighbors_change(graphs, tag):

    results = []

    for group1, group2 in CONTRASTS:
        G1 = graphs[group1]
        G2 = graphs[group2]

        nodes1 = set(G1.nodes)
        nodes2 = set(G2.nodes)
        shared = nodes1 & nodes2
        n_shared = len(shared)

        genes = sorted(set(G1.nodes) | set(G2.nodes))
        n_universe = len(genes)

        for gene in genes:

            gene_set1 = neighbors_by_state(gene, G1)
            gene_set2 = neighbors_by_state(gene, G2)

            if G1.has_node(gene) and G2.has_node(gene):
                status = "in_both"
            elif G2.has_node(gene):
                status = "in"
            else:
                status = "out"

            conserved = gene_set1 & gene_set2
            gained = gene_set2 - gene_set1
            lost = gene_set1 - gene_set2
            union = gene_set1 | gene_set2

            s1 = gene_set1 & shared
            s2 = gene_set2 & shared

            result = {
                "gene": gene,
                "group1": group1,
                "group2": group2,
                "status": status,
                "degree1": len(gene_set1),
                "degree2": len(gene_set2),
                "n_universe": n_universe,
                "n_conserved": len(conserved),
                "n_gained": len(gained),
                "n_lost": len(lost),
                "net_change_degree": len(gained) - len(lost),
                "edges_changed": len(gained) + len(lost),
                "jaccard": (len(conserved) / len(union) if union else None),
                "n_shared": n_shared, 
                "degree1_shared": len(s1),
                "degree2_shared": len(s2),
                "conserved_neighbors": sorted(conserved),
                "gained_neighbors": sorted(gained),
                "lost_neighbors": sorted(lost)
            }

            results.append(result)

    df_neighbors_change = pd.DataFrame(results)

    df_neighbors_change.to_csv(
        OUTDATADIR_GENE / f"gene_neighbors_change_{tag}.csv",
        index=False
    )

    print(f"Gene neighbors change analysis completed. Results saved to {OUTDATADIR_GENE / f'gene_neighbors_change_{tag}.csv'}.")

    return df_neighbors_change

def genes_in_out(df_neighbors_change, tag):
    in_out = (
        df_neighbors_change
        .groupby(["group1", "group2", "status"])
        .agg(
            count=("gene", "size"),
            genes=("gene", list)
        )
        .reset_index()
    )

    in_out.to_csv(
        OUTDATADIR_GENE / f"in_out_genes_{tag}.csv",
        index=False
    )
    
    return in_out

def neighbor_conservation_test(df_neighbors_change, tag):
    df_in_both= df_neighbors_change[df_neighbors_change["status"]=="in_both"].copy()
    
    df_in_both["expected"] = (
        df_in_both["degree1_shared"] * df_in_both["degree2_shared"]
        / (df_in_both["n_shared"] - 1)
    )
    
    df_in_both["p_value"] = hypergeom.sf(
        df_in_both["n_conserved"] - 1,
        df_in_both["n_shared"] - 1,
        df_in_both["degree1_shared"],
        df_in_both["degree2_shared"]
    )
    
    df_in_both.to_csv(
        OUTDATADIR_GENE / f"neighbor_conservation_{tag}.csv",
        index=False
    )
    
    return df_in_both

def intracommunity_metrics(communities, G, min_size):
    
    communities_minsz = {
        community: genes
        for community, genes in communities.items()
        if len(genes) >= min_size
    }
    
    results=[]
    
    for com, genes in communities_minsz.items():
        
        degrees_in = [
            sum(neighbor in genes for neighbor in G.neighbors(gene))
            for gene in genes
        ]

        results.append({
            "community": com,
            "size": len(genes),
            "mean_degree_in": np.mean(degrees_in),
            "std_degree_in": np.std(degrees_in, ddof=0)
        })

    return pd.DataFrame(results)


def gene_roles(G, communities, group, min_size):
    intracommunitym = intracommunity_metrics(communities, G, min_size)

    minsz_communities = intracommunitym["community"].tolist()

    results = pd.DataFrame()

    genes = list(G.nodes())

    results["group"] = [group] * len(genes)
    results["gene"] = genes
    results["degree"] = [G.degree(gene) for gene in genes]

    for gene in genes:

        community = None

        for com, community_genes in communities.items():
            if gene in community_genes:
                community = com
                break

        if community is None:
            raise ValueError(f"Gene {gene} has no assigned community.")

        degree = G.degree(gene)

        degree_in = sum(
            neighbor in communities[community]
            for neighbor in G.neighbors(gene)
        )

        z = np.nan

        if community in minsz_communities:

            mean_degree_in = intracommunitym.loc[
                intracommunitym["community"] == community, "mean_degree_in"
            ].item()

            std_degree_in = intracommunitym.loc[
                intracommunitym["community"] == community, "std_degree_in"
            ].item()

            if std_degree_in > 0:
                z = (degree_in - mean_degree_in) / std_degree_in

        P = np.nan

        if degree > 0:
            squared_sum = 0

            for com, community_genes in communities.items():
                neighbors_in_com = sum(
                    neighbor in community_genes
                    for neighbor in G.neighbors(gene)
                )

                squared_sum += (neighbors_in_com / degree) ** 2

            P = 1 - squared_sum

        row = results["gene"] == gene

        results.loc[row, "community"] = community
        results.loc[row, "community_size"] = len(communities[community])
        results.loc[row, "degree_in"] = degree_in
        results.loc[row, "z"] = z
        results.loc[row, "P"] = P

    return results

def roles_all_groups(graphs, min_size, tag):
    
    all_roles= []
    
    for group in GROUPS:
        
        communities= read_comunidades(group)
        
        df_roles= gene_roles(
            G =graphs[group],
            communities=communities,
            group= group,
            min_size=min_size
        )
    
        all_roles.append(df_roles)
        
    df_roles_all_groups = pd.concat(all_roles, ignore_index=True)

    df_roles_all_groups.to_csv(
        OUTDATADIR_GENE / f"gene_roles_{tag}.csv",
        index=False
    )

    return df_roles_all_groups

def merge_gene_metrics(df_neighbors_change, df_tested, df_roles, tag):

    df_in_out = df_neighbors_change[df_neighbors_change["status"] != "in_both"]

    df = pd.concat([df_tested, df_in_out], ignore_index=True)

    roles1 = df_roles[["group", "gene", "z", "P"]].rename(columns={"group": "group1", "z": "z1", "P": "P1"})
    roles2 = df_roles[["group", "gene", "z", "P"]].rename(columns={"group": "group2", "z": "z2", "P": "P2"})

    df = (
        df
        .merge(roles1, on=["gene", "group1"], how="left")
        .merge(roles2, on=["gene", "group2"], how="left")
    )

    cols = [
        "gene", "group1", "group2", "status",
        "degree1", "degree2", "z1", "P1", "z2", "P2",
        "n_shared", "degree1_shared", "degree2_shared",
        "n_conserved", "expected", "p_value",
        "n_gained", "n_lost", "jaccard",
        "conserved_neighbors", "gained_neighbors", "lost_neighbors"
    ]

    df = df[cols]

    df.to_csv(
        OUTDATADIR_GENE / f"gene_metrics_{tag}.csv",
        index=False
    )

    return df


if __name__ == "__main__":

    tag= "protein_coding"

    graphs= read_filteres_graphs_pc()

    df_gene_neighbors= gene_neighbors_change(graphs, tag)

    genes_in_out(df_gene_neighbors, tag)

    df_tested= neighbor_conservation_test(df_gene_neighbors, tag)

    df_roles= roles_all_groups(graphs, MIN_MODULE_SIZE, tag)

    merge_gene_metrics(df_gene_neighbors, df_tested, df_roles, tag=tag)
        