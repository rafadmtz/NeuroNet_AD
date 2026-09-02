import pandas as pd
import numpy as np
from pathlib import Path
import networkx as nx
import polars as pl
import plotly.graph_objects as go

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR_EGO = BASE / "output" / "red_mi_vip_completo" / "ego_groups"


GROUPS=["Not_AD", "Low", "Intermediate", "High"]

def read_graphs(nodes):
    edgelists = {
        group: pl.read_parquet(
            INDIR_EGO / f"ego01_{nodes}_{group}.parquet"
        )
        for group in GROUPS
    }

    graphs = {
        group: nx.from_pandas_edgelist(
            df.to_pandas(),
            source="source",
            target="target",
            edge_attr="MI"
        )
        for group, df in edgelists.items()
    }

    return graphs

def genes_incorporados(nodes):

    df_results = []

    graphs = read_graphs(nodes)

    # Sets de genes para facilitar las operaciones
    gene_sets = {
        group: set(graphs[group].nodes())
        for group in GROUPS
    }

    # Solo comparar estadios consecutivos
    for i in range(len(GROUPS) - 1):

        group = GROUPS[i]
        incorp = GROUPS[i + 1]

        genes_group = gene_sets[group]
        genes_incorp = gene_sets[incorp]

        # Nuevos respecto al estadio anterior
        incorporated = genes_incorp - genes_group

        # Se pierden en el siguiente estadio
        lost = genes_group - genes_incorp

        # Compartidos entre ambos estadios
        retained = genes_group & genes_incorp

        # Genes que permanecen desde este estadio hasta High
        future_groups = GROUPS[i:]

        persistent = set.intersection(
            *(gene_sets[g] for g in future_groups)
        )

        # Retenidos entre ambos estadios,
        # pero que no permanecen continuamente hasta High
        retained_nonpersistent = retained - persistent

        df_results.append({
            "group": group,
            "incorporado": incorp,

            # Totales
            "total_group": len(genes_group),
            "total_incorporado": len(genes_incorp),

            # Nuevos
            "incorporated_genes": list(incorporated),
            "num_incorporated": len(incorporated),

            # Perdidos
            "lost_genes": list(lost),
            "num_lost": len(lost),

            # Todos los retenidos
            "retained_genes": list(retained),
            "num_retained": len(retained),

            # Persistentes hasta High
            "persistent_genes": list(persistent),
            "num_persistent": len(persistent),

            # Retenidos pero no persistentes
            "retained_nonpersistent_genes": list(
                retained_nonpersistent
            ),
            "num_retained_nonpersistent": len(
                retained_nonpersistent
            )
        })

    df_results = pd.DataFrame(df_results)

    df_results.to_csv(
        INDIR_EGO / f"genes_incorporados_{nodes}.csv",
        index=False
    )
    
    return pl.DataFrame(df_results)

def new_not_reappeared_genes(nodes):
    df_results = []
    graphs = read_graphs(nodes)

    old_genes = set(graphs[GROUPS[0]].nodes())

    for next in GROUPS[1:]:
        next_genes = set(graphs[next].nodes())
        new_genes = next_genes - old_genes
        old_genes = old_genes | next_genes

        df_results.append({
            "group": next,
            "new_genes": list(new_genes),
        })

    df_results = pd.DataFrame(df_results)
    df_results.to_csv(
        INDIR_EGO / f"new_not_reappeared_genes_{nodes}.csv",
        index=False
    )
    
    return df_results



def plot_genes_incorporados(df_results, nodes):

    labels = GROUPS.copy()

    sources = []
    targets = []
    values = []
    link_colors = []
    link_labels = []

    # Crear nodos New / Lost
    for row in df_results.iter_rows(named=True):

        group = row["group"]
        incorp = row["incorporado"]

        labels.append(f"New in {incorp}")
        labels.append(f"Lost after {group}")

    # Eliminar duplicados
    labels = list(dict.fromkeys(labels))

    # Índices para Plotly
    group_idx = {
        label: i
        for i, label in enumerate(labels)
    }

    # Colores de nodos
    def get_node_color(label):

        stage_colors = {
            "Not_AD": "#4C78A8",
            "Low": "#F28E2B",
            "Intermediate": "#2A9D8F",
            "High": "#A06CD5"
        }

        if label in stage_colors:
            return stage_colors[label]

        elif label.startswith("New"):
            return "#70A95B"

        elif label.startswith("Lost"):
            return "#D95F59"

        return "#BAB0AC"

    node_colors = [
        get_node_color(label)
        for label in labels
    ]

    # Crear conexiones
    for row in df_results.iter_rows(named=True):

        group = row["group"]
        incorp = row["incorporado"]

        new_label = f"New in {incorp}"
        lost_label = f"Lost after {group}"

        # Persistent
        if row["num_persistent"] > 0:

            sources.append(group_idx[group])
            targets.append(group_idx[incorp])
            values.append(row["num_persistent"])

            link_colors.append(
                "rgba(31, 78, 121, 0.75)"
            )

            link_labels.append(
                f"Persistent {group} → High"
            )

        # Retained pero no persistent
        if row["num_retained_nonpersistent"] > 0:

            sources.append(group_idx[group])
            targets.append(group_idx[incorp])
            values.append(
                row["num_retained_nonpersistent"]
            )

            link_colors.append(
                "rgba(148, 0, 211, 0.85)"
            )

            link_labels.append(
                f"Retained {group} → {incorp}"
            )

        # Nuevos
        if row["num_incorporated"] > 0:

            sources.append(group_idx[new_label])
            targets.append(group_idx[incorp])
            values.append(row["num_incorporated"])

            link_colors.append(
                "rgba(112, 169, 91, 0.45)"
            )

            link_labels.append(
                f"New in {incorp}"
            )

        # Perdidos
        if row["num_lost"] > 0:

            sources.append(group_idx[group])
            targets.append(group_idx[lost_label])
            values.append(row["num_lost"])

            link_colors.append(
                "rgba(217, 95, 89, 0.45)"
            )

            link_labels.append(
                f"Lost after {group}"
            )

    # Sankey
    fig = go.Figure(
        go.Sankey(

            node=dict(
                label=labels,
                color=node_colors,
                pad=20,
                thickness=20
            ),

            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=link_colors,
                label=link_labels,

                hovertemplate=(
                    "%{label}<br>"
                    "<b>%{value} genes</b>"
                    "<extra></extra>"
                )
            )
        )
    )

    fig.update_layout(
        title=f"Gene turnover across AD stages {nodes}",
        font_size=12,
        width=1050,
        height=450,
        autosize=False
    )

    fig.show(
        config={
            "responsive": False,
            "toImageButtonOptions": {
                "format": "png",
                "filename": f"gene_turnover_{nodes}",
                "width": 1050,
                "height": 450,
                "scale": 3
            }
        }
    )


if __name__ == "__main__":
    nodes="SR_SUG_LUC"
    new_not_reappeared_genes(nodes)
    

            
            
            
            