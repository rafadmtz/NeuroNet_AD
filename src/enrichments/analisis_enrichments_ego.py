import pandas as pd
from pathlib import Path
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR_EGO = BASE / "output" / "red_mi_vip_completo" / "ego_groups"

GROUPS=["Not_AD", "Low", "Intermediate", "High"]


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

if __name__ == "__main__":
    nodes="SR_SUG_LUC"
    jaccard=jaccard_BP_enrich(nodes)
    plot_jaccard_BP(jaccard, nodes)