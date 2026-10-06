import pandas as pd
from pathlib import Path
import polars as pl
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = Path(f"{BASE}/output/protein_coding/mi_matrices")
OUTDIR = Path(f"{BASE}/output/protein_coding/mi_matrices/filtered_edgelists")

OUTDIR.mkdir(parents=True, exist_ok=True)

cutoffs = {
    "Not_AD": 0.2086,
    "Low": 0.2161,
    "Intermediate": 0.2401,
    "High": 0.2772
}


def filter_edges_by_cutoffs(cutoffs):
    edgelists = {
                group: pl.read_parquet(
                    INDIR / f"Vip_{group}_protein_coding_mi_edgelist.parquet"
                )
                for group in GROUPS
                }    
    
    filtered_edgelists={}
    
    for group in GROUPS:
        filtered_edgelists[group]= edgelists[group].filter(pl.col("weight") > cutoffs[group])
        filtered_edgelists[group].write_parquet(OUTDIR / f"Vip_{group}_protein_coding_mi_edgelist_filtered01.parquet")
        print(f"Filtered edgelist for {group} saved to {OUTDIR / f'Vip_{group}_protein_coding_mi_edgelist_filtered01.parquet'}")
        
if __name__ == "__main__":
    filter_edges_by_cutoffs(cutoffs)
        
    