import scanpy as sc
import anndata as ad
import h5py
from anndata._io.specs import read_elem
import pandas as pd
import numpy as np
import sys
import os
import pandas as pd
from joblib import Parallel, delayed
import glob

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = PROJECT_ROOT / "whole_taxonomy_MTG_AD.h5ad"
OUTPUT_DIR = PROJECT_ROOT / "output"

def load_scanpy(path: str) -> ad.AnnData:
    """
    Load scanpy.
    
    Returns:
        ad.AnnData: An AnnData object containing the test dataset.
    """
    # Load the example dataset from Scanpy
    adata = sc.read_h5ad(path, backed='r')
    # Return the AnnData object
    return adata


def load_metadata(path: str) -> pd.DataFrame:
    """
    Read and load only metadata  (obs).
    
    Returns:
        pd.DataFrame
    """
    with h5py.File(path) as f:
        obs = read_elem(f["obs"])
    return obs

def load_gene(path:str) -> pd.DataFrame:
    """
    Read and load only the genes.
    
    Returns:
        pd.DataFrame
    """
    with h5py.File(path) as f:
        genes = read_elem(f["var"])
    return genes


def save_matrix_patient_type(row, path, output_dir) -> None:
    """
    Read and load matrix for each patient and subclass. (84*24*2)
    Args:
        row: A row from the DataFrame containing "Donor ID" and "Subclass".
        adata: The AnnData object containing the dataset.

    Returns:
        pd.DataFrame
    """
    adata = load_scanpy(path)
    unique_adata = adata[(adata.obs["donor_id"] == row["donor_id"]) & (adata.obs["Subclass"] == row["Subclass"]), :]
    unique_adata = unique_adata.to_memory()
    unique_adata.var_names = unique_adata.var["feature_name"].values
    unique_adata.var_names_make_unique()                               
    df = unique_adata.to_df().T
    donor = row["donor_id"].replace(".", "-")
    subclass = row["Subclass"].replace(" ", "-").replace("/", "-")
    df.to_csv(output_dir / f"{donor}_{subclass}.tsv", sep="\t")
    print(f"Saved matrix for patient {donor} and subclass {subclass}.")
    


def main():
    n_threads = 6
    print("Starting Python prediction pipeline...")
    
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"Input file: {INPUT_FILE}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    ###
    #Load metadata and save csv
    ###
    load_metadata(INPUT_FILE).to_csv(OUTPUT_DIR / "metadata_test.csv", index=False)


    ###
    #Load genes and save csv
    load_gene(INPUT_FILE).to_csv(OUTPUT_DIR / "genes.txt", index=False, header=False)

    ###
    #Save matrix for each patient and subclass for MTG
    ###
    df_mtg = pd.read_csv(OUTPUT_DIR / "patients_subclass.csv")
    Parallel(n_jobs=n_threads)(
        delayed(save_matrix_patient_type)(row, INPUT_FILE, OUTPUT_DIR)
        for _, row in df_mtg.iterrows()
    )

if __name__ == "__main__":
    main()
