import polars as pl
import pandas as pd
from pathlib import Path
from numba_mi import pipeline
from datetime import datetime

GROUPS=["Not_AD", "Low", "Intermediate", "High"]

INDIR= Path("data/count_matrices/all_vip")

OUTDIR_MATRICES= (INDIR / "protein_coding")
OUTDIR_MATRICES.mkdir(parents=True, exist_ok=True)

OUTDIR_MI = Path("output/protein_coding/mi_matrices")
OUTDIR_MI.mkdir(parents=True, exist_ok=True)

coding_genes= pd.read_csv("data/gene_lists/vip/annotations/gene_annotation_protein_coding.csv")

coding_names = pd.unique(coding_genes["gene_name"].dropna())
coding_ids = pd.unique(coding_genes["gene_id"].dropna())

coding_names_ids = set(coding_ids) | set(coding_names)

def write_p_coding_tsv():
    for group in GROUPS:
        count_matrix= pl.read_csv(INDIR / f"Vip_{group}.tsv", separator= "\t")
        first_col = count_matrix.columns[0]
        count_matrix= count_matrix.filter(pl.col(first_col).is_in(coding_names_ids))
        count_matrix.write_csv( OUTDIR_MATRICES / f"Vip_{group}_protein_coding.tsv", separator= "\t")
    
def calculo_mi_numba():
    for group in GROUPS:
        
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]")
        print(f"Calculando MI para grupo: {group}")
        
        file_path= OUTDIR_MATRICES / f"Vip_{group}_protein_coding.tsv"
        mi_matrix = pipeline.whole_process(
                file_path=file_path,
                axis=1,
                sep="\t",
                index_col=0,
                num_threads=40,
            )
        mi_matrix.write_parquet(OUTDIR_MI / f"Vip_{group}_protein_coding_mi.parquet")
        
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]")
        print(f"Guardado MI para grupo: {group} en {OUTDIR_MI / f'Vip_{group}_protein_coding_mi.parquet'}")
        
if __name__ == "__main__":
    write_p_coding_tsv()
    calculo_mi_numba()