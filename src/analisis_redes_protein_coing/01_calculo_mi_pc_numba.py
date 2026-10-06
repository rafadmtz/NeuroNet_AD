import os
import time
from pathlib import Path
import polars as pl
from numba_mi import pipeline
import pandas as pd

BASE = "/home/rafadiaz/export/NeuroNet_AD"
INDIR = f"{BASE}/data/count_matrices/all_vip"
INDIR_PC = f"{INDIR}/count_matrices_protein_coding"
OUTDIR = f"{BASE}/output/protein_coding/mi_matrices"
GROUPS = ["Not_AD", "Low", "Intermediate", "High"]
GENES = Path("/export/space3/users/rafadiaz/NeuroNet_AD/data/gene_lists/vip/gene_annotation_from_anndata.csv")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(INDIR_PC, exist_ok=True)

df = pd.read_csv(GENES)
df = df[df["gene_biotype"] == "protein_coding"]
protein_coding = df["gene_name"].to_list()


def filter_protein_coding_for_group(group):
    count_matrix = f"{INDIR}/Vip_{group}.tsv"

    gene_count = pd.read_csv(count_matrix, sep="\t", index_col=0)

    genes = [g for g in protein_coding if g in gene_count.index]

    gene_count = gene_count.loc[genes]
    gene_count.index.name = "gene"

    print(f"Tamaño de {group}:")
    print(gene_count.shape)

    output_file = f"{INDIR_PC}/Vip_{group}_protein_coding.tsv"

    gene_count.to_csv(output_file, sep="\t", index_label="gene")

    print(f"Guardado en: {output_file}")


def matriz_a_edgelist(ruta):
    adj = pl.read_parquet(ruta)
    col_genes = adj.columns[0]

    edgelist = adj.unpivot(index=col_genes, variable_name="target", value_name="weight")
    edgelist = edgelist.rename({col_genes: "source"})
    edgelist = edgelist.filter(pl.col("source") < pl.col("target"))

    # red de seguridad: ningún pseudo-gen debe entrar a la red
    basura = ["Unnamed: 0", "gene", ""]
    antes = len(edgelist)
    edgelist = edgelist.filter(
        ~pl.col("source").is_in(basura) & ~pl.col("target").is_in(basura)
    )
    if antes != len(edgelist):
        print(f"  ATENCIÓN: se eliminaron {antes - len(edgelist)} aristas de pseudo-genes")

    ruta_edgelist = ruta.replace(".parquet", "_edgelist.parquet")
    edgelist.write_parquet(ruta_edgelist, compression="zstd")

    print(f"Edgelist guardado en: {ruta_edgelist}")
    return edgelist


def calculate_mi_protein_coding():
    for group in GROUPS:

        filter_protein_coding_for_group(group)

        file_path = f"{INDIR_PC}/Vip_{group}_protein_coding.tsv"

        output_path = f"{OUTDIR}/Vip_{group}_protein_coding_mi.parquet"

        print(f"{group}: calculando MI...")
        start = time.time()
        mi_matrix = pipeline.whole_process(
            file_path=file_path,
            axis=1,
            sep="\t",
            index_col=0,
            num_threads=40,
        )
        print(f"  completado en {time.time() - start:.1f}s — shape: {mi_matrix.shape}")

        mi_matrix.write_parquet(output_path, compression="zstd")
        print(f"  guardado en: {output_path}")

        matriz_a_edgelist(output_path)


if __name__ == "__main__":
    calculate_mi_protein_coding()