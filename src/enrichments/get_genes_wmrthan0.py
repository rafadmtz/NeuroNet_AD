from pathlib import Path
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt

INDIR = Path("data") / "count_matrices" / "all_vip"
GROUPS = ["High", "Not_AD", "Intermediate", "Low"]

def get_genes_dif0():
    

    all_genes = set()

    for group in GROUPS:
        df = pl.read_csv(
            INDIR / f"Vip_{group}.tsv",
            separator="\t"
        )

        # En Polars no hay índices:
        # la primera columna contiene los genes
        gene_col = df.columns[0]

        print(group)
        print()
        print(f"genes {df.height}")
        print()

        # Conservar genes que tengan al menos un valor != 0
        # en alguna célula
        df = df.filter(
            pl.any_horizontal(
                pl.exclude(gene_col) != 0
            )
        )

        print(f"genes despues de eliminar 0s {df.height}")
        print()

        # Agregar genes a un único conjunto
        all_genes.update(
            df.get_column(gene_col).to_list()
        )

    print(f"genes totales {len(all_genes)}")

    OUTFILE = Path(
        "data/gene_lists/vip/all_genes_with_morethan0.txt"
    )

    OUTFILE.parent.mkdir(parents=True, exist_ok=True)

    with open(OUTFILE, "w") as f:
        for gene in sorted(all_genes):
            f.write(f"{gene}\n")

    return all_genes

def counts_distribution():
    GROUPS = ["High", "Not_AD", "Intermediate", "Low"]

    means_list = {}

    for group in GROUPS:
        df = pl.read_csv(
            INDIR / f"Vip_{group}.tsv",
            separator="\t"
        )

        # En Polars no hay índices:
        # la primera columna contiene los genes
        gene_col = df.columns[0]

        means_list[group] = df.select(
            pl.mean_horizontal(pl.exclude(gene_col)).alias("mean")
        )
    
    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    axes = axes.flatten()

    for ax, group in zip(axes, GROUPS):

        means = means_list[group].get_column("mean").to_numpy()

        ax.hist(
            means,
            bins=50
        )

        ax.set_title(group)
        ax.set_yscale("log")
        ax.set_xlabel("Mean expression per gene")
        ax.set_ylabel("Number of genes")

    plt.tight_layout()

    OUTFIG = Path(__file__).parent / "counts_distribution.png"

    plt.savefig(
        OUTFIG,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Figura guardada en: {OUTFIG}")
                
INDIR_MI = Path("numba_mi_results_vip")


def genes_with_mi():

    genes_mi = set()
    genes_total = set()

    for group in GROUPS:

        df = pl.read_parquet(
            INDIR_MI / f"Vip_{group}_mi.parquet"
        )

        gene_col = df.columns[0]

        # Todos los genes presentes en la matriz
        genes_group_total = set(
            df.get_column(gene_col).to_list()
        )

        genes_total.update(genes_group_total)

        # Genes que tienen al menos un MI > 0
        genes_keep = (
            df.filter(
                pl.sum_horizontal(
                    pl.exclude(gene_col)
                ) > 0
            )
            .get_column(gene_col)
            .to_list()
        )

        # Agregarlos al set global
        genes_mi.update(genes_keep)

        print(group)
        print(f"Genes totales: {len(genes_group_total)}")
        print(f"Genes con MI > 0: {len(genes_keep)}")
        print()

    print("------------------------------")
    print(f"Genes únicos totales: {len(genes_total)}")
    print(f"Genes únicos con MI > 0: {len(genes_mi)}")
    print(f"Genes sin ninguna conexión MI: {len(genes_total - genes_mi)}")

    return genes_mi


           

if __name__ == "__main__":
    genes_with_mi()
    