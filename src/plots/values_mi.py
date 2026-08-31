import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = BASE / "output" / "red_mi_vip_completo" / "mi_matrices"
GROUPS = ["High"]

for group in GROUPS:
    mat = pl.read_parquet(INDIR / f"Vip_{group}_MI.parquet")
    id_col = mat.columns[0]

    arr = mat.drop(id_col).to_numpy()  # matriz numérica pura, sin la columna de nombres

    # índices del triángulo superior, sin la diagonal (k=1 excluye la diagonal)
    iu = np.triu_indices_from(arr, k=1)

    valores_mi = arr[iu]  # array 1D con un solo valor por par, sin diagonal, sin duplicados

    # a Series de polars si quieres seguir con .describe(), .hist(), etc.
    valores_mi = pl.Series("MI", valores_mi).drop_nulls()

    valores_mi.describe()
    valores_mi.hist()

    plt.hist(valores_mi.to_numpy(), bins=100)
    plt.yscale("log")
    plt.xlabel("MI")
    plt.ylabel("Frecuencia")
    plt.title("Distribución de valores de MI")
    plt.savefig(f"src/figures/Vip_{group}.png", dpi=150, bbox_inches="tight")
    plt.close()