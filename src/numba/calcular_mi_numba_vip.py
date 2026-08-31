import os
import time

import polars as pl
from numba_mi import pipeline

BASE = "/home/rafadiaz/export/NeuroNet_AD"
INDIR = f"{BASE}/data/count_matrices/all_vip"
OUTDIR = f"{BASE}/output/red_mi_vip_completo/mi_matrices"
GROUPS = ["Not_AD", "Low", "Intermediate", "High"]

# celda de esquina vacia que numba_mi convierte en un "gen" fantasma
# conectado con todo (bug de rep00)
GHOST = {"", "Unnamed: 0", "Unnamed:_0", "None", "null"}


def limpiar_fantasma(df):
    id_col = df.columns[0]
    df = df.drop([c for c in df.columns[1:] if c.strip() in GHOST])
    return df.filter(
        ~pl.col(id_col).cast(pl.Utf8).str.strip_chars().is_in(list(GHOST))
        & pl.col(id_col).is_not_null()
    )


os.makedirs(OUTDIR, exist_ok=True)

for grp in GROUPS:
    file_path = f"{INDIR}/Vip_{grp}.tsv"
    output_path = f"{OUTDIR}/Vip_{grp}_MI.parquet"

    if os.path.exists(output_path):
        print(f"{grp}: ya existe, se omite")
        continue

    print(f"{grp}: calculando MI...")
    start = time.time()
    mi_matrix = pipeline.whole_process(
        file_path=file_path,
        axis=1,
        sep="\t",
        index_col=0,
        num_threads=40,
    )
    print(f"  completado en {time.time() - start:.1f}s — shape: {mi_matrix.shape}")

    mi_matrix = limpiar_fantasma(mi_matrix)
    mi_matrix.write_parquet(output_path, compression="zstd")
    print(f"  guardado en: {output_path}")