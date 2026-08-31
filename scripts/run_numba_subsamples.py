import argparse
import glob
import os
import time

from numba_mi import pipeline

BASE = "/home/rafadiaz/export/NeuroNet_AD"
SUBDIR = f"{BASE}/output/subsamples_Vip"
MIDIR = f"{BASE}/output/mi_Vip"

parser = argparse.ArgumentParser()
parser.add_argument("--threads", type=int, default=-1)
args = parser.parse_args()

os.makedirs(MIDIR, exist_ok=True)

files = sorted(glob.glob(os.path.join(SUBDIR, "rep*_*.tsv")))
print(f"{len(files)} archivos de subsample encontrados\n")

for path in files:
    fname = os.path.basename(path)
    output_path = os.path.join(MIDIR, fname.replace(".tsv", "_MI.tsv"))

    if os.path.exists(output_path):
        print(f"[{fname}] ya existe, se omite")
        continue

    print(f"[{fname}] calculando MI...")
    start = time.time()
    mi_matrix = pipeline.whole_process(
        file_path=path,
        axis=1,
        sep="\t",
        index_col=0,
        num_threads=args.threads,
    )
    print(f"  Completado en {time.time() - start:.1f}s — shape: {mi_matrix.shape}")

    mi_matrix.write_csv(output_path, separator="\t")
    print(f"  Guardado en: {output_path}\n")