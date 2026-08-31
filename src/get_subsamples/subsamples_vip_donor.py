#!/usr/bin/env python
"""
Para cada replica: de la tabla de metadata, elige 7 donantes de cada grado
de severidad (ADNC). Busca el TSV de esos donantes en output/ y los
concatena, para quedar con 4 archivos (uno por grupo).

Uso:
    python subsample.py                 # 20 replicas
    python subsample.py --n-reps 30
    python subsample.py --rep 0         # solo una replica
"""

import argparse
import os

import numpy as np
import pandas as pd

BASE = "/home/rafadiaz/export/NeuroNet_AD"
META_CSV = f"{BASE}/data/metadata/overallAD.csv"
INPUT_DIR = f"{BASE}/output"
OUTDIR = f"{BASE}/output/subsamples_Vip"

CELL_TYPE = "Vip"
COL_DONOR = "Donor ID"
COL_ADNC = "Overall AD neuropathological Change"

GROUPS = ["Not AD", "Low", "Intermediate", "High"]
N_DONORS = 7

INPUT_ORIENTATION = "genes_x_cells"   # o "cells_x_genes"


def read_donor_tsv(donor):
    path = os.path.join(INPUT_DIR, f"{donor}_{CELL_TYPE}.tsv")
    df = pd.read_csv(path, sep="\t", index_col=0)
    if INPUT_ORIENTATION == "cells_x_genes":
        df = df.T
    return df


def filter_available(meta):
    """Se queda solo con donantes que sí tienen archivo {donor}_Vip.tsv."""
    has_file = meta[COL_DONOR].apply(
        lambda d: os.path.exists(os.path.join(INPUT_DIR, f"{d}_{CELL_TYPE}.tsv"))
    )
    faltantes = meta.loc[~has_file, COL_DONOR].tolist()
    if faltantes:
        print(f"Sin archivo {CELL_TYPE}, se excluyen ({len(faltantes)}): "
              f"{faltantes}")

    meta = meta[has_file].reset_index(drop=True)

    print("\nDonantes disponibles por grupo:")
    counts = meta[COL_ADNC].value_counts().reindex(GROUPS).fillna(0).astype(int)
    for g in GROUPS:
        flag = "" if counts[g] >= N_DONORS else "  <-- insuficiente"
        print(f"  {g:<13} {counts[g]:>3}{flag}")

    faltan = [g for g in GROUPS if counts[g] < N_DONORS]
    if faltan:
        raise SystemExit(f"\nSe piden {N_DONORS} donantes y no alcanzan en: {faltan}")

    return meta


def build_replicate(meta, rep, seed):
    rng = np.random.default_rng(seed + rep)

    for grp in GROUPS:
        donors = meta.loc[meta[COL_ADNC] == grp, COL_DONOR]
        chosen = rng.choice(donors.to_numpy(), size=N_DONORS, replace=False)
        print(f"  {grp:<13} {list(chosen)}")

        blocks = []
        for donor in chosen:
            df = read_donor_tsv(donor)
            df.columns = [f"{donor}|{c}" for c in df.columns]
            blocks.append(df)

        mat = pd.concat(blocks, axis=1).fillna(0)

        slug = grp.replace(" ", "_")
        out = os.path.join(OUTDIR, f"rep{rep:02d}_{slug}.tsv")
        mat.to_csv(out, sep="\t")
        print(f"    -> {out}  shape={mat.shape}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-reps", type=int, default=20)
    ap.add_argument("--rep", type=int, default=None,
                    help="genera solo esta replica")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    meta = pd.read_csv(META_CSV)[[COL_DONOR, COL_ADNC]].dropna()
    meta = filter_available(meta)

    reps = [args.rep] if args.rep is not None else range(args.n_reps)
    for r in reps:
        print(f"\n=== replica {r} ===")
        build_replicate(meta, r, args.seed)


if __name__ == "__main__":
    main()