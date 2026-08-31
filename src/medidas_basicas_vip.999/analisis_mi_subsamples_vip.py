import argparse
import time
from pathlib import Path

import networkx as nx
import polars as pl

MIDIR = Path("/home/rafadiaz/export/NeuroNet_AD/output/mi_Vip")
RESULTS_CSV = MIDIR / "metrics.csv"

GROUPS = ["Not_AD", "Low", "Intermediate", "High"]

parser = argparse.ArgumentParser()
parser.add_argument("--n-reps", type=int, default=20,
                     help="numero total de replicas esperadas")
parser.add_argument("--wait-seconds", type=int, default=60,
                     help="cuanto dormir si un archivo aun no existe")
parser.add_argument("--quantile", type=float, default=0.999)
args = parser.parse_args()


def wait_for_file(path):
    """Duerme mientras el archivo no exista (numba_mi puede seguir corriendo)."""
    while not path.exists():
        print(f"  {path.name} no existe todavia, durmiendo "
              f"{args.wait_seconds}s...")
        time.sleep(args.wait_seconds)


def already_done(results_csv, rep, group):
    if not results_csv.exists():
        return False
    done = pl.read_csv(results_csv)
    return ((done["rep"] == rep) & (done["group"] == group)).any()


def analyze(path, quantile):
    parquet_path = path.with_suffix(".parquet")

    if parquet_path.exists():
        df = pl.read_parquet(parquet_path)
    else:
        df = pl.read_csv(path, separator="\t")
        df.write_parquet(parquet_path)

    id_col = df.columns[0]
    df = df.rename({id_col: "Regulator"})

    # numba_mi convierte la celda vacia de la esquina del TSV en un "gen"
    # llamado "Unnamed: 0", que queda conectado a todo. Se quita como fila
    # y como columna. Filtrarlo aqui equivale a haberlo excluido desde el
    # inicio: la MI se calcula por pares, asi que su presencia no altera
    # los valores de los pares reales.
    FANTASMA = "Unnamed: 0"
    df = df.filter(pl.col("Regulator") != FANTASMA)
    if FANTASMA in df.columns:
        df = df.drop(FANTASMA)

    df_largo = df.unpivot(
        index="Regulator",
        variable_name="Target",
        value_name="MI",
    )

    df_largo = df_largo.filter(
        (pl.col("Regulator") != pl.col("Target")) & pl.col("MI").is_not_null()
    )

    umbral = df_largo.select(pl.col("MI").quantile(quantile)).item()
    df_filtrado = df_largo.filter(pl.col("MI") >= umbral)

    G = nx.from_pandas_edgelist(
        df_filtrado.to_pandas(),
        source="Regulator",
        target="Target",
        edge_attr="MI",
        create_using=nx.Graph(),
    )

    degree = dict(G.degree())
    nodo_central = max(degree, key=degree.get)

    kcore_num = nx.core_number(G)
    k_max = max(kcore_num.values())
    kcore_nodes = [n for n, k in kcore_num.items() if k == k_max]
    kcore_subgraph = G.subgraph(kcore_nodes)

    grados_kcore = dict(kcore_subgraph.degree())
    top_kcore = max(grados_kcore, key=grados_kcore.get)

    return {
        "umbral_mi": umbral,
        "n_nodos": G.number_of_nodes(),
        "n_aristas": G.number_of_edges(),
        "nodo_central": nodo_central,
        "grado_central": degree[nodo_central],
        "kcore_max": k_max,
        "kcore_n_nodos": len(kcore_nodes),
        "kcore_top_nodo": top_kcore,
        "kcore_top_grado": grados_kcore[top_kcore],
    }


def append_result(row):
    new_row = pl.DataFrame([row])
    if RESULTS_CSV.exists():
        existing = pl.read_csv(RESULTS_CSV)
        pl.concat([existing, new_row], how="diagonal_relaxed").write_csv(RESULTS_CSV)
    else:
        new_row.write_csv(RESULTS_CSV)


for rep in range(args.n_reps):
    for group in GROUPS:
        if already_done(RESULTS_CSV, rep, group):
            print(f"rep{rep:02d} {group}: ya en la tabla, se omite")
            continue

        path = MIDIR / f"rep{rep:02d}_{group}_MI.tsv"
        print(f"\nrep{rep:02d} {group}")
        wait_for_file(path)

        metrics = analyze(path, args.quantile)
        print(f"  nodos={metrics['n_nodos']:,}  aristas={metrics['n_aristas']:,}  "
              f"central={metrics['nodo_central']} (grado {metrics['grado_central']})  "
              f"kcore={metrics['kcore_max']}")

        row = {"rep": rep, "group": group}
        row.update(metrics)
        append_result(row)

print(f"\nTabla de resultados -> {RESULTS_CSV}")