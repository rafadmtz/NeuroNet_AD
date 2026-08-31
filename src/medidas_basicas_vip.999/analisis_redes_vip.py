from pathlib import Path
import polars as pl
import networkx as nx

ruta = Path("/home/rafadiaz/export/NeuroNet_AD/numba_mi_results_vip")
condiciones = ["Vip_High", "Vip_Intermediate", "Vip_Low", "Vip_Not_AD"]

for cond in condiciones:
    parquet_path = ruta / f"{cond}_mi.parquet"

    df = pl.read_parquet(parquet_path)

    # La primera columna es el nombre del gen (Regulator, sin nombre en el parquet)
    id_col = df.columns[0]
    df = df.rename({id_col: "Regulator"})

    # Convertir la matriz (formato ancho) a formato largo: Regulator, Target, MI
    df_largo = df.unpivot(
        index="Regulator",
        variable_name="Target",
        value_name="MI",
    )

    # Quitar diagonal (gen consigo mismo) y valores nulos
    df_largo = df_largo.filter(
        (pl.col("Regulator") != pl.col("Target")) & pl.col("MI").is_not_null()
    )

    # Como es matriz simétrica, cada par aparece 2 veces (A-B y B-A).
    # Quitamos duplicados quedándonos con un orden canónico del par.
    df_largo = df_largo.with_columns(
        pl.min_horizontal("Regulator", "Target").alias("_a"),
        pl.max_horizontal("Regulator", "Target").alias("_b"),
    ).unique(subset=["_a", "_b"]).drop(["_a", "_b"])

    umbral = df_largo.select(pl.col("MI").quantile(0.999)).item()
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

    print(f"\n=== {cond} ===")
    print(f"Umbral MI (q0.999): {umbral:.6f}")
    print(f"Nodos: {G.number_of_nodes():,}  |  Aristas: {G.number_of_edges():,}")
    print(f"Nodo central (mayor grado): {nodo_central}  (grado {degree[nodo_central]})")
    print(f"K-core maximo: {k_max}  ({len(kcore_nodes)} nodos)")
    print(f"Mayor grado en k-core: {top_kcore}  (grado {grados_kcore[top_kcore]})")
