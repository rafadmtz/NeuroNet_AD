import pandas as pd

def get_hvg_genes():
    condiciones = ["High", "Not_AD", "Intermediate", "Low"]
    N = 5000

    sets_hvg = {}

    for cond in condiciones:
        df = pd.read_csv(f"data/Vip_{cond}_metacells.tsv", sep="\t", index_col=0)
        cv = df.std(axis=1) / df.mean(axis=1)
        top = cv.sort_values(ascending=False).head(N).index
        sets_hvg[cond] = set(top)
        print(f"{cond}: {len(top)} HVGs")

    union = set.union(*sets_hvg.values())
    print(f"\nUnión total: {len(union)} genes")

    pd.Series(sorted(union)).to_csv("data/hvg_union.txt", index=False, header=False)
    print("Guardado en data/hvg_union.txt")

    print("\nSolapamiento entre condiciones:")
    for c1 in condiciones:
        for c2 in condiciones:
            if c1 < c2:
                overlap = len(sets_hvg[c1] & sets_hvg[c2])
                print(f"  {c1} ∩ {c2}: {overlap} genes")

def filter_hvg_genes():
    genes = pd.read_csv("data/hvg_union.txt", header=None)[0].tolist()
    print(f"HVGs en lista: {len(genes)}")

    for cond in ["High", "Not_AD", "Intermediate", "Low"]:
        df = pd.read_csv(f"data/Vip_{cond}_metacells.tsv", sep="\t", index_col=0)
        df_filtered = df.loc[df.index.isin(genes)]
        out = f"data/Vip_{cond}_hvg.tsv"
        df_filtered.to_csv(out, sep="\t")
        print(f"{cond}: {len(df)} → {len(df_filtered)} genes → {out}")