import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"

overallAD_df = pd.read_csv(DATA_DIR / "AD_state_fem_conteo_sorted.csv")

for condicion in ["Not AD", "Low", "Intermediate", "High"]:
    donors = overallAD_df[overallAD_df["Overall AD neuropathological Change"] == condicion]["Donor ID"]
    
    dfs = [pd.read_csv(OUTPUT_DIR / f"{donor}_Vip.tsv", sep="\t", index_col=0)
           for donor in donors
           if (OUTPUT_DIR / f"{donor}_Vip.tsv").exists()]
    
    if dfs:
        merged = pd.concat(dfs, axis=1, join="inner")
        nombre = condicion.replace(" ", "_")
        merged.to_csv(DATA_DIR / f"Vip_{nombre}.tsv", sep="\t")
        print(f"{condicion}: guardado con {merged.shape[0]} genes y {merged.shape[1]} células")
