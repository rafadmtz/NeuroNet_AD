import pandas as pd
from pathlib import Path
from itertools import product
import numpy as np

GROUPS= ["Not_AD", "Low", "Intermediate", "High"]
BASE = Path("/export/space3/users/silvanac/NeuroNet_AD")
INDIR = Path(f"{BASE}/output/protein_coding/mi_matrices")
OUTFIGDIR = Path(f"{BASE}/output/protein_coding/figures")
OUTDATADIR = Path(f"{BASE}/output/protein_coding/data")

def jaccard_modules(list1, list2):
    set1, set2 = set(list1), set(list2)
    return len(set1 & set2) / len(set1 | set2)
    
def resumen_similitud(matriz_similitud):
    
    row_sim_scores = []

    for fila in matriz_similitud:
        valores_overlap = fila[fila > 0]
        if valores_overlap.size == 0:
            row_sim_scores.append(0)
        else:
            row_sim_scores.append(np.mean(valores_overlap))

    similarity_score = np.mean(row_sim_scores)

    matching_modules = (matriz_similitud > 0).any(axis=1).sum()
    perfect_matches  = (matriz_similitud == 1).any(axis=1).sum()

    return {
        "similarity_score": similarity_score,
        "matching_modules": matching_modules,
        "perfect_matches": perfect_matches,
    }
    
def matrices_similitud_by_group(main_group):
    
    df_main= pd.read_csv(OUTDATADIR / f"infomap_comunidades_{main_group}.csv")
    
    num_main= df_main["modulo"].unique()
    
    comunidades_main = {i: df_main.loc[df_main["modulo"] ==i, "gen"].tolist()
                        for i in num_main}
    
    medidas_similitud= {
        "main": main_group
    }
    
    for group in GROUPS:
        if main_group == group:
            continue
        
        df_similitud= pd.read_csv(OUTDATADIR / f"infomap_comunidades_{group}.csv")
        
        num_similitud= df_similitud["modulo"].unique()
        
        comunidades_similitud = {i: df_similitud.loc[df_similitud["modulo"] ==i, "gen"].tolist()
                    for i in num_similitud}
        
        matriz_similitud = np.zeros((len(num_main), len(num_similitud)))
        
        for i, j in product(num_main, num_similitud):
            matriz_similitud[i-1][j-1]= jaccard_modules(comunidades_main[i], comunidades_similitud[j])
        
        medidas_similitud[group]=resumen_similitud(matriz_similitud)    
        
    return medidas_similitud

def df_similitud_estadios():
      
    filas = []

    for group in GROUPS:
        resultado = matrices_similitud_by_group(group)
        main = resultado["main"]

        for comparado, medidas in resultado.items():
            if comparado == "main":
                continue
            filas.append({
                "referencia": main,
                "comparado": comparado,
                **medidas
            })

    df_similitud = pd.DataFrame(filas)
    
    df_similitud.to_csv(OUTDATADIR / "df_similitud_estadios.csv", index=False)
        
    return df_similitud

if __name__ == "__main__":
    df_similitud_estadios()