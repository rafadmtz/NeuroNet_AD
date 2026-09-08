import pandas as pd
from pathlib import Path

gene_annotation = pd.read_csv("data/gene_lists/vip/all_genes_annotation.csv")

biotypes= gene_annotation["gene_biotype"].unique()

OUTPUT= Path("data/gene_lists/vip/annotations")
OUTPUT.mkdir(parents=True, exist_ok=True)

for bt in biotypes:
    gene_biotype= gene_annotation[gene_annotation["gene_biotype"]==bt]
    
    print(f"{bt}, total de {len(gene_biotype)}")
    
    gene_biotype.to_csv(OUTPUT / f"gene_annotation_{bt}.csv", index= False)
    