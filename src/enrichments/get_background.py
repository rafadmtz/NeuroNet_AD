
import gseapy as gp

go = gp.get_library(
    name="GO_Biological_Process_2023",
    organism="Human"
)

go_genes = set()

for genes in go.values():
    go_genes.update(genes)

with open("data/gene_lists/vip/all_genes_vip_ad.txt") as f:
    all_genes = {gene.strip() for gene in f}

background = all_genes & go_genes

print("Genes medidos:", len(all_genes))
print("Genes en GO BP:", len(go_genes))
print("Genes medidos con anotación GO BP:", len(background))
print("Genes medidos sin anotación GO BP:", len(all_genes - go_genes))