import scanpy as sc
import scipy.sparse as sp
from pathlib import Path
import time

t0 = time.time()

# ============================================================
# 0. Carga completa del AnnData
# ============================================================
data_path = Path("/export/space3/users/silvanac/NeuroNet_AD/data/whole_taxonomy_MTG_AD.h5ad")
print("Cargando AnnData completo...", flush=True)
adata = sc.read_h5ad(data_path)
print(f"Cargado en {time.time()-t0:.1f}s. Shape: {adata.shape}", flush=True)


# ============================================================
# 1. Estado de la matriz de expresión (adata.X)
# ============================================================
print("\n=== Matriz de expresión (adata.X) ===", flush=True)
print("Tipo:", type(adata.X))
print("Dtype:", adata.X.dtype)
print("Es sparse:", sp.issparse(adata.X))
print("Min:", adata.X.min(), flush=True)
print("Max:", adata.X.max(), flush=True)


# ============================================================
# 2. Layers y raw
# ============================================================
print("\n=== Layers y raw ===", flush=True)
print("Layers disponibles:", list(adata.layers.keys()))

if adata.raw is not None:
    print("adata.raw existe.")
    print("Min raw:", adata.raw.X.min(), "Max raw:", adata.raw.X.max())
else:
    print("adata.raw es None")


# ============================================================
# 3. Sparsity de la matriz
# ============================================================
print("\n=== Sparsity ===", flush=True)
if sp.issparse(adata.X):
    nnz = adata.X.nnz
    total = adata.X.shape[0] * adata.X.shape[1]
    print(f"Sparsity: {(1 - nnz/total)*100:.2f}%")
else:
    print("La matriz no es sparse (formato denso)")


# ============================================================
# 4. Genes con expresión nula/casi nula en el subset VIP
#    (usando protein_coding + lncRNA de feature_type)
# ============================================================
print("\n=== Expresión de genes filtrados en células VIP ===", flush=True)

genes_export = adata.var[['feature_name', 'feature_type']].copy().reset_index()
genes_export.columns = ['gene_id', 'gene_name', 'gene_biotype']
genes_filtrados = genes_export[genes_export['gene_biotype'].isin(['protein_coding', 'lncRNA'])]

vip = adata.obs[adata.obs['Subclass'] == 'Vip'].copy()
vip_idx = vip.index
gene_ids_filtrados = genes_filtrados['gene_id']

subset_vip = adata[vip_idx, gene_ids_filtrados]
gene_means = subset_vip.X.mean(axis=0)
gene_means = gene_means.A1 if sp.issparse(subset_vip.X) else gene_means

n_genes_cero = (gene_means == 0).sum()
print(f"Genes con expresión promedio 0 en VIP: {n_genes_cero} de {len(gene_means)}", flush=True)


# ============================================================
# 5. Balance de células VIP por donante
# ============================================================
print("\n=== Balance de células VIP por donante ===", flush=True)

vip_por_donante = vip.groupby('donor_id').size().reset_index(name='n_celulas_vip')
vip_por_donante = vip_por_donante.sort_values('n_celulas_vip')
print(vip_por_donante)

print(f"\n=== Script completo en {time.time()-t0:.1f}s ===", flush=True)