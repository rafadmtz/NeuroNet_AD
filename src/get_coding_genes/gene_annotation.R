library(AnnotationHub)
library(ensembldb)

annotation_hub <- AnnotationHub()

edatabase <- annotation_hub[["AH119325"]]

all_genes <- trimws(readLines("data/gene_lists/vip/all_genes_vip_ad.txt"))

annot <- genes(
    edatabase,
    columns = c(
        "gene_id",
        "gene_name",
        "gene_biotype"
    )
)

annot <- as.data.frame(annot)

annot_filtered <- annot[annot$gene_name %in% all_genes | annot$gene_id %in% all_genes, ]
annot_filtered <- unique(annot_filtered)
annot_filtered <- annot_filtered[, c("gene_id", "gene_name", "gene_biotype")]


write.csv(
  annot_filtered,
  "data/gene_lists/vip/all_genes_annotation.csv",
  row.names = FALSE
)

length(all_genes)
nrow(annot_filtered)
length(unique(annot_filtered$gene_id))
length(unique(annot_filtered$gene_name))
sum(duplicated(annot_filtered$gene_id))
sum(duplicated(annot_filtered$gene_name))