
library(org.Hs.eg.db)
library(AnnotationDbi)

all_genes <- readLines(
  "data/gene_lists/vip/all_genes_vip_ad.txt"
)

ensembl_genes <- all_genes[grepl("^ENSG", all_genes)]
symbol_genes <- all_genes[!grepl("^ENSG", all_genes)]

mapped <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_genes,
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)

mapped_symbols <- mapped$SYMBOL[!is.na(mapped$SYMBOL)]

all_genes_symbols <- unique(
  c(symbol_genes, mapped_symbols)
)

go_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
  keytype = "SYMBOL",
  columns = c("GO", "ONTOLOGY")
)

go_bp_genes <- unique(
  go_annotations$SYMBOL[
    go_annotations$ONTOLOGY == "BP" &
    !is.na(go_annotations$GO)
  ]
)

background <- intersect(
  all_genes_symbols,
  go_bp_genes
)

cat("Genes medidos:", length(all_genes), "\n")
cat("Genes medidos convertidos a SYMBOL:", length(all_genes_symbols), "\n")
cat("Genes anotados en GO BP:", length(go_bp_genes), "\n")
cat("Background final:", length(background), "\n")

writeLines(
  background,
  "data/gene_lists/vip/all_genes_with_GO_vip_ad.txt"
)