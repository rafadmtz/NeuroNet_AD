library(clusterProfiler)
library(org.Hs.eg.db)

GROUPS <- c("Low", "Intermediate", "High")

GENE_LIST_FILE <- "data/gene_lists/vip/all_genes_with_GO_vip_ad.txt"

INDIR <- "output/red_mi_vip_completo/ego_groups"

OUTDIR <- file.path(INDIR, "enrichment_GO_BP_new_genes")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# Background

background <- readLines(GENE_LIST_FILE)
background <- unique(background)

cat("Genes en background:", length(background), "\n\n")


# Leer tabla de genes nuevos

new_genes_df <- read.csv(
    file.path(INDIR, "new_not_reappeared_genes_SR_SUG_LUC.csv"),
    stringsAsFactors = FALSE
)

# Parsear la columna new_genes (guardada por pandas como string tipo "['A', 'B']")
parse_gene_list <- function(x) {
    x <- gsub("\\[|\\]|'", "", x)
    genes <- strsplit(x, ",\\s*")[[1]]
    genes[genes != ""]
}


# Enrichment por grupo (genes nuevos)

for (group in GROUPS) {

    cat("Grupo:", group, "\n")

    genes <- parse_gene_list(
        new_genes_df$new_genes[new_genes_df$group == group]
    )

    cat("Genes nuevos en este grupo:", length(genes), "\n")

    # Revisar si hay genes fuera del background
    genes_not_background <- setdiff(
        genes,
        background
    )

    cat(
        "Genes fuera del background:",
        length(genes_not_background),
        "\n"
    )



    # GO Biological Process

    ego <- enrichGO(
        gene          = genes,
        universe      = background,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 1,
        readable      = TRUE
    )



    # Resultados
    results <- as.data.frame(ego)

    cat(
        "Términos significativos:",
        nrow(results),
        "\n\n"
    )



    write.csv(
        results,
        file.path(
            OUTDIR,
            paste0(
                "GO_BP_new_genes_SR_SUG_LUC_",
                group,
                ".csv"
            )
        ),
        row.names = FALSE
    )
}