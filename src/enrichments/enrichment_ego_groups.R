library(arrow)
library(clusterProfiler)
library(org.Hs.eg.db)

GROUPS <- c("High", "Not_AD", "Intermediate", "Low")

GENE_LIST_FILE <- "data/gene_lists/vip/all_genes_with_GO_vip_ad.txt"

INDIR <- "output/red_mi_vip_completo/ego_groups"

OUTDIR <- file.path(INDIR, "enrichment_GO_BP")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# -------------------------
# Background
# -------------------------

background <- readLines(GENE_LIST_FILE)
background <- unique(background)

cat("Genes en background:", length(background), "\n\n")


# -------------------------
# Enrichment por ego group
# -------------------------

for (group in GROUPS) {

    cat("====================================\n")
    cat("Grupo:", group, "\n")

    # Leer red
    df <- read_parquet(
        file.path(
            INDIR,
            paste0("ego01_SR_SUG_LUC_", group, ".parquet")
        )
    )

    # Obtener todos los nodos del ego group
    genes <- unique(
        c(
            df$source,
            df$target
        )
    )

    cat("Genes en el ego group:", length(genes), "\n")

    # Revisar si hay genes del ego group fuera del background
    genes_not_background <- setdiff(
        genes,
        background
    )

    cat(
        "Genes fuera del background:",
        length(genes_not_background),
        "\n"
    )


    # -------------------------
    # GO Biological Process
    # -------------------------

    ego <- enrichGO(
        gene          = genes,
        universe      = background,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = TRUE
    )


    # -------------------------
    # Resultados
    # -------------------------

    results <- as.data.frame(ego)

    cat(
        "Términos significativos:",
        nrow(results),
        "\n\n"
    )


    # -------------------------
    # Guardar
    # -------------------------

    write.csv(
        results,
        file.path(
            OUTDIR,
            paste0(
                "GO_BP_ego01_SR_SUG_LUC_",
                group,
                ".csv"
            )
        ),
        row.names = FALSE
    )
}

