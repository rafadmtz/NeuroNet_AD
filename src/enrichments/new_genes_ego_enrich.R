library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

GROUPS <- c("Low", "Intermediate", "High")

GENE_LIST_FILE <- "data/gene_lists/vip/all_genes_with_GO_vip_ad.txt"

INDIR <- "output/red_mi_vip_completo/ego_groups"

OUTDIR <- file.path(INDIR, "enrichment_GO_BP_new_genes")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

FIGDIR <- file.path(INDIR, "figures")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)


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


# Acumulador para el dotplot comparativo
all_results_simplified <- list()


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
        pvalueCutoff  = 1,
        qvalueCutoff  = 1,
        readable      = TRUE
    )



    # Resultados
    results <- as.data.frame(ego)
    results <- results[order(results$p.adjust), ]

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


  
    # Simplify (deduplicacion semantica) para el dotplot


    ego_simpl <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

    cat("Términos tras simplify:", nrow(as.data.frame(ego_simpl)), "\n\n")

    all_results_simplified[[group]] <- as.data.frame(ego_simpl) %>%
        mutate(group = group)
}



# Dotplot comparativo (Low / Intermediate / High)


all_results <- bind_rows(all_results_simplified)

write.csv(
    all_results,
    file.path(OUTDIR, "GO_BP_new_genes_SR_SUG_LUC_simplified_all.csv"),
    row.names = FALSE
)

top8_by_group <- all_results %>%
    group_by(group) %>%
    slice_min(order_by = p.adjust, n = 8) %>%
    pull(Description) %>%
    unique()

shared_significant <- all_results %>%
    filter(p.adjust < 0.05) %>%
    group_by(Description) %>%
    filter(n_distinct(group) >= 2) %>%
    pull(Description) %>%
    unique()

terms_to_plot <- union(top8_by_group, shared_significant)

plot_df <- all_results %>%
    filter(Description %in% terms_to_plot) %>%
    mutate(
        sig = ifelse(p.adjust < 0.05, "p.adjust < 0.05", "nominal (no significativo)"),
        group = factor(group, levels = GROUPS)
    )

term_order <- plot_df %>%
    group_by(Description) %>%
    summarise(best_p = min(p.adjust)) %>%
    arrange(best_p) %>%
    pull(Description)

plot_df$Description <- factor(plot_df$Description, levels = rev(term_order))

p <- ggplot(plot_df, aes(x = group, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust, shape = sig)) +
    scale_shape_manual(values = c("p.adjust < 0.05" = 16, "nominal (no significativo)" = 1)) +
    scale_color_gradient(low = "#e74c3c", high = "#1f77b4") +
    labs(
        x = "Grupo de severidad",
        y = NULL,
        size = "Genes (Count)",
        color = "p.adjust",
        shape = "Significancia",
        title = "GO Biological Process — genes nuevos por estadio (SR_SUG_LUC)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.major.y = element_line(color = "grey90"),
        axis.text.y = element_text(size = 9)
    )

ggsave(
    file.path(FIGDIR, "dotplot_new_genes_GO_BP_simplified.png"),
    p, width = 9, height = 8, dpi = 300
)