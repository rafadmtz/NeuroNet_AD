library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)

select <- dplyr::select

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO/unitotal/simplificado")
INDIR_NS <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO/unitotal/sin_simplificar")
OUTFIGDIR <- file.path(BASE, "output/protein_coding/figures/ora_enrichmentGO")
COMUNIDADESDIR <- file.path(BASE, "output/protein_coding/data/comunidades_infomap")
OUTFIGDIRTREE <- file.path(BASE, "output/protein_coding/figures/ora_enrichmentGO/treeplots")

dir.create(OUTFIGDIRTREE, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTFIGDIR, showWarnings = FALSE, recursive = TRUE)

read_simplified_unitotal_enrichmentGO <- function(group){
    return(readRDS(file.path(INDIR, paste0("enrichmentGO_unitotal_simplified_", group, ".rds"))))
}

read_not_simplified_unitotal_enrichmentGO <- function(group){
    return(readRDS(file.path(INDIR_NS, paste0("enrichmentGO_unitotal_", group, ".rds"))))
}

dotplot_enrichmentGO <- function(n_terms, n_counts){
    for (group in GROUPS){
        enrichment <- read_simplified_unitotal_enrichmentGO(group)

        enrichment <- filter(enrichment, Count >=n_counts)

        plot_enrich <- dotplot(enrichment, showCategory = n_terms)

        plot_enrich <- plot_enrich +
            ggtitle(group) +
            scale_fill_gradient(
                low = "#D73027",
                high = "#4575B4",
                trans = "log10",
                limits = c(1e-20, 0.05),
                breaks = c(1e-20, 1e-15, 1e-10, 1e-5, 0.05),
                oob = scales::squish
            )

        ruta <- file.path(OUTFIGDIR, paste0("dotplot_",n_terms,"terminos_",group,".png"))

        ggsave(ruta, plot = plot_enrich, width = 12, height = 12)
    }
}

topx_comunidades_infomap <- function(group, x_coms){
    comunidades <- read.csv(file.path(COMUNIDADESDIR, paste0("infomap_comunidades_",group,".csv")))
    split_comunidades <- split(comunidades$gen, comunidades$modulo)

    topx <- split_comunidades %>% lengths() %>% sort(decreasing = TRUE) %>% head(x_coms) 

    print(topx)

    return(names(topx))
}

dotplot_topx_comunidades <- function(n_terms, top_coms){
    for (group in GROUPS){
        enrichment <- read_not_simplified_unitotal_enrichmentGO(group)

        enrichment <- filter(enrichment, Count >= 3,
                        Cluster %in% topx_comunidades_infomap(group, x_coms=top_coms))

        plot_enrich <- dotplot(enrichment, showCategory = n_terms)

        ruta <- file.path(OUTFIGDIR, paste0("dotplot_",top_coms,"_terminos",n_terms, "_",group,".png"))

        plot_enrich <- plot_enrich + ggtitle(group)

        ggsave(ruta, plot = plot_enrich, width = 12, height = 12)
    }
}

# tree_plot_by_state <- function(group, d, top_coms){

#     enrichment <- read_simplified_unitotal_enrichmentGO(group)

#     enrichment <- filter(enrichment, Count >= 3,
#                         Cluster %in% topx_comunidades_infomap(group, x_coms=top_coms))

#     enrichment_sim <- pairwise_termsim(enrichment, method = "Wang", semData = d)

#     tree_plot <- treeplot(enrichment_sim, showCategory = 30, nCluster = 5, cluster_panel = "dotplot")

#     ruta <- file.path(OUTFIGDIRTREE, paste0("treeplot_",group,".png"))

#     tree_plot <- tree_plot + ggtitle(group)

#     ggsave(ruta, plot = tree_plot, width = 12, height = 12)
# }

get_all_conditions_csv <- function(){

    all_conditions <- data.frame()

    for (group in GROUPS){

        enrichment <- read_not_simplified_unitotal_enrichmentGO(group)

        enrichment_df <- as.data.frame(enrichment)

        enrichment_df <- enrichment_df %>%
                        filter(Count >= 3) %>%
                        select(Cluster, ID, Description, p.adjust, Count, geneID) %>%
                        mutate(Grupo = group)

        all_conditions <- bind_rows(all_conditions, enrichment_df)
    }
    return(all_conditions)
}

main <- function(){
    dotplot_enrichmentGO(n_terms= 3, n_counts=3)

    dotplot_topx_comunidades(n_terms= 15, top_coms= 3)

    d <- godata('org.Hs.eg.db', ont="BP")

    # for (group in GROUPS){
    #     tree_plot_by_state(group= group, d= d, top_coms= 3)
    # }
}

main()