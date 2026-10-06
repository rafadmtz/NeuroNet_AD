library(clusterProfiler)
library(org.Hs.eg.db)
library(arrow)
library(GOSemSim)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR      <- file.path(BASE, "output/protein_coding/data/comunidades_infomap")
OUTDATADIR <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO")
GENEINDIR <- file.path(BASE, "output/protein_coding/mi_matrices")

genes_universo_por_red <- function(group){
    genes <- read.csv(file.path(INDIR, paste0("infomap_comunidades_", group, ".csv")))

    genes <- sub("^MT\\.", "MT-", genes$gen)

    return(unique(genes))
}

genes_universo_medidos_total <- function(group) {
    path <- file.path(GENEINDIR, paste0("Vip_", group, "_protein_coding_mi.parquet"))

    genes <- setdiff(open_dataset(path)$schema$names, c("", "Unnamed: 0"))
    genes <- sub("^MT\\.", "MT-", genes)

    return(unique(genes))
}

ora_enrichmentGO <- function(min_genes, funcion_universo){

    resultados_por_grupo <- list()

    for (group in GROUPS){
        df_comunidades <- read.csv(file.path(INDIR, paste0("infomap_comunidades_", group, ".csv")))

        df_comunidades$gen <- sub("^MT\\.", "MT-", df_comunidades$gen)

        lista_comunidades <- split(df_comunidades$gen, df_comunidades$modulo)

        lista_comunidades <- lista_comunidades[lengths(lista_comunidades) >= min_genes]

        resultados <- list()

        universo <- funcion_universo(group)

        for(i in names(lista_comunidades)){
            resultados[[i]] <- enrichGO(
                gene= lista_comunidades[[i]],
                OrgDb= org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                universe= universo,
                pvalueCutoff = 0.05,
                qvalueCutoff = 1
                )
        }
        resultados_por_grupo[[group]] <- resultados
    }
    return(resultados_por_grupo)    
}

compare_clusters_GO <- function(universo, group, min_genes) {

    df_comunidades <- read.csv(file.path(INDIR, paste0("infomap_comunidades_", group, ".csv")))

    df_comunidades$gen <- sub("^MT\\.", "MT-", df_comunidades$gen)

    lista_comunidades <- split(df_comunidades$gen, df_comunidades$modulo)

    lista_comunidades <- lista_comunidades[lengths(lista_comunidades) >= min_genes]

    result <- compareCluster(
        geneClusters  = lista_comunidades,
        fun           = "enrichGO",
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 1,
        universe      = universo
    )

    return(result)
}

simplify_results <- function(resultado_go, cutoff = 0.7, sem_data = NULL) {
    simplified_results <- clusterProfiler::simplify(
        resultado_go,
        cutoff     = cutoff,
        by         = "p.adjust",
        select_fun = min,
        measure    = "Wang",
        semData    = sem_data
    )
}

main <- function(){
    go_bp_data <- godata(annoDb = "org.Hs.eg.db", ont = "BP")

    for (group in GROUPS){

        universo_total <- genes_universo_medidos_total(group)
        universo_red <- genes_universo_por_red(group)

        resultados_go_total <- compare_clusters_GO(universo_total, group, min_genes= 5)
        resultados_go_red <- compare_clusters_GO(universo_red, group, min_genes= 5)

        simplified_total <- simplify_results(resultados_go_total, sem_data = go_bp_data)
        simplified_red <- simplify_results(resultados_go_red, sem_data = go_bp_data)

        saveRDS(resultados_go_total, file.path(OUTDATADIR, paste0("enrichmentGO_unitotal_", group, ".rds")))
        saveRDS(resultados_go_red,   file.path(OUTDATADIR, paste0("enrichmentGO_unired_", group, ".rds")))
        saveRDS(simplified_total,    file.path(OUTDATADIR, paste0("enrichmentGO_unitotal_simplified_", group, ".rds")))
        saveRDS(simplified_red,      file.path(OUTDATADIR, paste0("enrichmentGO_unired_simplified_", group, ".rds")))

        archivo_total <- file.path(OUTDATADIR, paste0("enrichmentGO_unitotal_", group, ".csv"))
        archivo_red <- file.path(OUTDATADIR, paste0("enrichmentGO_unired_", group, ".csv"))
        archivo_total_simpli <- file.path(OUTDATADIR, paste0("enrichmentGO_unitotal_simplified_", group, ".csv"))
        archivo_red_simpli <- file.path(OUTDATADIR, paste0("enrichmentGO_unired_simplified_", group, ".csv"))

        write.csv(resultados_go_total, archivo_total, row.names = FALSE)
        write.csv(resultados_go_red, archivo_red, row.names = FALSE)
        write.csv(simplified_total, archivo_total_simpli, row.names = FALSE)
        write.csv(simplified_red, archivo_red_simpli, row.names = FALSE)
    }
}

main()