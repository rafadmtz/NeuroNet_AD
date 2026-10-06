library (dplyr)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO/unitotal/sin_simplificar")
OUTFIGDIR <- file.path(BASE, "output/protein_coding/figures/ora_enrichmentGO")

jaccard_enrichments_NS <- function(group1, group2){
    
    enrich1 <- read.csv(file.path(INDIR, paste0("enrichmentGO_unitotal_", group1, ".csv")))
    enrich2 <- read.csv(file.path(INDIR, paste0("enrichmentGO_unitotal_", group2, ".csv")))

    v_1 <- c(enrich1$ID)
    v_2 <- c(enrich2$ID)

    cat(group1, length(unique(v_1)), "terminos,", group2, length(unique(v_2)), "terminos\n")

    return(length(intersect(v_1, v_2)) / length(union(v_1, v_2)))
}

jacc_com <- function(com1, com2){
    return(length(intersect(com1, com2)) / length(union(com1, com2)))
}

jaccard_comunidades <- function(group1, group2){
    enrich1 <- read.csv(file.path(INDIR, paste0("enrichmentGO_unitotal_", group1, ".csv")))
    enrich2 <- read.csv(file.path(INDIR, paste0("enrichmentGO_unitotal_", group2, ".csv")))

    comunidades1 <- split(enrich1$ID, enrich1$Cluster)
    comunidades2 <- split(enrich2$ID, enrich2$Cluster)

    cat(group1, "vs", group2, "\n")
    cat("comunidades", group1, ":", length(comunidades1), "\n")
    cat("comunidades", group2, ":", length(comunidades2), "\n")

    m <- matrix(NA_real_, nrow = length(comunidades1), ncol = length(comunidades2)
                ,dimnames = list(names(comunidades1), names(comunidades2)))

    for (i in names(comunidades1)){
        for (j in names(comunidades2)){
            m[i,j] <- jacc_com(comunidades1[[i]], comunidades2[[j]])
        }   
    }
    return(m)
}

matriz_jaccard_funcional_red <- function(){
    
    m <- matrix(NA_real_, nrow = length(GROUPS), ncol = length(GROUPS)
                ,dimnames = list(GROUPS ,GROUPS))

    for (i in seq_along(GROUPS)){
        for (j in seq_along(GROUPS)){
            m[i,j] <- jaccard_enrichments_NS(GROUPS[i],GROUPS[j])
        }
    }

    cat("jaccard entre redes:\n")
    print(round(m, 3))

    return(m)
}

matrices_jaccard_comunidades <- function(referencia = "Not_AD") {
    matrices <- list()

    for (group in setdiff(GROUPS, referencia)) {
        matrices[[group]] <- jaccard_comunidades(referencia, group)
    }

    return(matrices)
}

resumen_matriz_funcional <- function(m, umbral = 0.5){

    j_max <- apply(m, 1, max)

    n_comunidades <- nrow(m)
    n_exactas <- sum(apply(m, 1, function(fila) any(fila == 1)))
    n_umbral <- sum(j_max >= umbral)

    return(data.frame(
        comunidades_ref = n_comunidades,
        con_J_1 = n_exactas,
        con_J_umbral = n_umbral,
        J_max_promedio = mean(j_max)
    ))
}

resumen_jaccard_comunidades <- function(matrices, umbral = 0.5){
    resumen <- list()

    for (group in names(matrices)){
        resumen[[group]] <- resumen_matriz_funcional(matrices[[group]], umbral)
        resumen[[group]]$comparado <- group
    }

    resumen <- do.call(rbind, resumen)

    cat("resumen jaccard entre comunidades:\n")
    print(resumen, row.names = FALSE)

    return(resumen)
}

function(){

    m_redes <- matriz_jaccard_funcional_red()
    write.csv(m_redes, file.path(INDIR, "jaccard_funcional_redes.csv"))

    matrices <- matrices_jaccard_comunidades()

    for (group in names(matrices)){
        write.csv(matrices[[group]], file.path(INDIR, paste0("jaccard_funcional_comunidades_Not_AD_vs_", group, ".csv")))
    }

    resumen_com <- resumen_jaccard_comunidades(matrices)
    write.csv(resumen_com, file.path(INDIR, "resumen_jaccard_comunidades.csv"), row.names = FALSE)
}

main()