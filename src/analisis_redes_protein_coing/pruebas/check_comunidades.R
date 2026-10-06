GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR      <- file.path(BASE, "output/protein_coding/data/comunidades_infomap")
OUTDATADIR <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO")

UMBRALES <- c(2, 3, 5, 8, 10, 15, 20)

check_umbrales <- function(group) {
    df <- read.csv(file.path(INDIR, paste0("infomap_comunidades_", group, ".csv")))
    tamanos <- table(df$modulo)

    cat("\n==", group, "==\n")
    cat("comunidades:", length(tamanos), "| genes:", nrow(df), "\n")
    cat("tamaño de comunidad -> min:", min(tamanos), "| mediana:", median(tamanos),
        "| max:", max(tamanos), "\n")
    cat("distribución de tamaños (tamaño: número de comunidades):\n")
    print(table(as.integer(tamanos)))

    resumen <- data.frame(
        grupo              = group,
        min_genes          = UMBRALES,
        comunidades        = sapply(UMBRALES, function(u) sum(tamanos >= u)),
        pct_comunidades    = sapply(UMBRALES, function(u) round(100 * mean(tamanos >= u), 1)),
        genes_cubiertos    = sapply(UMBRALES, function(u) sum(tamanos[tamanos >= u])),
        pct_genes_cubiertos = sapply(UMBRALES, function(u) round(100 * sum(tamanos[tamanos >= u]) / nrow(df), 1))
    )
    print(resumen, row.names = FALSE)

    resumen
}

main <- function() {
    resumenes <- lapply(GROUPS, check_umbrales)
    tabla <- do.call(rbind, resumenes)

    write.csv(tabla, file.path(OUTDATADIR, "check_umbrales_comunidades.csv"), row.names = FALSE)
}

main()