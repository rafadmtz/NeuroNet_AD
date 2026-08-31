library(arrow)
library(igraph)

#comunidades para las redes con node to link ratio de 01

BASE <- "/home/rafadiaz/export/NeuroNet_AD"

INDIR <- file.path(
    BASE,
    "output/red_mi_vip_completo/edgelist_vip_completo"
)

OUTDIR <- file.path(
    BASE,
    "output/red_mi_vip_completo/edgelist_vip_completo"
)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

set.seed(2)

for (group in GROUPS) {

    cat("Procesando:", group, "\n")

    # Leer edgelist
    infile <- file.path(
        INDIR,
        paste0("edgelist_", group, "_MI_nodelink01.parquet")
    )

    edges <- read_parquet(
        infile,
        col_select = c("source", "target")
    )

    # Crear red no ponderada
    g <- graph_from_data_frame(
        edges,
        directed = FALSE
    )

    # Infomap no ponderado
    infomap <- cluster_infomap(g)

    # Gen -> comunidad
    result <- data.frame(
        gene = names(membership(infomap)),
        community = membership(infomap)
    )

    # Guardar CSV
    outfile <- file.path(
        OUTDIR,
        paste0("infomap_", group, ".csv")
    )

    write.csv(
        result,
        outfile,
        row.names = FALSE
    )

    cat(
        "  Nodos:", vcount(g),
        "| Aristas:", ecount(g),
        "| Comunidades:", length(infomap),
        "\n"
    )
}