library(arrow)
library(igraph)

BASE <- "/home/rafadiaz/export/NeuroNet_AD"

INDIR <- file.path(
    BASE,
    "output/red_mi_vip_completo/edgelist_vip_completo"
)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

for (group in GROUPS) {

    cat("Procesando:", group, "\n")

    # Edgelist (igual que en el script original)
    infile <- file.path(
        INDIR,
        paste0("edgelist_", group, "_MI_nodelink01.parquet")
    )

    edges <- read_parquet(
        infile,
        col_select = c("source", "target")
    )

    g <- graph_from_data_frame(
        edges,
        directed = FALSE
    )

    # Leer comunidades ya calculadas
    comm_file <- file.path(
        INDIR,
        paste0("infomap_", group, ".csv")
    )

    comm <- read.csv(comm_file)

    # Alinear orden de comm$community con V(g)$name
    membership_vec <- comm$community[match(V(g)$name, comm$gene)]

    Q <- modularity(g, membership_vec)

    cat("  Modularity Q:", Q, "\n")
}