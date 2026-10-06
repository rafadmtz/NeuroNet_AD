library(igraph)
library(arrow)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR      <- file.path(BASE, "output/protein_coding/mi_matrices")
OUTDATADIR <- file.path(BASE, "output/protein_coding/data")

run_infomap <- function(){

    infomap_communities <- list()

    for (group in GROUPS) {
        path <- file.path(
        INDIR, "filtered_edgelists",
        paste0("Vip_", group, "_protein_coding_mi_edgelist_filtered01.parquet")
        )
        cat(group, ":", path, "\n")

        edgelist <- as.data.frame(read_parquet(path))

        G <- graph_from_data_frame(d= edgelist, directed = FALSE)

        infomap_communities[[group]] <- cluster_infomap(G, nb.trials = 1000)

    }
    cat("Q values por estadio\n")
    Q_values <- sapply(infomap_communities, function(x) modularity(x))
    print(round(Q_values, 3))

    cat("longitud de comunidades:\n")
    print(sapply(infomap_communities, length))

    for (group in GROUPS) {
        df <- data.frame(
                gen    = names(membership(infomap_communities[[group]])),
                modulo = as.integer(membership(infomap_communities[[group]]))
            )

        archivo <- file.path(OUTDATADIR, paste0("infomap_comunidades_", group, ".csv"))
                write.csv(df, archivo, row.names = FALSE)
    }

}

run_infomap()