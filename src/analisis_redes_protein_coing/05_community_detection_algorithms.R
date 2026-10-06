library(igraph)
library(arrow)
library(aricode)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR      <- file.path(BASE, "output/protein_coding/mi_matrices")
OUTFIGDIR  <- file.path(BASE, "output/protein_coding/figures")
OUTDATADIR <- file.path(BASE, "output/protein_coding/data")


ami_comparacion <- function(){
  ranking_ami <- list()
  for (group in GROUPS) {
    path <- file.path(
      INDIR, "filtered_edgelists",
      paste0("Vip_", group, "_protein_coding_mi_edgelist_filtered01.parquet")
    )
    cat(group, ":", path, "\n")

    edgelist <- as.data.frame(read_parquet(path))

    G <- graph_from_data_frame(d= edgelist, directed = FALSE)

    comunidades <- list() 

    set.seed(2)

    comunidades[["infomap"]] <- cluster_infomap(G, nb.trials = 1000)

    comunidades[["fast_greedy"]] <- cluster_fast_greedy(G)

    comunidades[["louvain"]] <- cluster_louvain(G)

    comunidades[["walktrap"]] <- cluster_walktrap(G)

    comunidades[["eigen"]] <- cluster_leading_eigen(G)

    comunidades[["prop"]] <- cluster_label_prop(G)

    comunidades[["betweenness"]] <- cluster_edge_betweenness(G)

    ami_matrix <- matrix(NA_real_, nrow = length(comunidades), ncol = length(comunidades)
                  , dimnames = list(names(comunidades), names(comunidades)))

    cat("Q values para cada algoritmo\n")
    Q_values <- sapply(comunidades, function(x) modularity(x))
    print(round(Q_values, 3))

    cat("longitud de comunidades:\n")
    print(sapply(comunidades, length))

    community_vectors <- lapply(comunidades, function(x) as.integer(membership(x)))

    stopifnot(all(sapply(community_vectors, length) == vcount(G)))

    for (p in combn(names(community_vectors), 2, simplify = FALSE)) {
      ami_value <- AMI(community_vectors[[p[1]]], community_vectors[[p[2]]])     
      ami_matrix[p[1], p[2]] <- ami_value
      ami_matrix[p[2], p[1]] <- ami_value
    }

    diag(ami_matrix) <- NA

    cat("matriz de ami\n")
    print(round(ami_matrix, 2))

    ranking_ami[[group]] <- sort(rowMeans(ami_matrix, na.rm = TRUE), decreasing = TRUE)

    cat("valores de ami para cada algoritmo:\n")
    print(round(ranking_ami[[group]], 3))
  }

  ami_tabla <- sapply(ranking_ami, function(x) x[sort(names(x))])

  ranking_final <- sort(rowMeans(ami_tabla), decreasing = TRUE)

  cat("AMI promedio de cada algoritmo en las 4 condiciones:\n")
  print(round(ranking_final, 3))

  write.csv(ami_tabla, file.path(OUTDATADIR, "ami_medio_por_grupo.csv"))
}



comparar_estabilidad <- function(GROUPS, INDIR, n_semillas = 10, seed_base = 2) {

  estabilidad <- list()

  for (group in GROUPS) {
    path <- file.path(
      INDIR, "filtered_edgelists",
      paste0("Vip_", group, "_protein_coding_mi_edgelist_filtered01.parquet")
    )
    cat(group, ":", path, "\n")

    edgelist <- as.data.frame(read_parquet(path))
    G <- graph_from_data_frame(d = edgelist, directed = FALSE)

    algoritmos <- list(
      louvain = cluster_louvain,
      prop    = cluster_label_prop,
      infomap = function(g) cluster_infomap(g, nb.trials = 100)
    )

    # lista de corridas por algoritmo
    corridas <- list(louvain = list(), prop = list(), infomap = list())

    for (i in seq_len(n_semillas)) {
      set.seed(seed_base + i)
      for (nombre in names(algoritmos)) {
        comunidades <- as.integer(membership(algoritmos[[nombre]](G)))
        corridas[[nombre]][[i]] <- comunidades
      }
    }

    # AMI promedio entre las n_semillas corridas de cada algoritmo
    ami_por_algoritmo <- sapply(names(corridas), function(nombre) {
      vals <- c()
      for (p in combn(seq_len(n_semillas), 2, simplify = FALSE)) {
        vals <- c(vals, AMI(corridas[[nombre]][[p[1]]], corridas[[nombre]][[p[2]]]))
      }
      mean(vals)
    })

    estabilidad[[group]] <- sort(ami_por_algoritmo, decreasing = TRUE)

    cat(group, "- estabilidad AMI (", n_semillas, "semillas):\n")
    print(round(estabilidad[[group]], 3))
  }

  estabilidad
}

ami_comparacion()
estabilidad <- comparar_estabilidad(GROUPS, INDIR)