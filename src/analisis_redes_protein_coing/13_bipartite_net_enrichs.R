library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(GOSemSim)
library(igraph)
library(tidyr)

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter


GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE       <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR_NS <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO/unitotal/sin_simplificar")
OUTFIGDIR <- file.path(BASE, "output/protein_coding/figures/ora_enrichmentGO")
OUTDATADIR <- file.path(BASE, "output/protein_coding/data/ora_enrichmentGO/bipartite")
COMUNIDADESDIR <- file.path(BASE, "output/protein_coding/data/comunidades_infomap")

dir.create(OUTFIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDATADIR, showWarnings = FALSE, recursive = TRUE)

HIGHLIGHT_THEMES <- c("regulation of RNA splicing" = "#e31a1c",
                      "regulation of mRNA catabolic process" = "#6a3d9a",
                      "RNA transport" = "#fb9a99",
                      "regulation of ubiquitin-dependent protein catabolic process" = "#1f78b4",
                      "protein polyubiquitination" = "#cab2d6",
                      "endosomal transport" = "#33a02c",
                      "vesicle localization" = "#b2df8a",
                      "regulation of DNA repair" = "#ff7f00",
                      "nuclear chromosome segregation" = "#fdbf6f",
                      "regulation of sister chromatid segregation" = "#ffff99")

PUBU_COLORS <- c("#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf",
                 "#3690c0", "#0570b0", "#045a8d", "#023858")

MITO_THEMES <- c("mitochondrial respiratory chain complex assembly",
                 "nucleoside phosphate biosynthetic process",
                 "respiratory electron transport chain",
                 "response to hypoxia")

MITO_LABEL <- "Mitochondrial energy metabolism"
OTHER_LABEL <- "Other themes"

PUBU_COLORS <- c("#fff7fb", "#ece7f2", "#d0d1e6", "#a6bddb", "#74a9cf",
                 "#3690c0", "#0570b0", "#045a8d", "#023858")

read_not_simplified_unitotal_enrichmentGO <- function(group){
    return(readRDS(file.path(INDIR_NS, paste0("enrichmentGO_unitotal_", group, ".rds"))))
}

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

get_go_themes <- function(all_conditions, d, k = 20){

    ids <- unique(all_conditions$ID)

    sim <- termSim(ids, ids, semData = d, method = "Wang")

    valid <- !apply(is.na(sim), 1, all)

    cat("Total terms:", length(ids), "\n")
    cat("Terms without similarity:", sum(!valid), "\n")

    sim <- sim[valid, valid]

    hc <- hclust(as.dist(1 - sim), method = "ward.D2")

    clusters <- cutree(hc, k = k)

    themes <- data.frame(ID = names(clusters), theme = clusters)

    themes <- themes %>%
            left_join(distinct(all_conditions, ID, Description), by = "ID")

    themes <- get_theme_labels(themes, sim)

    return(list(themes = themes, hc = hc, sim = sim))
}

get_theme_labels <- function(themes, sim){

    medoids <- sapply(split(themes$ID, themes$theme), function(x){
        x[which.max(rowMeans(sim[x, x, drop = FALSE]))]
    })

    labels <- data.frame(theme = as.integer(names(medoids)),
                         theme_label = themes$Description[match(medoids, themes$ID)])

    themes <- themes %>%
            select(-any_of("theme_label")) %>%
            left_join(labels, by = "theme")

    return(themes)
}

reassign_themes <- function(themes, check){

    themes <- themes %>%
            left_join(select(check, ID, best_theme), by = "ID") %>%
            mutate(theme = ifelse(is.na(best_theme), theme, best_theme)) %>%
            select(-best_theme)

    return(themes)
}

check_theme_assignment <- function(sim, themes){

    mean_sim <- sapply(split(themes$ID, themes$theme), function(x){
        rowMeans(sim[themes$ID, x, drop = FALSE])
    })

    best_theme <- as.integer(colnames(mean_sim)[apply(mean_sim, 1, which.max)])

    check <- themes %>%
            mutate(best_theme = best_theme,
                   own_sim = mean_sim[cbind(ID, as.character(theme))],
                   best_sim = apply(mean_sim, 1, max)) %>%
            filter(theme != best_theme)

    return(check)
}

get_theme_edges <- function(all_conditions, themes){

    cat("Rows without theme:", sum(!all_conditions$ID %in% themes$ID), "\n")

    edges <- all_conditions %>%
            inner_join(select(themes, ID, theme, theme_label), by = "ID") %>%
            group_by(Grupo, Cluster, theme, theme_label) %>%
            summarise(weight = n(), .groups = "drop")

    return(edges)
}

build_bipartite <- function(edges, group){

    group_edges <- edges %>%
                filter(Grupo == group) %>%
                mutate(from = paste0("C", Cluster), to = paste0("T", theme))

    vertices <- bind_rows(
        distinct(group_edges, name = from) %>% mutate(label = name, type = TRUE),
        distinct(group_edges, name = to, label = theme_label) %>% mutate(type = FALSE)
    )

    b <- graph_from_data_frame(select(group_edges, from, to, weight),
                               directed = FALSE,
                               vertices = vertices)

    return(b)
}

bipartite_metrics <- function(b, group){

    is_community <- V(b)$type
    k <- degree(b)

    metrics <- data.frame(Grupo = group,
                          communities = sum(is_community),
                          themes = sum(!is_community),
                          edges = ecount(b),
                          mean_themes_per_community = mean(k[is_community]),
                          mean_communities_per_theme = mean(k[!is_community]),
                          shared_themes = sum(k[!is_community] > 1))

    return(metrics)
}

theme_community_table <- function(edges){

    table <- edges %>%
            group_by(theme, theme_label, Grupo) %>%
            summarise(n_communities = n_distinct(Cluster), .groups = "drop") %>%
            pivot_wider(names_from = Grupo, values_from = n_communities, values_fill = 0) %>%
            select(theme, theme_label, all_of(GROUPS))

    return(as.data.frame(table))
}

plot_bipartites <- function(bipartites, themes, radius = 10, shrink = 0.55,
                            min_dist = 1.8, max_r = 0.7, curvature = 0.08){

    missing <- setdiff(c(names(HIGHLIGHT_THEMES), MITO_THEMES), themes$theme_label)
    cat("Highlighted themes not found:", length(missing), "\n")
    if (length(missing) > 0) print(missing)

    display_colors <- c(HIGHLIGHT_THEMES,
                        setNames("#b15928", MITO_LABEL),
                        setNames("#a6cee3", OTHER_LABEL))
    display_order <- names(display_colors)

    to_display <- function(label){
        case_when(label %in% names(HIGHLIGHT_THEMES) ~ label,
                  label %in% MITO_THEMES ~ MITO_LABEL,
                  TRUE ~ OTHER_LABEL)
    }

    to_level <- function(display){
        case_when(display == OTHER_LABEL ~ "other",
                  display == MITO_LABEL ~ "mito",
                  TRUE ~ "highlight")
    }

    spread_points <- function(xy, min_dist, max_r, iters = 500){
        n <- nrow(xy)
        if (n < 2) return(xy)
        for (it in seq_len(iters)){
            moved <- FALSE
            for (i in 1:(n - 1)) for (j in (i + 1):n){
                d <- xy[j, ] - xy[i, ]
                dist <- sqrt(sum(d^2))
                if (dist < min_dist){
                    if (dist < 1e-6){
                        d <- runif(2, -1, 1)
                        dist <- sqrt(sum(d^2))
                    }
                    shift <- (min_dist - dist) / 2 * d / dist
                    xy[i, ] <- xy[i, ] - shift
                    xy[j, ] <- xy[j, ] + shift
                    moved <- TRUE
                }
            }
            r <- sqrt(rowSums(xy^2))
            out <- r > max_r
            xy[out, ] <- xy[out, ] * max_r / r[out]
            if (!moved) break
        }
        xy
    }

    # ---- theme order: highlighted + mito by co-occurrence, "Other" in a compact arc at the bottom ----
    theme_info <- bind_rows(lapply(bipartites, function(b){
            data.frame(name = V(b)$name[!V(b)$type], label = V(b)$label[!V(b)$type])
        })) %>%
        distinct(name, .keep_all = TRUE) %>%
        mutate(display = to_display(label),
               level = to_level(display))

    incidence <- bind_rows(lapply(names(bipartites), function(g){
            b <- bipartites[[g]]
            is_comm <- setNames(V(b)$type, V(b)$name)
            el <- as_edgelist(b)
            data.frame(theme = ifelse(is_comm[el[, 1]], el[, 2], el[, 1]),
                       community = paste(g, ifelse(is_comm[el[, 1]], el[, 1], el[, 2])),
                       weight = E(b)$weight)
        })) %>%
        tidyr::pivot_wider(names_from = community, values_from = weight,
                           values_fill = 0, values_fn = sum)

    mat <- as.matrix(incidence[, -1])
    rownames(mat) <- incidence$theme

    main_themes <- theme_info$name[theme_info$level != "other"]
    m <- mat[main_themes, , drop = FALSE]
    m <- m / sqrt(rowSums(m^2))
    hc <- hclust(as.dist(1 - m %*% t(m)), method = "average")
    main_order <- rownames(m)[hc$order]

    other_order <- theme_info %>%
        dplyr::filter(level == "other") %>%
        arrange(label) %>%
        pull(name)

    n_main <- length(main_order)
    spacing <- c(rep(1, n_main), rep(0.4, length(other_order)))
    gap <- 1.5
    pos <- cumsum(c(0, head(spacing, -1)))
    if (length(other_order) > 0) pos[(n_main + 1):length(pos)] <- pos[(n_main + 1):length(pos)] + gap
    total <- max(pos) + gap

    other_mid <- if (length(other_order) > 0) mean(pos[(n_main + 1):length(pos)]) else total / 2
    theta0 <- -pi / 2 + 2 * pi * other_mid / total            # "Other" arc centered at the bottom
    angles <- theta0 - 2 * pi * pos / total

    theme_xy <- data.frame(name = c(main_order, other_order),
                           x = radius * cos(angles),
                           y = radius * sin(angles)) %>%
        left_join(dplyr::select(theme_info, name, label, display, level), by = "name")

    # ---- per-stage layout ----
    community_nodes <- data.frame()
    present_themes <- data.frame()
    links <- data.frame()

    for (group in names(bipartites)){

        b <- bipartites[[group]]
        is_comm <- setNames(V(b)$type, V(b)$name)

        el <- as_edgelist(b)
        e <- data.frame(from = el[, 1], to = el[, 2], weight = E(b)$weight) %>%
            mutate(community = ifelse(is_comm[from], from, to),
                   theme = ifelse(is_comm[from], to, from)) %>%
            left_join(dplyr::select(theme_xy, name, tx = x, ty = y), by = c("theme" = "name"))

        comm_xy <- e %>%
            group_by(community) %>%
            summarise(x = sum(weight * tx) / sum(weight),
                      y = sum(weight * ty) / sum(weight), .groups = "drop")

        comm_xy <- data.frame(community = V(b)$name[V(b)$type]) %>%
            left_join(comm_xy, by = "community") %>%
            mutate(x = coalesce(x * shrink, 0), y = coalesce(y * shrink, 0))

        set.seed(42)
        xyc <- spread_points(as.matrix(comm_xy[, c("x", "y")]),
                             min_dist = min_dist, max_r = max_r * radius)
        comm_xy$x <- xyc[, 1]
        comm_xy$y <- xyc[, 2]

        group_comm <- data.frame(name = V(b)$name, label = V(b)$label, type = V(b)$type,
                                 strength = strength(b)) %>%
            dplyr::filter(type) %>%
            left_join(dplyr::rename(comm_xy, name = community), by = "name") %>%
            mutate(Grupo = group)

        group_links <- e %>%
            dplyr::select(community, theme, weight) %>%
            left_join(dplyr::select(comm_xy, community, x, y), by = "community") %>%
            left_join(dplyr::select(theme_xy, theme = name, xend = x, yend = y, display, level),
                      by = "theme") %>%
            mutate(Grupo = group)

        group_themes <- theme_xy %>%
            dplyr::filter(name %in% V(b)$name) %>%
            mutate(Grupo = group)

        community_nodes <- bind_rows(community_nodes, group_comm)
        present_themes <- bind_rows(present_themes, group_themes)
        links <- bind_rows(links, group_links)
    }

    # themes absent in a stage, drawn as hollow placeholders
    absent_themes <- tidyr::crossing(theme_xy, Grupo = names(bipartites)) %>%
        anti_join(dplyr::select(present_themes, name, Grupo), by = c("name", "Grupo"))

    to_factor <- function(df) mutate(df, Grupo = factor(Grupo, levels = GROUPS))
    community_nodes <- to_factor(community_nodes)
    present_themes <- to_factor(present_themes)
    absent_themes <- to_factor(absent_themes)
    links <- to_factor(links)

    # ---- edge widths and transparency by level ----
    w_lim <- range(links$weight)
    w_breaks <- unique(round(seq(w_lim[1], w_lim[2], length.out = 4)))
    lw_range <- c(0.5, 5)

    links <- links %>%
        mutate(lw = scales::rescale(weight, to = lw_range, from = w_lim),
               lw = ifelse(level == "other", pmax(lw * 0.6, 0.3), lw),
               lw_casing = lw + 0.9,
               alpha = case_when(level == "highlight" ~ scales::rescale(weight, to = c(0.55, 1), from = w_lim),
                                 level == "mito" ~ scales::rescale(weight, to = c(0.3, 0.55), from = w_lim),
                                 TRUE ~ 0.45)) %>%
        arrange(weight)

    links_other <- dplyr::filter(links, level == "other")
    links_mito <- dplyr::filter(links, level == "mito")
    links_high <- dplyr::filter(links, level == "highlight")

    # label only communities connected to highlighted themes
    label_nodes <- community_nodes %>%
        semi_join(distinct(links_high, community, Grupo), by = c("name" = "community", "Grupo"))

    node_size <- c(highlight = 8, mito = 7, other = 4)

    p <- ggplot() +

        # absent themes
        geom_point(data = absent_themes,
                   aes(x = x, y = y, size = I(node_size[level])),
                   shape = 21, fill = "white", color = "grey80", stroke = 0.6) +

        # level 3: "Other" edges, grey and thin
        geom_curve(data = links_other,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = lw, alpha = alpha),
                   color = "grey70", curvature = curvature, lineend = "round", show.legend = FALSE) +

        # level 2: mitochondria, white casing + transparent color
        geom_curve(data = links_mito,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = lw_casing),
                   color = "white", curvature = curvature, lineend = "round", show.legend = FALSE) +
        geom_curve(data = links_mito,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = lw, color = display, alpha = alpha),
                   curvature = curvature, lineend = "round", show.legend = FALSE) +

        # level 1: highlighted themes, white casing + saturated color
        geom_curve(data = links_high,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = lw_casing),
                   color = "white", curvature = curvature, lineend = "round", show.legend = FALSE) +
        geom_curve(data = links_high,
                   aes(x = x, y = y, xend = xend, yend = yend, linewidth = lw, color = display, alpha = alpha),
                   curvature = curvature, lineend = "round") +

        # present themes: white halo + colored node with border
        geom_point(data = present_themes,
                   aes(x = x, y = y, size = I(node_size[level] + 2)),
                   shape = 21, fill = "white", color = "white") +
        geom_point(data = dplyr::filter(present_themes, level == "other"),
                   aes(x = x, y = y, fill = display, size = I(node_size["other"])),
                   shape = 21, color = "grey50", stroke = 0.5) +
        geom_point(data = dplyr::filter(present_themes, level == "mito"),
                   aes(x = x, y = y, fill = display, size = I(node_size["mito"])),
                   shape = 21, color = "grey35", stroke = 0.8) +
        geom_point(data = dplyr::filter(present_themes, level == "highlight"),
                   aes(x = x, y = y, fill = display, size = I(node_size["highlight"])),
                   shape = 21, color = "grey20", stroke = 1) +

        # communities
        geom_point(data = community_nodes,
                   aes(x = x, y = y, size = strength),
                   shape = 23, fill = "grey15", color = "white", stroke = 0.6) +
        ggrepel::geom_text_repel(data = label_nodes,
                                 aes(x = x, y = y, label = label),
                                 size = 3.2, fontface = "bold", color = "grey10",
                                 bg.color = "white", bg.r = 0.15,
                                 box.padding = 0.5, point.padding = 0.4,
                                 min.segment.length = 0.3, segment.color = "grey50",
                                 segment.size = 0.3, max.overlaps = Inf, seed = 42) +

        scale_color_manual(values = display_colors, guide = "none") +
        scale_alpha_identity() +
        scale_linewidth_identity(name = "GO terms per edge", guide = "legend",
                                 breaks = scales::rescale(w_breaks, to = lw_range, from = w_lim),
                                 labels = w_breaks) +
        scale_fill_manual(values = display_colors, breaks = display_order, name = "GO theme") +
        scale_size(range = c(3, 9), name = "GO terms per community") +
        scale_x_continuous(expand = expansion(mult = 0.06)) +
        scale_y_continuous(expand = expansion(mult = 0.06)) +
        coord_equal(clip = "off") +
        facet_wrap(~ Grupo, ncol = 2, scales = "fixed") +
        guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 5,
                                                                  color = "grey20", stroke = 0.9)),
               size = guide_legend(order = 2, override.aes = list(shape = 23, fill = "grey15",
                                                                  color = "white")),
               linewidth = guide_legend(order = 3, override.aes = list(color = "grey40", alpha = 1))) +
        theme_void() +
        theme(legend.position = "right",
              legend.key.width = unit(1.2, "cm"),
              strip.text = element_text(size = 14, face = "bold", margin = margin(b = 10)),
              panel.spacing = unit(2.5, "lines"),
              plot.margin = margin(20, 20, 20, 20))

    return(p)
}

plot_theme_heatmap <- function(edges, themes){

    theme_levels <- themes %>%
                distinct(theme, theme_label) %>%
                arrange(theme)

    community_order <- edges %>%
                    group_by(Grupo, Cluster) %>%
                    summarise(total = sum(weight), .groups = "drop") %>%
                    arrange(factor(Grupo, levels = GROUPS), desc(total)) %>%
                    mutate(community = paste0(Grupo, "_", Cluster))

    plot_df <- edges %>%
            mutate(community = factor(paste0(Grupo, "_", Cluster), levels = community_order$community),
                   theme_label = factor(theme_label, levels = rev(theme_levels$theme_label)),
                   Grupo = factor(Grupo, levels = GROUPS)) %>%
            dplyr::select(Grupo, community, theme_label, weight) %>%
            complete(nesting(Grupo, community), theme_label) %>%
            mutate(dark = !is.na(weight) & weight >= 10)

    p <- ggplot(plot_df, aes(x = community, y = theme_label, fill = weight)) +
        geom_tile(color = "grey85", linewidth = 0.3) +
        geom_text(data = dplyr::filter(plot_df, !is.na(weight)),
                  aes(label = weight, color = dark), size = 2.8) +
        scale_fill_gradientn(colours = PUBU_COLORS[3:9], na.value = PUBU_COLORS[1], name = "GO terms") +
        scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
        scale_x_discrete(labels = function(x) sub(".*_", "C", x), expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        facet_grid(. ~ Grupo, scales = "free_x", space = "free_x") +
        labs(x = "Community", y = NULL) +
        theme_minimal(base_size = 11) +
        theme(panel.grid = element_blank(),
              panel.border = element_rect(fill = NA, color = "grey40", linewidth = 0.6),
              panel.spacing = unit(0.8, "lines"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              strip.text = element_text(face = "bold", size = 12))

    return(p)
}

main <- function(){
    d <- godata(annoDb = "org.Hs.eg.db", ont = "BP")

    all_conditions <- get_all_conditions_csv()

    result <- get_go_themes(all_conditions, d, k = 20)

    themes <- result$themes

    check <- check_theme_assignment(result$sim, themes)
    cat("Terms closer to another theme:", nrow(check), "\n")
    print(check)

    write.csv(check, file.path(OUTDATADIR, "go_themes_reassigned.csv"), row.names = FALSE)

    themes <- reassign_themes(themes, check)
    themes <- get_theme_labels(themes, result$sim)

    write.csv(themes, file.path(OUTDATADIR, "go_themes.csv"), row.names = FALSE)
    saveRDS(result$hc, file.path(OUTDATADIR, "go_themes_hclust.rds"))

    edges <- get_theme_edges(all_conditions, themes)
    write.csv(edges, file.path(OUTDATADIR, "community_theme_edges.csv"), row.names = FALSE)

    bipartites <- list()
    metrics <- data.frame()

    for (group in GROUPS){
        b <- build_bipartite(edges, group)
        bipartites[[group]] <- b
        metrics <- bind_rows(metrics, bipartite_metrics(b, group))
        write_graph(b, file.path(OUTDATADIR, paste0("bipartite_", group, ".graphml")), format = "graphml")
    }

    saveRDS(bipartites, file.path(OUTDATADIR, "bipartites.rds"))
    write.csv(metrics, file.path(OUTDATADIR, "bipartite_metrics.csv"), row.names = FALSE)

    theme_table <- theme_community_table(edges)
    write.csv(theme_table, file.path(OUTDATADIR, "theme_communities_by_group.csv"), row.names = FALSE)

    print(metrics)
    print(theme_table)

    p <- plot_bipartites(bipartites, themes)

    ggsave(file.path(OUTFIGDIR, "bipartite_communities_go.png"), p, width = 16, height = 13, dpi = 300, bg = "white")
    ggsave(file.path(OUTFIGDIR, "bipartite_community_theme.pdf"), plot = p, width = 16, height = 12)

    p_heatmap <- plot_theme_heatmap(edges, themes)

    ggsave(file.path(OUTFIGDIR, "heatmap_community_theme.png"), plot = p_heatmap, width = 16, height = 8, dpi = 300, bg = "white")
    ggsave(file.path(OUTFIGDIR, "heatmap_community_theme.pdf"), plot = p_heatmap, width = 16, height = 8)

  for (group in GROUPS){
        comunidades <- read.csv(file.path(COMUNIDADESDIR, paste0("infomap_comunidades_", group, ".csv")))
        print(filter(comunidades, gen %in% c("PTPRD", "IL1RAPL1", "ALCAM", "SPOCK1")) %>% mutate(Grupo = group))
    }
}

main()