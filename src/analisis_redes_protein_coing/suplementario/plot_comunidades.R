GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

COLORES <- c(Not_AD = "#4A78B0", Low = "#E3C042", Intermediate = "#D4843A", High = "#CD5544")

BASE      <- "/export/space3/users/silvanac/NeuroNet_AD"
INDIR     <- file.path(BASE, "output/protein_coding/data/comunidades_infomap")
OUTFIGDIR <- file.path(BASE, "output/protein_coding/figures/comunidades_infomap")

dir.create(OUTFIGDIR, showWarnings = FALSE, recursive = TRUE)

main <- function(){

    tamanos <- list()

    for (group in GROUPS){
        cat("grupo:", group, "\n")
        comunidades <- read.csv(file.path(INDIR, paste0("infomap_comunidades_", group, ".csv")))
        split_comunidades <- split(comunidades$gen, comunidades$modulo)

        tamanos[[group]] <- sort(lengths(split_comunidades), decreasing = TRUE)
    }

    y_max <- max(unlist(tamanos))

    png(file.path(OUTFIGDIR, "tamano_comunidades_infomap.png"), width = 2400, height = 1800, res = 300)
    par(mfrow = c(2, 2), mar = c(1, 4, 3, 1), mgp = c(2, 0.6, 0), oma = c(2, 0, 3, 0))

    for (group in GROUPS){
        t <- tamanos[[group]]
        pct_top3 <- round(100 * sum(t[1:3]) / sum(t), 1)

        colores <- rep(adjustcolor(COLORES[[group]], alpha.f = 0.35), length(t))
        colores[1:3] <- COLORES[[group]]

        barplot(t, col = colores, border = NA, names.arg = NA, ylim = c(0, y_max),
                main = paste0(group, " (", pct_top3, "% de genes)"),
                ylab = "número de genes")
        abline(h = 5, lty = 2, col = "grey40")

        cat(group, "comunidades:", length(t), "top 3:", pct_top3, "%\n")
    }

    mtext("Top 3 comunidades por red y % de genes que abarcan", side = 3, outer = TRUE, line = 1, font = 2, cex = 1.2)
    mtext("comunidades (ordenadas por tamaño)", side = 1, outer = TRUE, line = 0.5)

    dev.off()
}

main()