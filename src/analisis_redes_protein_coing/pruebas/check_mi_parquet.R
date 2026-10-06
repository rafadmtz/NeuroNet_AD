library(arrow)

GROUPS <- c("Not_AD", "Low", "Intermediate", "High")

BASE      <- "/export/space3/users/silvanac/NeuroNet_AD"
GENEINDIR <- file.path(BASE, "output/protein_coding/mi_matrices")   # AJUSTA si tus parquet están en otra carpeta

for (group in GROUPS) {
  path <- file.path(GENEINDIR, paste0("Vip_", group, "_protein_coding_mi.parquet"))
  cat("\n==", group, "==\n")
  cat("archivo:", path, "\n")
  cat("existe:", file.exists(path), "\n")

  ds <- open_dataset(path)
  columnas <- ds$schema$names

  cat("número de columnas:", length(columnas), "\n")
  cat("número de filas:", nrow(ds), "\n")
  cat("primeras 10 columnas:\n")
  print(head(columnas, 10))
  cat("últimas 5 columnas:\n")
  print(tail(columnas, 5))

  cat("¿genes conocidos entre las columnas? (SUGP2, BIN1, CANX):\n")
  print(c("SUGP2", "BIN1", "CANX") %in% columnas)

  cat("genes mitocondriales (formato MT- o MT.):\n")
  print(grep("^MT[-.]", columnas, value = TRUE))
}