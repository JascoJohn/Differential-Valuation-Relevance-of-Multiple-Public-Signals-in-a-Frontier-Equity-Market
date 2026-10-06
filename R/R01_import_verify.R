# =============================================================
# R01_import_verify.R — Import and verify the locked dataset
# =============================================================
options(width = 200)
sink("logs/R01_import_verify.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(readxl); library(digest)

f <- "data/raw/Tanzania_DSE_AnalyticalDataset_v1.3.1_LOCKED.xlsx"
stopifnot(file.exists(f))
cat("File:", f, "\n")
cat("Size (bytes):", file.info(f)$size, "\n")
cat("File SHA-256:", digest(file = f, algo = "sha256"), "\n")
cat("Panel-content checksum recorded in the workbook (VERSION_LOG, compute_panel_checksum.py):\n",
    "  92fc00fdf74a8333a6e823429f44e57af2c6943eb726d543d590df6e50681043\n")
cat("\nSheets:\n"); print(excel_sheets(f))

panel <- as.data.frame(read_excel(f, sheet = "ANALYTICAL_PANEL"))
cat("\nDimensions:", dim(panel), "\n"); str(panel); print(summary(panel))

# Structural checks against the dataset's documented design (5 banks x 2016-2025)
stopifnot(nrow(panel) == 50, ncol(panel) == 22,
          all(table(panel$bank_id) == 10), all(sort(unique(panel$year)) == 2016:2025))
cat("\nStructure verified: 5 banks x 10 years, 22 variables.\n")

panel <- panel[order(panel$bank_id, panel$year), ]
write.csv(panel, "data/derived/panel_imported.csv", row.names = FALSE, na = "")
saveRDS(panel, "data/derived/panel_imported.rds")
cat("Imported panel CSV SHA-256:", digest(file = "data/derived/panel_imported.csv", algo = "sha256"), "\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()