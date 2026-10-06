# =============================================================
# master.R — Reproduce the full analysis from the locked dataset
# Run from the project root:  source("R/master.R")
# First run creates reference/; later runs verify against it.
# =============================================================
options(width = 200)
library(digest)
sink("logs/master.log", split = TRUE)
cat("Master run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(R.version.string, "\n")
print(renv::status())

# --- 1. Remove all derived data and outputs ---
unlink(list.files("data/derived", full.names = TRUE))
unlink(list.files("output/tables", full.names = TRUE))
cat("Derived data and output tables removed.\n")

# --- 2. Run the pipeline (each script writes its own log) ---
scripts <- c("R01_import_verify", "R02_construct", "R03_descriptives_diagnostics",
             "R04_H1", "R05_H2", "R06_robustness", "R07_power")
for (s in scripts) {
  cat("\n### Running", s, "at", format(Sys.time(), "%H:%M:%S"), "\n")
  source(file.path("R", paste0(s, ".R")), local = new.env(), echo = TRUE, max.deparse.length = Inf)
}

# --- 3. Hashes of data and outputs ---
tabs  <- sort(list.files("output/tables", pattern = "\\.csv$"))
files <- c(file.path("data/derived", c("panel_imported.csv", "panel_analysis.csv")),
           file.path("output/tables", tabs))
hashes <- data.frame(file = files,
                     sha256 = vapply(files, function(f) digest(file = f, algo = "sha256"), ""),
                     row.names = NULL)
write.csv(hashes, "output/OUTPUT_HASHES.csv", row.names = FALSE)
cat("\nOutput hashes:\n"); print(hashes, row.names = FALSE)

stopifnot(hashes$sha256[1] == "0e7db0e18ce5c906f7af7eef5e71bccc3b310bae5fc3382170df961bbf35775a",
          hashes$sha256[2] == "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
cat("Data hashes verified.\n")

# --- 4. Numerical comparison with the reference outputs ---
if (!dir.exists("reference")) {
  dir.create("reference")
  file.copy(file.path("output/tables", tabs), "reference")
  cat("\nReference outputs created from this run (", length(tabs), "tables).",
      "Run master.R again to verify reproduction.\n")
} else {
  rtabs <- sort(list.files("reference", pattern = "\\.csv$"))
  stopifnot(identical(tabs, rtabs))
  for (f in tabs) {
    a <- read.csv(file.path("output/tables", f)); b <- read.csv(file.path("reference", f))
    stopifnot(identical(dim(a), dim(b)), identical(names(a), names(b)))
    num <- vapply(a, is.numeric, TRUE)
    stopifnot(identical(a[!num], b[!num]), identical(is.na(a[num]), is.na(b[num])))
    dmax <- if (any(num)) max(abs(as.matrix(a[num]) - as.matrix(b[num])), na.rm = TRUE) else 0
    if (!is.finite(dmax)) dmax <- 0
    cat(sprintf("%-40s max abs difference = %.2e  %s\n", f, dmax, if (dmax < 1e-8) "PASS" else "FAIL"))
    stopifnot(dmax < 1e-8)
  }
  cat("\nALL OUTPUTS REPRODUCED.\n")
}

cat("Master run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()