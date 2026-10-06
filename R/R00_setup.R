# =============================================================
# R00_setup.R — Environment and reproducibility setup (FINAL, R)
# Run from the project root:
#   source("R/R00_setup.R", echo = TRUE, max.deparse.length = Inf)
# =============================================================

# --- Folder structure ---
dirs <- c("data/raw", "data/derived", "R", "logs", "output/tables")
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- Log (output written to console and file) ---
sink("logs/R00_setup.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")

# --- Packages (installed into the project library) ---
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(pkgType = "binary",
        install.packages.compile.from.source = "never")
pkgs <- c("readxl", "fixest", "clubSandwich", "sandwich",
          "lmtest", "plm", "digest")
install.packages(pkgs)

# Wild cluster bootstrap: CRAN first, then the author's r-universe repository
ok <- tryCatch({ install.packages("fwildclusterboot"); TRUE },
               warning = function(w) FALSE, error = function(e) FALSE)
if (!ok || !requireNamespace("fwildclusterboot", quietly = TRUE)) {
  install.packages("fwildclusterboot",
                   repos = c(SA = "https://s3alfisc.r-universe.dev",
                             CRAN = "https://cloud.r-project.org"))
}

# --- Record exact versions ---
for (p in c(pkgs, "fwildclusterboot", "renv")) {
  cat(sprintf("%-18s %s\n", p, as.character(packageVersion(p))))
}
cat("RNG kind:", paste(RNGkind(), collapse = " / "), "\n")

# --- Pin all installed versions in renv.lock ---
renv::snapshot(type = "all", prompt = FALSE)

print(sessionInfo())
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()