# =============================================================
# R04_H1.R — H1: valuation relevance of lagged signals (Eqs. 2a-2c, 3)
# =============================================================
options(width = 200)
sink("logs/R04_H1.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest); library(fixest)
source("R/wcb_functions.R")

stopifnot(digest(file = "data/derived/panel_analysis.csv", algo = "sha256") ==
            "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
p <- readRDS("data/derived/panel_analysis.rds")
d <- p[p$core, ]; d <- d[order(d$bank_id, d$year), ]
d$bank_id <- factor(d$bank_id); d$year_f <- factor(d$year)
B <- 9999; SEED <- 20261004

specs <- list(
  Eq2a = c("Lcar", "Lnpl", "Lroa", "Lfull"),
  Eq2b = c("Lcar", "Lnpl", "Lroe", "Lfull"),
  Eq2c = c("Lcar", "Lnpl", "Leps", "Lfull"),
  Eq3  = c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull"))

rows <- list(); fits <- list()
for (s in names(specs)) {
  f <- as.formula(paste("ln_mtb ~", paste(specs[[s]], collapse = " + "), "+ year_f | bank_id"))
  m <- feols(f, data = d, cluster = ~bank_id); fits[[s]] <- m
  cat("\n================ ", s, " ================\n"); print(summary(m))
  ct <- coeftable(m)
  for (v in specs[[s]]) {
    w <- wcb_test(m, d, list(setNames(1, v)), B = B, seed = SEED, ci = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      Model = s, Variable = v, Coef = ct[v, 1], SE_cluster = ct[v, 2], t_cluster = ct[v, 3],
      p_cluster = ct[v, 4], p_wcb = w["p"], CI_lo = w["ci_lo"], CI_hi = w["ci_hi"],
      N = nobs(m), Within_R2 = unname(r2(m, "wr2")), row.names = NULL)
  }
}

# Determinism check: repeating a bootstrap test must give an identical p-value
p1 <- wcb_test(fits$Eq3, d, list(c(Lnpl = 1)), B = B, seed = SEED)["p"]
p2 <- wcb_test(fits$Eq3, d, list(c(Lnpl = 1)), B = B, seed = SEED)["p"]
cat(sprintf("\nDeterminism check (Eq3, NPL): %.6f vs %.6f\n", p1, p2)); stopifnot(identical(p1, p2))

res <- do.call(rbind, rows)
cat("\n================ H1 results (bootstrap p = principal inference) ================\n")
print(format(res, digits = 4), row.names = FALSE)
write.csv(res, "output/tables/T3_H1_baseline.csv", row.names = FALSE)

cmp <- with(res, tapply(sprintf("%.4f [%.3f]", Coef, p_wcb), list(Variable, Model), identity))
cmp <- cmp[c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull"), ]
cat("\nCoefficient [bootstrap p]\n"); print(noquote(cmp))
cat("\nWithin R2:\n"); print(noquote(sapply(fits, function(m) sprintf("%.3f", r2(m, "wr2")))))

cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()