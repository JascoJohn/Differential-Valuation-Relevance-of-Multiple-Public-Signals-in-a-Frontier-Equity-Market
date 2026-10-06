# =============================================================
# R05_H2.R — H2: differential valuation relevance (Eqs. 4-6)
# =============================================================
options(width = 200)
sink("logs/R05_H2.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest); library(fixest)
source("R/wcb_functions.R")

stopifnot(digest(file = "data/derived/panel_analysis.csv", algo = "sha256") ==
            "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
p <- readRDS("data/derived/panel_analysis.rds")
d <- p[p$core, ]; d <- d[order(d$bank_id, d$year), ]
d$bank_id <- factor(d$bank_id); d$year_f <- factor(d$year)
B <- 9999; SEED <- 20261004; G <- nlevels(d$bank_id)

clw <- function(m, rows) {                         # clustered Wald F(q, G-1), secondary
  R <- Rm(names(coef(m)), rows); b <- coef(m); V <- vcov(m); q <- nrow(R); Rb <- R %*% b
  Fs <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb) / q
  c(F = Fs, p = pf(Fs, q, G - 1, lower.tail = FALSE), Rb = if (q == 1) drop(Rb) else NA)
}

models <- list(
  Eq4 = c("cw_car", "cw_negnpl", "Lfull"),
  Eq5 = c("cw_roa", "cw_roe", "cw_eps", "Lfull"),
  Eq6 = c("cwbsi", "cwfpi", "cwfull"))
fits <- list(); coefrows <- list()
for (s in names(models)) {
  f <- as.formula(paste("ln_mtb ~", paste(models[[s]], collapse = " + "), "+ year_f | bank_id"))
  m <- feols(f, data = d, cluster = ~bank_id); fits[[s]] <- m
  cat("\n================ ", s, " ================\n"); print(summary(m))
  ct <- coeftable(m)
  for (v in models[[s]]) {
    w <- wcb_test(m, d, list(setNames(1, v)), B = B, seed = SEED, ci = TRUE)
    coefrows[[length(coefrows) + 1]] <- data.frame(Model = s, Variable = v, Coef = ct[v, 1],
                                                   SE_cluster = ct[v, 2], t_cluster = ct[v, 3], p_cluster = ct[v, 4], p_wcb = w["p"],
                                                   CI_lo = w["ci_lo"], CI_hi = w["ci_hi"], N = nobs(m), Within_R2 = unname(r2(m, "wr2")), row.names = NULL)
  }
}
coefs <- do.call(rbind, coefrows)
cat("\n================ H2 model coefficients (per common within-bank SD) ================\n")
print(format(coefs, digits = 4), row.names = FALSE)
write.csv(coefs, "output/tables/T4a_H2_coefficients.csv", row.names = FALSE)

tests <- list(
  list("Eq4", "H2a", "CAR = -NPL",          list(c(cw_car = 1, cw_negnpl = -1)), "principal"),
  list("Eq5", "H2b", "ROA = ROE = EPS",     list(c(cw_roa = 1, cw_roe = -1), c(cw_roa = 1, cw_eps = -1)), "principal (joint)"),
  list("Eq5", "H2b", "ROA = ROE",           list(c(cw_roa = 1, cw_roe = -1)), "pairwise"),
  list("Eq5", "H2b", "ROA = EPS",           list(c(cw_roa = 1, cw_eps = -1)), "pairwise"),
  list("Eq5", "H2b", "ROE = EPS",           list(c(cw_roe = 1, cw_eps = -1)), "pairwise"),
  list("Eq6", "H2c", "BSI = FPI = FULLREP", list(c(cwbsi = 1, cwfpi = -1), c(cwbsi = 1, cwfull = -1)), "principal (joint)"),
  list("Eq6", "H2c", "BSI = FPI",           list(c(cwbsi = 1, cwfpi = -1)), "pairwise"),
  list("Eq6", "H2c", "BSI = FULLREP",       list(c(cwbsi = 1, cwfull = -1)), "pairwise"),
  list("Eq6", "H2c", "FPI = FULLREP",       list(c(cwfpi = 1, cwfull = -1)), "pairwise"))

trows <- list()
for (tst in tests) {
  m <- fits[[tst[[1]]]]; rr <- tst[[4]]; q <- length(rr)
  w  <- wcb_test(m, d, rr, B = B, seed = SEED, ci = (q == 1))
  cl <- clw(m, rr)
  trows[[length(trows) + 1]] <- data.frame(Model = tst[[1]], Hypothesis = tst[[2]], Restriction = tst[[3]],
                                           Role = tst[[5]], q = q, Estimate = cl["Rb"],
                                           Statistic = if (q == 1) sign(cl["Rb"]) * sqrt(cl["F"]) else cl["F"],
                                           Stat_type = if (q == 1) "t (clustered)" else "F (clustered)",
                                           p_wcb = w["p"], CI_lo = w["ci_lo"], CI_hi = w["ci_hi"], p_cluster = cl["p"], row.names = NULL)
}
tt <- do.call(rbind, trows)
for (h in c("H2b", "H2c")) {
  jp <- tt$p_wcb[tt$Hypothesis == h & tt$Role == "principal (joint)"]
  tt$Role[tt$Hypothesis == h & tt$Role == "pairwise"] <-
    if (jp < 0.05) "pairwise (joint rejected)" else "pairwise (exploratory: joint not rejected)"
}

# Determinism check: repeated tests give identical p-values
a1 <- wcb_test(fits$Eq4, d, list(c(cw_car = 1, cw_negnpl = -1)), B = B, seed = SEED)["p"]
a2 <- wcb_test(fits$Eq4, d, list(c(cw_car = 1, cw_negnpl = -1)), B = B, seed = SEED)["p"]
j1 <- wcb_test(fits$Eq6, d, list(c(cwbsi = 1, cwfpi = -1), c(cwbsi = 1, cwfull = -1)), B = B, seed = SEED)["p"]
j2 <- wcb_test(fits$Eq6, d, list(c(cwbsi = 1, cwfpi = -1), c(cwbsi = 1, cwfull = -1)), B = B, seed = SEED)["p"]
cat(sprintf("\nDeterminism check: H2a %.6f vs %.6f; H2c joint %.6f vs %.6f\n", a1, a2, j1, j2))
stopifnot(identical(a1, a2), identical(j1, j2))

cat("\n================ H2 equality tests (bootstrap p = principal inference) ================\n")
print(format(tt, digits = 4), row.names = FALSE)
write.csv(tt, "output/tables/T4b_H2_tests.csv", row.names = FALSE)

cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()