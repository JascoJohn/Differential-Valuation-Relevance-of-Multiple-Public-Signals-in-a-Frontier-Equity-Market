# =============================================================
# R03_descriptives_diagnostics.R — Descriptives and diagnostics
# =============================================================
options(width = 200)
sink("logs/R03_descriptives_diagnostics.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest); library(fixest); library(plm)

stopifnot(digest(file = "data/derived/panel_analysis.csv", algo = "sha256") ==
            "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
p <- readRDS("data/derived/panel_analysis.rds")
d <- p[p$core, ]; d <- d[order(d$bank_id, d$year), ]

# A. Descriptives
vars <- c("mtb", "ln_mtb", "Lcar", "Lnpl", "Lldr", "Lroa", "Lroe", "Leps", "Lfull", "Lsize")
desc <- t(sapply(vars, function(v) { x <- d[[v]]
c(N = sum(!is.na(x)), Mean = mean(x), SD = sd(x), Min = min(x), Median = median(x), Max = max(x)) }))
cat("\nDescriptive statistics, core sample\n"); print(round(desc, 3))
write.csv(round(desc, 6), "output/tables/T1_descriptives.csv")

xts <- t(sapply(vars[-1], function(v) { x <- d[[v]]; bm <- tapply(x, d$bank_id, mean)
c(Overall_SD = sd(x), Between_SD = sd(bm), Within_SD = sd(x - bm[d$bank_id] + mean(x))) }))
cat("\nBetween / within variation\n"); print(round(xts, 4))
write.csv(round(xts, 6), "output/tables/T1b_between_within.csv")

cv <- vars[-1]
cat("\nCorrelations (pooled)\n"); print(round(cor(d[, cv]), 3))
dm <- sapply(cv, function(v) d[[v]] - ave(d[[v]], d$bank_id))
cat("\nCorrelations (within-bank)\n"); print(round(cor(dm), 3))
write.csv(round(cor(dm), 6), "output/tables/T1c_within_correlations.csv")

shares <- sapply(c("car", "npl", "roa", "roe", "eps"), function(v) {
  sq <- p[[paste0("d_", v)]]^2; tapply(sq, p$bank_id, sum, na.rm = TRUE) / sum(sq, na.rm = TRUE) })
cat("\nWithin-bank variance shares by bank\n"); print(round(shares, 3))
write.csv(round(shares, 6), "output/tables/T1d_variance_shares.csv")

# B. Multicollinearity (Eq. 3 regressors, with bank and year dummies)
xv <- c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull")
vif <- sapply(xv, function(v) {
  f <- as.formula(paste(v, "~", paste(setdiff(xv, v), collapse = " + "), "+ factor(bank_id) + factor(year)"))
  1 / (1 - summary(lm(f, data = d))$r.squared) })
cat("\nVariance inflation factors, Eq. 3\n"); print(round(vif, 2))

# C. Residual diagnostics on Eq. 3 (two-way fixed effects)
m <- feols(as.formula(paste("ln_mtb ~", paste(xv, collapse = " + "), "| bank_id + year")), data = d)
e <- resid(m)

# Modified Wald test for groupwise heteroskedasticity (Greene 2000)
s2i <- tapply(e^2, d$bank_id, mean); Ti <- tapply(e, d$bank_id, length)
Vi  <- tapply(e^2, d$bank_id, function(z) sum((z - mean(z))^2)) / (Ti * (Ti - 1))
W   <- sum((s2i - mean(e^2))^2 / Vi)

# Pesaran (2004) CD test for cross-sectional dependence (pairwise correlations, unbalanced panel)
M <- tapply(e, list(d$year, d$bank_id), sum); N <- ncol(M); s <- 0; sa <- 0; k <- 0
for (i in 1:(N - 1)) for (j in (i + 1):N) {
  ok <- !is.na(M[, i]) & !is.na(M[, j]); rho <- cor(M[ok, i], M[ok, j])
  s <- s + sqrt(sum(ok)) * rho; sa <- sa + abs(rho); k <- k + 1 }
CD <- sqrt(2 / (N * (N - 1))) * s

# Wooldridge (2002) test for first-order serial correlation (Drukker 2003 implementation)
full <- p[order(p$bank_id, p$year), ]
for (yr in 2018:2025) full[[paste0("yr", yr)]] <- as.numeric(full$year == yr)
vv <- c("ln_mtb", xv, paste0("yr", 2018:2025))
D <- full[, c("bank_id", "year", "core")]
for (v in vv) D[[v]] <- full[[v]] - ave(full[[v]], full$bank_id, FUN = function(z) c(NA, head(z, -1)))
D <- D[D$core & complete.cases(D[, vv]), ]
X <- as.matrix(D[, setdiff(vv, "ln_mtb")]); D$u <- as.vector(D$ln_mtb - X %*% qr.solve(X, D$ln_mtb))
D$ul <- ave(D$u, D$bank_id, FUN = function(z) c(NA, head(z, -1)))
D$yl <- ave(D$year, D$bank_id, FUN = function(z) c(NA, head(z, -1)))
D$ul[!is.na(D$yl) & D$year - D$yl != 1] <- NA
E <- D[!is.na(D$ul), ]; rho1 <- sum(E$ul * E$u) / sum(E$ul^2); r <- E$u - rho1 * E$ul
G <- length(unique(E$bank_id)); n <- nrow(E)
meat <- sum(tapply(E$ul * r, E$bank_id, sum)^2)
V <- (G / (G - 1)) * ((n - 1) / (n - 1)) * meat / sum(E$ul^2)^2
Fw <- (rho1 + 0.5)^2 / V

diag_tab <- data.frame(
  Test = c("Modified Wald (groupwise heteroskedasticity)", "Pesaran CD (cross-sectional dependence)",
           "Wooldridge (first-order serial correlation)"),
  Statistic = c(sprintf("chi2(%d) = %.2f", N, W), sprintf("CD = %.3f; mean |rho| = %.3f", CD, sa / k),
                sprintf("F(1,%d) = %.3f", G - 1, Fw)),
  p = c(pchisq(W, N, lower.tail = FALSE), 2 * pnorm(-abs(CD)), pf(Fw, 1, G - 1, lower.tail = FALSE)))
cat("\nResidual diagnostics, Eq. 3\n"); print(diag_tab, row.names = FALSE)

# D. Hausman test (Amemiya variance components; Swamy-Arora not estimable with five banks)
pd <- pdata.frame(d, index = c("bank_id", "year"))
fh <- as.formula(paste("ln_mtb ~", paste(xv, collapse = " + "), "+ factor(year)"))
hz <- tryCatch(phtest(plm(fh, data = pd, model = "within"),
                      plm(fh, data = pd, model = "random", random.method = "amemiya")),
               error = function(err) { cat("Hausman not computable:", conditionMessage(err), "\n"); NULL })
cat("\nHausman test, fixed vs random effects\n"); if (!is.null(hz)) print(hz)
diag_tab <- rbind(diag_tab, data.frame(Test = "Hausman (FE vs RE, Amemiya)",
                                       Statistic = if (is.null(hz)) "not computable" else sprintf("chi2(%d) = %.2f", hz$parameter, hz$statistic),
                                       p = if (is.null(hz)) NA else hz$p.value))
write.csv(diag_tab, "output/tables/T2_diagnostics.csv", row.names = FALSE)

cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()