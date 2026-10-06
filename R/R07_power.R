# =============================================================
# R07_power.R — Monte Carlo size and power (final specification)
# =============================================================
options(width = 200)
sink("logs/R07_power.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest); library(fixest)
source("R/wcb_functions.R")

stopifnot(digest(file = "data/derived/panel_analysis.csv", algo = "sha256") ==
            "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
p <- readRDS("data/derived/panel_analysis.rds")
d <- p[p$core, ]; d <- d[order(d$bank_id, d$year), ]
d$bank_id <- factor(d$bank_id); d$year_f <- factor(d$year); G <- nlevels(d$bank_id)

NSIM <- 500; BSIM <- 399
fml <- function(y, x) as.formula(paste(y, "~", paste(x, collapse = " + "), "+ year_f | bank_id"))
clw <- function(m, rows) {
  R <- Rm(names(coef(m)), rows); b <- coef(m); V <- vcov(m); q <- nrow(R); Rb <- R %*% b
  Fs <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb) / q
  pf(Fs, q, G - 1, lower.tail = FALSE)
}

# Data-generating process: estimated model with the tested coefficients set to the null
# (zero for NPL; their common mean for equality tests) plus lambda x the estimated departure.
# Errors: normal with bank-specific SD of the estimated residuals.
cells <- list(
  list(test = "NPL (Eq. 3)", x = c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull"), cmp = "Lnpl",
       rows = list(c(Lnpl = 1)), lams = c(0, 1), null = "zero"),
  list(test = "H2a", x = c("cw_car", "cw_negnpl", "Lfull"), cmp = c("cw_car", "cw_negnpl"),
       rows = list(c(cw_car = 1, cw_negnpl = -1)), lams = c(0, 1, 2), null = "equal"),
  list(test = "H2b", x = c("cw_roa", "cw_roe", "cw_eps", "Lfull"), cmp = c("cw_roa", "cw_roe", "cw_eps"),
       rows = list(c(cw_roa = 1, cw_roe = -1), c(cw_roa = 1, cw_eps = -1)), lams = c(0, 1, 2), null = "equal"),
  list(test = "H2c", x = c("cwbsi", "cwfpi", "cwfull"), cmp = c("cwbsi", "cwfpi", "cwfull"),
       rows = list(c(cwbsi = 1, cwfpi = -1), c(cwbsi = 1, cwfull = -1)), lams = c(0, 1, 2), null = "equal"))

set.seed(20261004)
out <- list()
for (cc in cells) {
  m0 <- feols(fml("ln_mtb", cc$x), data = d, cluster = ~bank_id)
  b  <- coef(m0)[cc$cmp]; Xc <- as.matrix(d[, cc$cmp])
  bbar <- if (cc$null == "zero") 0 else mean(b)
  base <- fitted(m0) - drop(Xc %*% b) + bbar * rowSums(Xc)
  dev  <- drop(Xc %*% (b - bbar))
  sdb  <- ave(resid(m0), d$bank_id, FUN = sd)
  for (lam in cc$lams) {
    pw <- pc <- numeric(NSIM)
    for (r in seq_len(NSIM)) {
      dd <- d; dd$ystar <- base + lam * dev + sdb * rnorm(nrow(d))
      m <- feols(fml("ystar", cc$x), data = dd, cluster = ~bank_id)
      pw[r] <- wcb_test(m, dd, cc$rows, B = BSIM, seed = NULL, depvar = "ystar")["p"]
      pc[r] <- clw(m, cc$rows)
    }
    out[[length(out) + 1]] <- data.frame(Test = cc$test, Lambda = lam,
                                         Reject_WCB_05 = mean(pw < .05), Reject_cluster_05 = mean(pc < .05),
                                         Reject_WCB_10 = mean(pw < .10), Reject_cluster_10 = mean(pc < .10),
                                         MC_SE_05 = sqrt(mean(pw < .05) * (1 - mean(pw < .05)) / NSIM), Sims = NSIM)
    cat(sprintf("%-12s lambda = %d: WCB rejects %.3f, clustered rejects %.3f (5%% level)\n",
                cc$test, lam, mean(pw < .05), mean(pc < .05)))
  }
}
pw_tab <- do.call(rbind, out)
cat("\nLambda 0 = size (null true); 1 = effect as estimated; 2 = twice the estimate\n")
print(format(pw_tab, digits = 3), row.names = FALSE)
write.csv(pw_tab, "output/tables/T9_power.csv", row.names = FALSE)

cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()