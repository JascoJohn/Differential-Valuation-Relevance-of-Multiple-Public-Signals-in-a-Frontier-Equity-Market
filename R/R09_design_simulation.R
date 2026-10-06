# =============================================================
# R09_design_simulation.R — Exploratory design simulation (not part of the official pipeline)
# Power of H2a-H2c tests with 5-25 banks under two scenarios for how new banks vary.
# =============================================================
options(width = 200)
dir.create("output/exploratory", showWarnings = FALSE)
sink("logs/R09_design_simulation.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest); library(fixest); library(MASS)
source("R/wcb_functions.R")

stopifnot(digest(file = "data/derived/panel_analysis.csv", algo = "sha256") ==
            "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
p <- readRDS("data/derived/panel_analysis.rds")
d <- p[p$core, ]; d <- d[order(d$bank_id, d$year), ]
d$bank_id <- factor(d$bank_id); d$year_f <- factor(d$year)

NSIM <- 300; BSIM <- 399; GS <- c(5, 10, 15, 20, 25); LAMS <- c(0, 1)
fml <- function(y, x) as.formula(paste(y, "~", paste(x, collapse = " + "), "+ year_f | bank_id"))
clw <- function(m, rows, G) {
  R <- Rm(names(coef(m)), rows); b <- coef(m); V <- vcov(m); q <- nrow(R); Rb <- R %*% b
  pf(as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb) / q, q, G - 1, lower.tail = FALSE)
}

tests <- list(
  H2a = list(x = c("cw_car", "cw_negnpl"), rows = list(c(cw_car = 1, cw_negnpl = -1)), ctrl = "Lfull"),
  H2b = list(x = c("cw_roa", "cw_roe", "cw_eps"),
             rows = list(c(cw_roa = 1, cw_roe = -1), c(cw_roa = 1, cw_eps = -1)), ctrl = "Lfull"),
  H2c = list(x = c("cwbsi", "cwfpi", "cwfull"),
             rows = list(c(cwbsi = 1, cwfpi = -1), c(cwbsi = 1, cwfull = -1)), ctrl = NULL))

# DGP parameters from the final models
par <- lapply(tests, function(tt) {
  m <- feols(fml("ln_mtb", c(tt$x, tt$ctrl)), data = d, cluster = ~bank_id)
  b <- coef(m)[tt$x]
  list(b = b, bbar = mean(b), sdb = sapply(split(resid(m), d$bank_id), sd))
})
cat("\nDGP coefficients (per common within-bank SD):\n"); print(lapply(par, `[[`, "b"))

sigs <- c("cw_car", "cw_negnpl", "cw_roa", "cw_roe", "cw_eps", "cw_full")
Rcor <- cor(d[, sigs])                       # within-bank correlations (signals are bank-demeaned)
cat("\nWithin-bank signal correlations used in the Even scenario:\n"); print(round(Rcor, 3))

addcomp <- function(s) {
  b <- (s$cw_car + s$cw_negnpl) / 2; f <- (s$cw_roa + s$cw_roe + s$cw_eps) / 3
  s$cwbsi <- b / sd(b); s$cwfpi <- f / sd(f); s$cwfull <- s$cw_full; s
}
# Resampled: the 5 actual banks plus (G - 5) banks drawn from their observed trajectories
gen_resampled <- function(G) {
  src <- c(levels(d$bank_id), sample(levels(d$bank_id), G - 5, replace = TRUE))
  s <- do.call(rbind, lapply(seq_along(src), function(k) {
    x <- d[d$bank_id == src[k], c("year", sigs, "cwbsi", "cwfpi", "cwfull")]
    x$src <- src[k]; x$bank_id <- k; x }))
  s$bank_id <- factor(s$bank_id); s$year_f <- factor(s$year); s
}
# Even: every bank has its own variation (9 years), observed within-bank correlations
gen_even <- function(G, Tn = 9) {
  s <- data.frame(bank_id = factor(rep(seq_len(G), each = Tn)), year = rep(2016 + seq_len(Tn), G))
  Z <- mvrnorm(nrow(s), rep(0, length(sigs)), Rcor); colnames(Z) <- sigs
  for (v in sigs) { z <- Z[, v] - ave(Z[, v], s$bank_id); s[[v]] <- z / sd(z) }
  s$year_f <- factor(s$year); addcomp(s)
}

set.seed(20261004)
out <- list()
for (sc in c("Resampled", "Even")) for (G in GS) for (h in names(tests)) for (lam in LAMS) {
  tt <- tests[[h]]; pr <- par[[h]]; pw <- pc <- numeric(NSIM)
  for (r in seq_len(NSIM)) {
    s <- if (sc == "Resampled") gen_resampled(G) else gen_even(G)
    X <- as.matrix(s[, tt$x])
    mu <- drop(X %*% (pr$bbar + lam * (pr$b - pr$bbar)))
    sde <- if (sc == "Resampled") pr$sdb[s$src] else sample(pr$sdb, G, replace = TRUE)[as.integer(s$bank_id)]
    s$y <- as.numeric(mu + as.numeric(sde) * rnorm(nrow(s)))
    m <- feols(fml("y", tt$x), data = s, cluster = ~bank_id)
    pw[r] <- wcb_test(m, s, tt$rows, B = BSIM, seed = NULL, depvar = "y")["p"]
    pc[r] <- clw(m, tt$rows, G)
  }
  out[[length(out) + 1]] <- data.frame(Scenario = sc, Banks = G, Test = h, Lambda = lam,
                                       Reject_WCB_05 = mean(pw < .05), Reject_cluster_05 = mean(pc < .05),
                                       MC_SE = sqrt(mean(pw < .05) * (1 - mean(pw < .05)) / NSIM), Sims = NSIM)
  cat(sprintf("%-9s G=%2d %s lambda=%d: WCB %.3f | clustered %.3f   (%s)\n",
              sc, G, h, lam, mean(pw < .05), mean(pc < .05), format(Sys.time(), "%H:%M:%S")))
}
res <- do.call(rbind, out)
write.csv(res, "output/exploratory/design_simulation.csv", row.names = FALSE)

cat("\n=== Size (lambda = 0): WCB rejection rate at 5% ===\n")
print(xtabs(Reject_WCB_05 ~ Scenario + Test + Banks, data = res[res$Lambda == 0, ]))
cat("\n=== Power (lambda = 1, differences as estimated): WCB rejection rate at 5% ===\n")
print(xtabs(Reject_WCB_05 ~ Scenario + Test + Banks, data = res[res$Lambda == 1, ]))

cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()