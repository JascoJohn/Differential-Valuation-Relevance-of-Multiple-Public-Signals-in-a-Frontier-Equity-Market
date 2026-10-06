# =============================================================
# R06_robustness.R — Robustness and sensitivity (final specification)
# =============================================================
options(width = 200)
sink("logs/R06_robustness.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest); library(fixest); library(clubSandwich)
source("R/wcb_functions.R")

stopifnot(digest(file = "data/derived/panel_analysis.csv", algo = "sha256") ==
            "de5b2268e868832710b01c711ae081dc3248d933feab3fef31f9c224f40ec1f6")
p <- readRDS("data/derived/panel_analysis.rds")
d <- p[p$core, ]; d <- d[order(d$bank_id, d$year), ]
B <- 9999; SEED <- 20261004
banknames <- c(B01 = "CRDB", B02 = "NMB", B03 = "DCB", B05 = "Mkombozi", B06 = "Mwalimu")

# --- Helpers ---
buildcw <- function(dd) {                                   # rebuild v2 scaling on a (sub)sample
  for (v in c("car", "npl", "ldr", "roa", "roe", "eps", "full")) {
    x <- dd[[paste0("L", v)]]; dv <- x - ave(x, dd$bank_id); dd[[paste0("cw_", v)]] <- dv / sd(dv) }
  dd$cw_negnpl <- -dd$cw_npl
  b <- (dd$cw_car + dd$cw_negnpl) / 2; f <- (dd$cw_roa + dd$cw_roe + dd$cw_eps) / 3
  dd$cwbsi <- b / sd(b); dd$cwfpi <- f / sd(f); dd$cwfull <- dd$cw_full
  dd
}
prep <- function(dd) { dd$bank_id <- droplevels(factor(dd$bank_id)); dd$year_f <- droplevels(factor(dd$year)); dd }
fitm <- function(dd, x) feols(as.formula(paste("ln_mtb ~", paste(x, collapse = " + "), "+ year_f | bank_id")),
                              data = dd, cluster = ~bank_id)
rowc <- function(block, model, m, dd, term)
  data.frame(Block = block, Model = model, Term = term, Coef = unname(coef(m)[term]),
             p_wcb = unname(wcb_test(m, dd, list(setNames(1, term)), B = B, seed = SEED)["p"]),
             N = nobs(m), row.names = NULL)
rowt <- function(block, model, m, dd, label, rows)
  data.frame(Block = block, Model = model, Term = label, Coef = NA,
             p_wcb = unname(wcb_test(m, dd, rows, B = B, seed = SEED)["p"]), N = nobs(m), row.names = NULL)
H2a <- list(c(cw_car = 1, cw_negnpl = -1))
H2b <- list(c(cw_roa = 1, cw_roe = -1), c(cw_roa = 1, cw_eps = -1))
H2c <- list(c(cwbsi = 1, cwfpi = -1), c(cwbsi = 1, cwfull = -1))

runset <- function(dd, block, ctrl = character(0)) {        # Eq. 3 and Eqs. 4-6 with H2 tests
  dd <- prep(dd); out <- list()
  m3 <- fitm(dd, c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull", ctrl))
  for (v in c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull")) out[[length(out) + 1]] <- rowc(block, "Eq3", m3, dd, v)
  m4 <- fitm(dd, c("cw_car", "cw_negnpl", "Lfull", ctrl))
  for (v in c("cw_car", "cw_negnpl")) out[[length(out) + 1]] <- rowc(block, "Eq4", m4, dd, v)
  out[[length(out) + 1]] <- rowt(block, "Eq4", m4, dd, "H2a: CAR = -NPL", H2a)
  m5 <- fitm(dd, c("cw_roa", "cw_roe", "cw_eps", "Lfull", ctrl))
  for (v in c("cw_roa", "cw_roe", "cw_eps")) out[[length(out) + 1]] <- rowc(block, "Eq5", m5, dd, v)
  out[[length(out) + 1]] <- rowt(block, "Eq5", m5, dd, "H2b joint", H2b)
  m6 <- fitm(dd, c("cwbsi", "cwfpi", "cwfull", ctrl))
  for (v in c("cwbsi", "cwfpi", "cwfull")) out[[length(out) + 1]] <- rowc(block, "Eq6", m6, dd, v)
  out[[length(out) + 1]] <- rowt(block, "Eq6", m6, dd, "H2c joint", H2c)
  do.call(rbind, out)
}
show <- function(x) print(format(x, digits = 4), row.names = FALSE)

# Consistency: rebuilding the scaling on the full core sample reproduces the stored variables
dchk <- buildcw(d)
stopifnot(max(abs(dchk$cw_car - d$cw_car)) < 1e-12, max(abs(dchk$cwbsi - d$cwbsi)) < 1e-12)
cat("Scaling rebuild reproduces stored variables.\n")

res <- list()
cat("\n######## Main specification (reference) ########\n");  res$main  <- runset(d, "Main");                   show(res$main)
cat("\n######## R1: Size added ########\n");                  res$size  <- runset(d, "Size added", "Lsize");   show(res$size)

cat("\n######## R2: LDR added ########\n")
dl <- prep(d); out <- list()
m <- fitm(dl, c("Lcar", "Lnpl", "Lldr", "Lroa", "Lroe", "Leps", "Lfull"))
for (v in c("Lcar", "Lnpl", "Lldr", "Lroa", "Lroe", "Leps", "Lfull")) out[[length(out) + 1]] <- rowc("LDR added", "Eq3", m, dl, v)
m <- fitm(dl, c("cw_car", "cw_negnpl", "cw_ldr", "Lfull"))
for (v in c("cw_car", "cw_negnpl", "cw_ldr")) out[[length(out) + 1]] <- rowc("LDR added", "Eq4", m, dl, v)
out[[length(out) + 1]] <- rowt("LDR added", "Eq4", m, dl, "H2a: CAR = -NPL", H2a)
res$ldr <- do.call(rbind, out); show(res$ldr)

cat("\n######## R3: Residualized performance ########\n")
dr <- prep(d); out <- list()
for (k in c("roa", "roe", "eps")) {
  o <- setdiff(c("roa", "roe", "eps"), k)
  r <- resid(lm(as.formula(paste0("cw_", k, " ~ cw_", o[1], " + cw_", o[2], " + bank_id")), data = dr))
  dr[[paste0("r_", k)]] <- r / sd(r)
}
for (k in c("r_roa", "r_roe", "r_eps")) { m <- fitm(dr, c(k, "Lfull")); out[[length(out) + 1]] <- rowc("Residualized", paste(k, "alone"), m, dr, k) }
m <- fitm(dr, c("r_roa", "r_roe", "r_eps", "Lfull"))
for (k in c("r_roa", "r_roe", "r_eps")) out[[length(out) + 1]] <- rowc("Residualized", "Joint", m, dr, k)
out[[length(out) + 1]] <- rowt("Residualized", "Joint", m, dr, "rROA = rROE = rEPS",
                               list(c(r_roa = 1, r_roe = -1), c(r_roa = 1, r_eps = -1)))
res$resid <- do.call(rbind, out); show(res$resid)

cat("\n######## R4: Winsorized signals (5/95) ########\n")
dw <- prep(d); m <- fitm(dw, c("Wcar", "Wnpl", "Wroa", "Wroe", "Weps", "Lfull"))
res$wins <- do.call(rbind, lapply(c("Wcar", "Wnpl", "Wroa", "Wroe", "Weps", "Lfull"), function(v) rowc("Winsorized", "Eq3", m, dw, v)))
show(res$wins)

cat("\n######## R5: Stale-price observations excluded ########\n")
res$stale <- runset(buildcw(d[!d$stale, ]), "Stale excluded"); show(res$stale)

cat("\n######## R6: Eq. 5 augmented with banking-risk signals ########\n")
da <- prep(d); out <- list()
m <- fitm(da, c("cw_roa", "cw_roe", "cw_eps", "cw_car", "cw_negnpl", "Lfull"))
for (v in c("cw_roa", "cw_roe", "cw_eps", "cw_car", "cw_negnpl")) out[[length(out) + 1]] <- rowc("Eq5 augmented", "Eq5+", m, da, v)
out[[length(out) + 1]] <- rowt("Eq5 augmented", "Eq5+", m, da, "H2b joint", H2b)
out[[length(out) + 1]] <- rowt("Eq5 augmented", "Eq5+", m, da, "H2a: CAR = -NPL", H2a)
res$aug <- do.call(rbind, out); show(res$aug)

cat("\n######## R7: Leave-one-bank-out (scaling rebuilt; 4 clusters) ########\n")
lobo <- list()
for (b in names(banknames)) {
  lab <- paste("Drop", banknames[b]); cat("\n---", lab, "---\n")
  lobo[[b]] <- runset(buildcw(d[d$bank_id != b, ]), lab); show(lobo[[b]])
}
res$lobo <- do.call(rbind, lobo)

cat("\n######## R8: Dominant-bank exclusions (from R7) ########\n")
pick <- function(tab, model, term) tab[tab$Model == model & tab$Term == term, c("Coef", "p_wcb", "N")]
dom <- function(signal, model, term, b, h2) {
  f <- pick(res$main, model, term); x <- pick(lobo[[b]], model, term)
  data.frame(Signal = signal, Excluded = banknames[b], N = x$N, Coef_full = f$Coef, p_full = f$p_wcb,
             Coef_excl = x$Coef, p_excl = x$p_wcb,
             H2_test = h2, H2_p_full = pick(res$main, model, h2)$p_wcb, H2_p_excl = pick(lobo[[b]], model, h2)$p_wcb,
             H2_conclusion_changes = (pick(res$main, model, h2)$p_wcb < .05) != (pick(lobo[[b]], model, h2)$p_wcb < .05),
             row.names = NULL)
}
domtab <- rbind(dom("CAR", "Eq4", "cw_car", "B06", "H2a: CAR = -NPL"),
                dom("ROA", "Eq5", "cw_roa", "B06", "H2b joint"),
                dom("ROE", "Eq5", "cw_roe", "B05", "H2b joint"),
                dom("EPS", "Eq5", "cw_eps", "B02", "H2b joint"))
nplt <- do.call(rbind, lapply(names(banknames), function(b) dom("NPL", "Eq4", "cw_negnpl", b, "H2a: CAR = -NPL")))
cat("\nDominant-bank exclusions:\n"); show(domtab)
cat("\nNPL, each bank excluded:\n"); show(nplt)

cat("\n######## R9: Panel-corrected SEs (Beck-Katz, casewise covariance) ########\n")
pcse <- function(dd, xv, tests) {
  dd <- prep(dd)
  X <- model.matrix(as.formula(paste("~", paste(xv, collapse = " + "), "+ bank_id + year_f")), dd)
  y <- dd$ln_mtb; XtXi <- solve(crossprod(X)); b <- drop(XtXi %*% crossprod(X, y)); e <- drop(y - X %*% b)
  E <- tapply(e, list(dd$year, as.character(dd$bank_id)), sum)
  Ec <- E[complete.cases(E), , drop = FALSE]; Sig <- crossprod(Ec) / nrow(Ec)
  meat <- matrix(0, ncol(X), ncol(X))
  for (t in unique(dd$year)) {
    ix <- which(dd$year == t); bk <- as.character(dd$bank_id[ix])
    meat <- meat + t(X[ix, , drop = FALSE]) %*% Sig[bk, bk, drop = FALSE] %*% X[ix, , drop = FALSE]
  }
  V <- XtXi %*% meat %*% XtXi; se <- sqrt(diag(V))
  coefs <- data.frame(Term = xv, Coef = b[xv], SE_pcse = se[xv], z = b[xv] / se[xv],
                      p = 2 * pnorm(-abs(b[xv] / se[xv])), row.names = NULL)
  tst <- do.call(rbind, lapply(names(tests), function(nm) {
    R <- Rm(colnames(X), tests[[nm]]); Rb <- R %*% b
    W <- as.numeric(t(Rb) %*% solve(R %*% V %*% t(R)) %*% Rb)
    data.frame(Test = nm, chi2 = W, df = nrow(R), p = pchisq(W, nrow(R), lower.tail = FALSE)) }))
  cat("Common years used for the cross-bank covariance:", rownames(Ec), "\n")
  list(coefs = coefs, tests = tst)
}
pc3 <- pcse(d, c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull"), list())
pc4 <- pcse(d, c("cw_car", "cw_negnpl", "Lfull"), list("H2a: CAR = -NPL" = H2a))
pc5 <- pcse(d, c("cw_roa", "cw_roe", "cw_eps", "Lfull"), list("H2b joint" = H2b))
pc6 <- pcse(d, c("cwbsi", "cwfpi", "cwfull"), list("H2c joint" = H2c))
cat("\nEq. 3 (PCSE):\n"); show(pc3$coefs)
pcse_tests <- rbind(pc4$tests, pc5$tests, pc6$tests); cat("\nH2 tests (PCSE):\n"); show(pcse_tests)

cat("\n######## R10: CR2 small-sample inference (Satterthwaite t; HTZ Wald) ########\n")
cr2 <- function(dd, xv, cons) {
  dd <- prep(dd)
  m <- lm(as.formula(paste("ln_mtb ~", paste(xv, collapse = " + "), "+ bank_id + year_f")), data = dd)
  V <- vcovCR(m, cluster = dd$bank_id, type = "CR2")
  ct <- as.data.frame(coef_test(m, vcov = V, test = "Satterthwaite", coefs = xv))
  pc <- grep("^p_", names(ct), value = TRUE)[1]; dc <- grep("^df", names(ct), value = TRUE)[1]
  coefs <- data.frame(Term = xv, Coef = coef(m)[xv], df = ct[[dc]], p_CR2 = ct[[pc]], row.names = NULL)
  tst <- if (is.null(cons)) NULL else {
    w <- as.data.frame(Wald_test(m, constraints = constrain_equal(cons), vcov = V, test = "HTZ"))
    data.frame(Test = paste(cons, collapse = " = "), F = w$Fstat, df_num = w$df_num, df_denom = w$df_denom, p = w$p_val)
  }
  list(coefs = coefs, tests = tst)
}
c3 <- cr2(d, c("Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull"), NULL)
c4 <- cr2(d, c("cw_car", "cw_negnpl", "Lfull"), c("cw_car", "cw_negnpl"))
c5 <- cr2(d, c("cw_roa", "cw_roe", "cw_eps", "Lfull"), c("cw_roa", "cw_roe", "cw_eps"))
c6 <- cr2(d, c("cwbsi", "cwfpi", "cwfull"), c("cwbsi", "cwfpi", "cwfull"))
cat("\nEq. 3 (CR2):\n"); show(c3$coefs)
cr2_tests <- rbind(c4$tests, c5$tests, c6$tests); cat("\nH2 tests (CR2/HTZ):\n"); show(cr2_tests)

# --- Save ---
allres <- do.call(rbind, res[c("main", "size", "ldr", "resid", "wins", "stale", "aug", "lobo")])
write.csv(allres, "output/tables/T5_robustness_bootstrap.csv", row.names = FALSE)
write.csv(rbind(domtab, nplt), "output/tables/T6_dominant_bank_exclusions.csv", row.names = FALSE)
write.csv(rbind(cbind(Spec = "Eq3", pc3$coefs)), "output/tables/T7a_pcse_eq3.csv", row.names = FALSE)
write.csv(pcse_tests, "output/tables/T7b_pcse_H2.csv", row.names = FALSE)
write.csv(c3$coefs, "output/tables/T8a_cr2_eq3.csv", row.names = FALSE)
write.csv(cr2_tests, "output/tables/T8b_cr2_H2.csv", row.names = FALSE)

cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()