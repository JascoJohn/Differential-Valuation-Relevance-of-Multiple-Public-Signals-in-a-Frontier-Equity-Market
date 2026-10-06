# =============================================================
# R02_construct.R — Analysis variables (final specification)
# =============================================================
options(width = 200)
sink("logs/R02_construct.log", split = TRUE)
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
library(digest)

stopifnot(digest(file = "data/derived/panel_imported.csv", algo = "sha256") ==
            "0e7db0e18ce5c906f7af7eef5e71bccc3b310bae5fc3382170df961bbf35775a")
p <- readRDS("data/derived/panel_imported.rds")
p <- p[order(p$bank_id, p$year), ]

reldif <- function(a, b) abs(a - b) / (abs(b) + 1)
chk <- function(ok, msg) { if (!all(ok)) stop("FAILED: ", msg); cat("PASS:", msg, "\n") }

# 1. Integrity of derived columns
chk(is.na(p$mtb)   | reldif(p$mtb,   p$market_cap / p$book_equity)  < 1e-6, "mtb = market_cap / book_equity")
chk(is.na(p$roa)   | reldif(p$roa,   p$net_profit / p$total_assets) < 1e-6, "roa = net_profit / total_assets")
chk(is.na(p$roe)   | reldif(p$roe,   p$net_profit / p$book_equity)  < 1e-6, "roe = net_profit / book_equity")
chk(is.na(p$ln_ta) | reldif(p$ln_ta, log(p$total_assets))           < 1e-6, "ln_ta = ln(total_assets)")
chk(is.na(p$full_report) | p$full_report %in% c(0, 1), "full_report in {0,1}")

# 2. Units (percent) and dependent variable
for (v in c("car", "npl", "ldr", "roa", "roe")) p[[paste0(v, "_p")]] <- 100 * p[[v]]
p$ln_mtb <- log(p$mtb)

# 3. One-year lags (prior year of the same bank)
lagmap <- c(car = "car_p", npl = "npl_p", ldr = "ldr_p", roa = "roa_p", roe = "roe_p",
            eps = "eps", full = "full_report", size = "ln_ta", mcap = "market_cap")
prev <- p[, c("bank_id", "year", unname(lagmap))]
prev$year <- prev$year + 1
names(prev)[-(1:2)] <- paste0("L", names(lagmap))
p <- merge(p, prev, by = c("bank_id", "year"), all.x = TRUE, sort = FALSE)
p <- p[order(p$bank_id, p$year), ]; rownames(p) <- NULL

# 4. Core estimation sample
p$core <- complete.cases(p[, c("ln_mtb", "Lcar", "Lnpl", "Lroa", "Lroe", "Leps", "Lfull")])
cat("\nCore sample N =", sum(p$core), "\n"); print(table(bank = p$bank_id, core = p$core))
cat("FULLREP = 0 observations in core:", sum(p$Lfull[p$core] == 0), "\n")

# 5. Common within-bank SD standardization (core sample)
sdw <- c()
for (v in c("car", "npl", "ldr", "roa", "roe", "eps", "full")) {
  x <- p[[paste0("L", v)]]; x[!p$core] <- NA
  d <- x - ave(x, p$bank_id, FUN = function(z) mean(z, na.rm = TRUE))
  sdw[v] <- sd(d, na.rm = TRUE)
  p[[paste0("d_", v)]] <- d; p[[paste0("cw_", v)]] <- d / sdw[v]
}
cat("\nCommon within-bank SDs:\n"); print(round(sdw, 4))
p$cw_negnpl <- -p$cw_npl

# 6. Composites on the same scale
p$bsi2 <- (p$cw_car + p$cw_negnpl) / 2
p$fpi2 <- (p$cw_roa + p$cw_roe + p$cw_eps) / 3
p$cwbsi <- p$bsi2 / sd(p$bsi2, na.rm = TRUE)
p$cwfpi <- p$fpi2 / sd(p$fpi2, na.rm = TRUE)
p$cwfull <- p$cw_full
for (v in c("cw_car", "cw_negnpl", "cw_roa", "cw_roe", "cw_eps", "cwbsi", "cwfpi", "cwfull")) {
  bm <- tapply(p[[v]][p$core], p$bank_id[p$core], mean)
  chk(all(abs(bm) < 1e-10) && abs(sd(p[[v]], na.rm = TRUE) - 1) < 1e-10, paste(v, ": bank means 0, SD 1"))
}

# 7. Robustness variables
for (v in c("car", "npl", "roa", "roe", "eps")) {          # winsorized at 5th/95th percentiles (core)
  x <- p[[paste0("L", v)]]; q <- quantile(x[p$core], c(.05, .95), type = 2)
  w <- pmin(pmax(x, q[1]), q[2]); w[!p$core] <- NA; p[[paste0("W", v)]] <- w
  cat(sprintf("W%-4s winsorized observations: %d\n", v, sum(p$core & w != x)))
}
p$stale   <- !is.na(p$market_cap) & !is.na(p$Lmcap) & reldif(p$market_cap, p$Lmcap) < 1e-9
p$core_ns <- p$core & !p$stale
cat("\nStale-price core observations (market cap unchanged from prior year):\n")
print(p[p$stale & p$core, c("bank_id", "year", "market_cap", "mtb")])

# 8. Save
write.csv(p, "data/derived/panel_analysis.csv", row.names = FALSE, na = "")
saveRDS(p, "data/derived/panel_analysis.rds")
cat("Analysis CSV SHA-256:", digest(file = "data/derived/panel_analysis.csv", algo = "sha256"), "\n")
cat("Run finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()