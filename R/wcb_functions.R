# =============================================================
# wcb_functions.R — Restricted wild cluster bootstrap (Webb six-point weights)
# Linear restrictions R b = r; single (t / Wald) and joint (Wald) tests;
# confidence intervals by test inversion. Bank fixed effects removed by
# within-bank demeaning (nested in the clusters). R's RNG: set.seed() fixes draws.
# =============================================================

webb <- c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2))

Rm <- function(cn, rows) {
  R <- matrix(0, length(rows), length(cn), dimnames = list(NULL, cn))
  for (i in seq_along(rows)) R[i, names(rows[[i]])] <- rows[[i]]
  R
}

wcb_setup <- function(m, data, B, seed, cluster = "bank_id", depvar = "ln_mtb") {
  X  <- model.matrix(m, type = "rhs"); g <- factor(data[[cluster]])
  dm <- function(z) z - ave(z, g)
  Xd <- apply(X, 2, dm); yd <- dm(data[[depvar]])
  XtXi <- solve(crossprod(Xd)); bh <- drop(XtXi %*% crossprod(Xd, yd))
  
  stopifnot(max(abs(bh - coef(m)[colnames(Xd)])) < 1e-8)
  
  if (!is.null(seed)) set.seed(seed)
  
  V <- matrix(
    sample(webb, nlevels(g) * B, replace = TRUE),
    nlevels(g), B
  )
  
  list(
    Xd = Xd,
    yd = yd,
    idx = split(seq_len(nrow(Xd)), g),
    XtXi = XtXi,
    bh = bh,
    Vn = V[as.integer(g), , drop = FALSE],
    cn = colnames(Xd)
  )
}

wcb_core <- function(S, R, r) {
  q <- nrow(R)
  A <- R %*% S$XtXi
  
  wald <- function(beta, U) {
    cv <- R %*% beta - r
    
    a <- lapply(
      S$idx,
      function(ix)
        A %*% crossprod(
          S$Xd[ix, , drop = FALSE],
          U[ix, , drop = FALSE]
        )
    )
    
    if (q == 1)
      return(as.vector(
        cv^2 / Reduce(`+`, lapply(a, function(z) z^2))
      ))
    
    v11 <- Reduce(`+`, lapply(a, function(z) z[1, ]^2))
    v22 <- Reduce(`+`, lapply(a, function(z) z[2, ]^2))
    v12 <- Reduce(`+`, lapply(a, function(z) z[1, ] * z[2, ]))
    
    (
      v22 * cv[1, ]^2 -
        2 * v12 * cv[1, ] * cv[2, ] +
        v11 * cv[2, ]^2
    ) /
      (v11 * v22 - v12^2) / q
  }
  
  bh <- matrix(S$bh)
  
  W0 <- wald(
    bh,
    S$yd - S$Xd %*% bh
  )
  
  bt <- bh -
    S$XtXi %*%
    t(R) %*%
    solve(R %*% S$XtXi %*% t(R)) %*%
    (R %*% bh - r)
  
  fit <- drop(S$Xd %*% bt)
  ut <- S$yd - fit
  
  Ys <- fit + ut * S$Vn
  
  Bs <- S$XtXi %*%
    crossprod(S$Xd, Ys)
  
  mean(
    wald(
      Bs,
      Ys - S$Xd %*% Bs
    ) >= W0
  )
}

wcb_test <- function(
    m,
    data,
    rows,
    B = 9999,
    seed = 20261004,
    ci = FALSE,
    level = 0.95,
    npts = 401,
    depvar = "ln_mtb"
) {
  
  S <- wcb_setup(
    m,
    data,
    B,
    seed,
    depvar = depvar
  )
  
  R <- Rm(S$cn, rows)
  
  out <- c(
    p = wcb_core(
      S,
      R,
      rep(0, nrow(R))
    ),
    ci_lo = NA,
    ci_hi = NA
  )
  
  if (ci && nrow(R) == 1) {
    
    est <- drop(R %*% S$bh)
    
    se <- sqrt(
      drop(
        R %*%
          vcov(m)[S$cn, S$cn] %*%
          t(R)
      )
    )
    
    grid <- seq(
      est - 10 * se,
      est + 10 * se,
      length.out = npts
    )
    
    pg <- sapply(
      grid,
      function(b0)
        wcb_core(S, R, b0)
    )
    
    acc <- grid[pg > 1 - level]
    
    out["ci_lo"] <-
      if (pg[1] > 1 - level) -Inf else min(acc)
    
    out["ci_hi"] <-
      if (pg[npts] > 1 - level) Inf else max(acc)
  }
  
  out
}
