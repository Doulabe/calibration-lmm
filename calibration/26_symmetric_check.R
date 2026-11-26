# ==============================================================================
# 26_symmetric_check.R
# Numerical support (calibration of the symmetric treatments).
#   A. Exact identities on random instances of s ~ N(xi beta, diag(1 + delta l_i)):
#      (i)  LR_m = 2{max ltilde_m - max l0_m} = min_delta {D_m(delta) + T_m(delta)},
#           T_m(delta) = m log{1 + t(delta)^2/(nu-1)}, for m = n (ML, as in the pilot) and m = nu;
#           Breusch ordering T_m(delta_tilde_1) >= LR_m >= T_m(delta_hat_0);
#      (ii) log m1(delta) = ltilde_nu(delta) - log(1+v)/2 - E(delta)/2, with
#           0 <= E <= Q^2 / {(1+v)(1 - Q^2/nu)}  (conjugate BF, sigma^2 integrated).
#   B. On six pilot configurations: Delta(u) = 2 log MBF - LR as a function of the effect u,
#      for the marginal BF with fixed tau^2 = 0.04 and with tau ~ half-Cauchy(0, 0.2).
#      The Laplace argument predicts Delta(u) nearly constant in u when the prior scale is
#      large compared with the effect, and a fall at large u for fixed tau^2.
#   C. First-order predictions at threshold c (rho(c) = kappa_j +
#      c^2/(2 nu_j)) against the exact size ratios of 06_treatments_pilot.R (24 configurations).
# Usage (from the repository root): nice -n 15 Rscript calibration/26_symmetric_check.R [cores=3]
# Outputs: results/symmetric_check_{identities,laplace,firstorder}.csv
# ==============================================================================
args <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
B_DRAW <- 200L
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
RES <- "calibration/results"

## ---------------------------------------------------------------- A. exact identities
set.seed(26)
idrows <- list()
for (rep in 1:4) {
  nu <- sample(20:80, 1); n <- nu + 1L
  ell <- c(rexp(nu - 6, 1 / 4), rep(0, 6)); xi <- rnorm(nu) * sqrt(ell + 0.3)
  d0 <- runif(1, 0.05, 2); w0 <- 1 / (1 + d0 * ell); beta <- c(0, 0.4, 1.5)[(rep - 1) %% 3 + 1]
  s <- xi * beta + rnorm(nu) / sqrt(w0)
  L  <- function(d) sum(log1p(d * ell)); Sf <- function(d) sum(xi^2 / (1 + d * ell))
  Af <- function(d) sum(s^2 / (1 + d * ell)); Nf <- function(d) sum(xi * s / (1 + d * ell))
  RSSf <- function(d) Af(d) - Nf(d)^2 / Sf(d)
  t2f <- function(d) (nu - 1) * (Af(d) - RSSf(d)) / RSSf(d); Q2f <- function(d) nu * (1 - RSSf(d) / Af(d))
  argmax <- function(f) { o <- optimize(function(z) f(exp(z)), c(-16, 9.3), maximum = TRUE, tol = 1e-12)
    if (f(0) >= o$objective) c(d = 0, v = f(0)) else c(d = exp(o$maximum), v = o$objective) }
  for (m in c(n, nu)) {
    lt <- function(d) -0.5 * (L(d) + m * log(RSSf(d))); l0 <- function(d) -0.5 * (L(d) + m * log(Af(d)))
    Tm <- function(d) m * log1p(t2f(d) / (nu - 1))
    M1 <- argmax(lt); M0 <- argmax(l0); LR <- 2 * (M1[["v"]] - M0[["v"]])
    mn <- -argmax(function(d) -(2 * (M1[["v"]] - lt(d)) + Tm(d)))[["v"]]
    idrows[[length(idrows) + 1]] <- data.frame(instance = rep, m = if (m == n) "n (ML)" else "nu",
      beta = beta, LR = LR, min_D_plus_T = mn, abs_diff = abs(LR - mn),
      wald_side = Tm(M1[["d"]]), score_side = Tm(M0[["d"]]),
      breusch_ok = Tm(M1[["d"]]) >= LR - 1e-8 && LR >= Tm(M0[["d"]]) - 1e-8)
  }
  ## (ii) marginal-likelihood decomposition, conjugate BF
  e2 <- 0; ebound <- TRUE
  for (tau2 in c(0.04, 1, 25)) for (d in c(0, 0.01, 0.3, 2, 30)) {
    v <- tau2 * Sf(d); Q2 <- Q2f(d)
    lbf <- -0.5 * log1p(v) - (nu / 2) * log1p(-v * Q2 / (nu * (1 + v)))
    lm1 <- -0.5 * (L(d) + nu * log(Af(d))) + lbf
    Tn <- -nu * log1p(-Q2 / nu); Tv <- -nu * log1p(-v * Q2 / (nu * (1 + v))); E <- Tn - Tv
    e2 <- max(e2, abs(lm1 - (-0.5 * (L(d) + nu * log(RSSf(d))) - 0.5 * log1p(v) - 0.5 * E)))
    ebound <- ebound && E >= -1e-12 && E <= Q2 / ((1 + v) * (1 - Q2 / nu)) + 1e-10
  }
  idrows[[length(idrows)]]$mbf_decomp_err <- e2; idrows[[length(idrows)]]$E_bound_ok <- ebound
}
ID <- do.call(rbind, lapply(idrows, function(z) { for (k in c("mbf_decomp_err", "E_bound_ok"))
  if (is.null(z[[k]])) z[[k]] <- NA; z }))
write.csv(ID, file.path(RES, "symmetric_check_identities.csv"), row.names = FALSE)
cat("A. identities\n"); print(ID, row.names = FALSE, digits = 4)

## ---------------------------------------------------------------- shared pilot machinery
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))              # ML LRT grid (as in 06)
QG <- 10^seq(-4, 3, length.out = 61); TAU2 <- 0.04          # MBF quadrature (as in 06)
lwq <- { tq <- log(QG); dq <- tq[2] - tq[1]
  log(c(dq / 2, rep(dq, length(tq) - 2), dq / 2)) + 0.5 * tq - log(pi) - log1p(QG) }
TAU0 <- 0.2; T2 <- 10^seq(-6, 4, length.out = 41)
lw_t2 <- { t <- log(T2); dt <- t[2] - t[1]
  w <- log(c(dt / 2, rep(dt, length(t) - 2), dt / 2)) - log(pi * TAU0) - log1p(T2 / TAU0^2) + 0.5 * t
  w - (max(w) + log(sum(exp(w - max(w))))) }
rowLSE <- function(M) { mx <- apply(M, 1, max); mx + log(rowSums(exp(M - mx))) }
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L, 213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, j = match(g$lead_pos, pos))
}
make_sp <- function(z) {
  if (z$kind == "syn") { set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    list(sp = spectral_setup(Z, j = 25L), id = sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)) } else {
    L <- load_geu(z$i); list(sp = spectral_setup(L$Z, j = L$j), id = sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)) }
}

## ---------------------------------------------------------------- B. Laplace: 2 log MBF - LR along u
UG <- seq(0, 30, by = 0.5)
laplace_config <- function(z) {
  M <- make_sp(z); sp <- M$sp; delta0 <- z$h2 / (1 - z$h2)
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L; B <- B_DRAW
  set.seed(5100 + round(100 * z$h2) + nchar(M$id))
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0; Rk <- s - outer(xi, b); Rp <- spp - xp * b
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  WQ <- 1 / (1 + outer(ell, QG)); ldQ <- colSums(log(1 + outer(ell, QG))); SQg <- colSums(xi^2 * WQ) + xp^2
  vQ <- TAU2 * SQg
  Df <- Dh <- LRm <- matrix(NA_real_, B, length(UG))
  for (g in seq_along(UG)) {
    bv <- UG[g] / sqrt(S0); sk <- Rk + xi * bv; s2p <- Rp + xp * bv; sk2 <- sk * sk
    A2 <- crossprod(sk2, W2) + (s2p^2 + ssr); B2 <- crossprod(xi * sk, W2) + xp * s2p
    RSS2 <- A2 - sweep(B2^2, 2, S2g, "/")
    LRm[, g] <- apply(-sweep(n * log(RSS2), 2, ld2, "+"), 1, max) - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max)
    AQ <- crossprod(sk2, WQ) + (s2p^2 + ssr); BQ <- crossprod(xi * sk, WQ) + xp * s2p
    q <- sweep(BQ^2 / AQ, 2, SQg, "/")                         # Q^2 / nu
    l0 <- -0.5 * sweep(nu * log(AQ), 2, ldQ, "+")
    lw <- sweep(l0, 2, lwq, "+")
    lbf_f <- sweep(-(nu / 2) * log1p(-sweep(q, 2, vQ / (1 + vQ), "*")), 2, 0.5 * log1p(vQ), "-")
    lbf_h <- matrix(-Inf, B, length(QG))
    for (k in seq_along(QG)) {
      vv <- T2 * SQg[k]
      Lk <- outer(q[, k], vv / (1 + vv)); Lk <- -(nu / 2) * log1p(-Lk)
      Lk <- sweep(Lk, 2, -0.5 * log1p(vv) + lw_t2, "+")
      lbf_h[, k] <- rowLSE(Lk)
    }
    den <- rowLSE(lw)
    Df[, g] <- 2 * (rowLSE(lw + lbf_f) - den) - LRm[, g]
    Dh[, g] <- 2 * (rowLSE(lw + lbf_h) - den) - LRm[, g]
  }
  rng <- function(D, umax) { k <- UG <= umax; apply(D[, k, drop = FALSE], 1, function(v) max(v) - min(v)) }
  data.frame(config = M$id, nu_delta = indices(sp, delta0)[["nu_delta"]],
    prior = c("tau2 fixe 0.04", "tau demi-Cauchy"),
    delta_u0_med = c(median(Df[, 1]), median(Dh[, 1])),
    range_u8_med = c(median(rng(Df, 8)), median(rng(Dh, 8))),
    range_u8_q90 = c(quantile(rng(Df, 8), 0.9), quantile(rng(Dh, 8), 0.9)),
    range_u30_med = c(median(rng(Df, 30)), median(rng(Dh, 30))),
    lr_range_u8_med = median(rng(LRm, 8)))
}
cfgB <- list(list(kind = "geu", i = 1L, h2 = 0.05), list(kind = "geu", i = 2L, h2 = 0.05),
             list(kind = "geu", i = 8L, h2 = 0.05), list(kind = "geu", i = 1L, h2 = 0.20),
             list(kind = "syn", n = 200L, nb = 5L, h2 = 0.05), list(kind = "syn", n = 1000L, nb = 20L, h2 = 0.20))
LB <- mclapply(cfgB, function(z) tryCatch(laplace_config(z), error = function(e) conditionMessage(e)),
               mc.cores = NCORES, mc.preschedule = FALSE)
bad <- vapply(LB, is.character, TRUE); if (any(bad)) print(unlist(LB[bad]))
LBd <- do.call(rbind, LB[!bad])
write.csv(LBd, file.path(RES, "symmetric_check_laplace.csv"), row.names = FALSE)
cat("\nB. 2 log MBF - LR along u\n"); print(LBd, row.names = FALSE, digits = 3)

## ---------------------------------------------------------------- C. first-order predictions vs exact
T6 <- read.csv(file.path(RES, "treatments_pilot.csv"))
cfgC <- list()
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20)) cfgC[[length(cfgC) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2)
for (i in seq_len(nrow(GENES))) for (h2 in c(0.05, 0.20)) cfgC[[length(cfgC) + 1]] <- list(kind = "geu", i = i, h2 = h2)
AL <- c(1e-3, 1e-5, 5e-8)
FO <- do.call(rbind, lapply(cfgC, function(z) {
  M <- make_sp(z); ix <- indices(M$sp, z$h2 / (1 - z$h2)); nu <- M$sp$nu0
  do.call(rbind, lapply(AL, function(a) {
    cc <- qt(1 - a / 2, nu - 1); rho <- ix[["kappa"]] + cc^2 / (2 * ix[["nu_delta"]])
    ex <- function(m) T6$size_ratio[T6$config == M$id & T6$method == m & abs(T6$alpha - a) < a * 1e-6]
    data.frame(config = M$id, nu_delta = ix[["nu_delta"]], kappa = ix[["kappa"]], alpha = a, rho_c = rho,
      wald_pred = 2 * pt(-cc / sqrt(1 + rho), nu - 1) / a, wald_exact = ex("wald_h1"),
      score_pred = 2 * pt(-cc * sqrt(1 + rho), nu - 1) / a, score_exact = ex("score_h0"),
      lrt_pred = 1, lrt_exact = ex("lrt_ml"))
  }))
}))
write.csv(FO, file.path(RES, "symmetric_check_firstorder.csv"), row.names = FALSE)
cat("\nC. first-order predictions vs exact (median over configurations)\n")
print(aggregate(cbind(wald_pred, wald_exact, score_pred, score_exact, lrt_exact) ~ alpha, FO, median), digits = 3)
cat("log10 ratio pred/exact, Wald: median", round(median(log10(FO$wald_pred / FO$wald_exact)), 2),
    " IQR", round(quantile(log10(FO$wald_pred / FO$wald_exact), c(.25, .75)), 2), "\n")
