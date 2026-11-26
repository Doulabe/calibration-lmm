# ==============================================================================
# 47_protected_plugin_pilot.R
#
# Feasibility pilot for new rows of Table 2 of the article: the protected and the
# plug-in tests built on the Wald and likelihood-ratio statistics.
#   protected  reject when T >= c* = sup_delta c_alpha(delta), where c_alpha(delta) is
#              the exact (1 - alpha) null quantile of T at true ratio delta. Then
#              p_sup(T) <= alpha exactly when T >= c*, so the size is <= alpha.
#   plug-in    reject when T >= c_alpha(delta_hat(r)), delta_hat being the estimate with
#              the candidate (REML for Wald, ML for LR), a function of r only.
# Exact quantiles use B draws of the error contrasts at each delta of a grid, with common
# random numbers across the grid. The LR is tabulated on a grid of u: at each delta of the
# estimation grid the null residual sum of squares is quadratic in u, so no cross-product
# is recomputed per u. The calibration draws at delta0 replicate 36_calibration_enlarged.R
# exactly, so the nominal Wald and LR sizes are checked against results/calib_enlarged.
# Usage (from the repository root): nice -n 10 Rscript calibration/47_protected_plugin_pilot.R
# Output: calibration/results/protected_plugin_pilot.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
B_DRAW <- 1000L
source("calibration/00_functions.R")
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))          # estimation grid of 36
AL    <- c(1e-3, 1e-5, 5e-8)
DG    <- c(0, 10^seq(-4, 4, by = 0.5), 1e6, 1e10)       # true ratios for the supremum
UG    <- seq(-20, 20, by = 0.02)                        # u grid for the LR

pick_grid <- function(O, grid) {                          # as in 36
  g <- max.col(O, ties.method = "first"); dh <- grid[g]
  int <- g > 2 & g < length(grid)
  if (any(int)) {
    ii <- which(int); lg <- log(grid)
    x1 <- lg[g[ii] - 1]; x2 <- lg[g[ii]]; x3 <- lg[g[ii] + 1]
    f1 <- O[cbind(ii, g[ii] - 1)]; f2 <- O[cbind(ii, g[ii])]; f3 <- O[cbind(ii, g[ii] + 1)]
    den <- (x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1)
    num <- (x2 - x1)^2 * (f2 - f3) - (x2 - x3)^2 * (f2 - f1)
    xs <- ifelse(abs(den) > 1e-14, x2 - 0.5 * num / den, x2)
    dh[ii] <- exp(pmin(pmax(xs, x1), x3))
  }
  dh
}
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {                                 # as in 36
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, pos = pos)
}

## error contrasts at true ratio delta from standard normals Zk, zp and ssr (as calib_config of 36)
contrasts <- function(sp, delta, Zk, zp, ssr) {
  d <- 1 + delta * sp$ell; S0 <- sum(sp$xi^2 / d) + sp$xp^2
  s <- Zk * sqrt(d)
  b <- (colSums((sp$xi / d) * s) + sp$xp * zp) / S0
  list(Rk = s - outer(sp$xi, b), Rp = zp - sp$xp * b, ssr = ssr, S0 = S0)
}

## Wald coefficients (t = wa u + wg given r) and LR tabulated on UG
pieces <- function(sp, C) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- nu + 1L; B <- ncol(C$Rk)
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  Rk2 <- C$Rk * C$Rk
  A_r <- crossprod(Rk2, W2) + (C$Rp^2 + C$ssr); B_r <- crossprod(xi * C$Rk, W2) + xp * C$Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  dh1 <- pick_grid(-0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"), GRID2)
  L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+")
  ml1 <- L1ml[cbind(seq_len(B), max.col(L1ml, "first"))]; dh1ml <- pick_grid(L1ml, GRID2)
  iD <- 1 / (1 + outer(ell, dh1))
  A <- colSums(Rk2 * iD) + C$Rp^2 + C$ssr; Bx <- colSums(xi * C$Rk * iD) + xp * C$Rp
  S2 <- colSums(xi^2 * iD) + xp^2; sig_h <- sqrt((A - Bx^2 / S2) / (nu - 1))
  S2m <- matrix(S2g, B, length(GRID2), byrow = TRUE); ldm <- matrix(ld2, B, length(GRID2), byrow = TRUE)
  LRt <- matrix(0, B, length(UG))
  for (k in seq_along(UG)) {
    bv <- UG[k] / sqrt(C$S0)
    O <- -n * log(A_r + (2 * bv) * B_r + (bv * bv) * S2m) - ldm
    LRt[, k] <- ml1 - O[cbind(seq_len(B), max.col(O, "first"))]
  }
  list(dh1 = dh1, dh1ml = dh1ml, wa = sqrt(S2 / C$S0) / sig_h, wg = (Bx / S2) * sqrt(S2) / sig_h, LRt = LRt)
}

## rejection probabilities with u ~ N(lam, 1); c is a scalar or one threshold per draw
tail_w <- function(c, lam, P) mean(pnorm(lam - (c - P$wg) / P$wa) + pnorm((-c - P$wg) / P$wa - lam))
tail_lr <- function(c, lam, LRt) {
  K <- length(UG); z <- UG - lam
  PU <- pnorm(z, lower.tail = FALSE); PL <- pnorm(z)
  cell <- ifelse(z[-K] >= 0, PU[-K] - PU[-1], PL[-1] - PL[-K])
  G <- LRt >= c
  L1 <- G[, -K, drop = FALSE]; L2 <- G[, -1, drop = FALSE]
  tot <- drop((L1 & L2) %*% cell) + G[, 1] * PL[1] + G[, K] * PU[K]
  tr <- which(xor(L1, L2), arr.ind = TRUE)
  if (nrow(tr)) {
    b <- tr[, 1]; k <- tr[, 2]
    v1 <- LRt[cbind(b, k)]; v2 <- LRt[cbind(b, k + 1)]
    cb <- if (length(c) == 1) rep(c, length(b)) else c[b]
    zs <- UG[k] + (cb - v1) / (v2 - v1) * (UG[k + 1] - UG[k]) - lam
    m <- ifelse(L1[tr],
                ifelse(z[k] >= 0, PU[k] - pnorm(zs, lower.tail = FALSE), pnorm(zs) - PL[k]),
                ifelse(zs >= 0, pnorm(zs, lower.tail = FALSE) - PU[k + 1], PL[k + 1] - pnorm(zs)))
    rs <- rowsum(m, b); tot[as.integer(rownames(rs))] <- tot[as.integer(rownames(rs))] + rs[, 1]
  }
  mean(tot)
}
qsolve <- function(f, alpha, hi, it = 40L) {
  if (f(hi) > alpha) return(NA_real_)
  lo <- 0
  for (i in seq_len(it)) { mid <- (lo + hi) / 2; if (f(mid) > alpha) lo <- mid else hi <- mid }
  (lo + hi) / 2
}
interp_c <- function(Q, dh) {                             # c_alpha(delta_hat), log scale in delta
  sapply(seq_len(ncol(Q)), function(k) {
    out <- rep(Q[1, k], length(dh)); pos <- dh > 0
    out[pos] <- approx(log10(DG[-1]), Q[-1, k], log10(dh[pos]), rule = 2)$y
    out })
}

run_cfg <- function(sp, delta0, seed_calib, id) {
  B <- B_DRAW; nu <- sp$nu0
  tq <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  t0 <- Sys.time()
  set.seed(seed_calib + 1L)                               # common random numbers on the delta grid
  Zk <- matrix(rnorm(sp$r * B), sp$r, B); zp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  QW <- QL <- matrix(NA_real_, length(DG), length(AL)); tp <- tq_ <- 0
  for (g in seq_along(DG)) {
    ta <- Sys.time(); P <- pieces(sp, contrasts(sp, DG[g], Zk, zp, ssr))
    tp <- tp + as.numeric(difftime(Sys.time(), ta, units = "secs")); ta <- Sys.time()
    for (k in seq_along(AL)) {
      QW[g, k] <- qsolve(function(c) tail_w(c, 0, P), AL[k], hi = 1e3)
      QL[g, k] <- qsolve(function(c) tail_lr(c, 0, P$LRt), AL[k], hi = 400)
    }
    tq_ <- tq_ + as.numeric(difftime(Sys.time(), ta, units = "secs"))
  }
  set.seed(seed_calib)                                    # calibration draws of 36
  Zk0 <- matrix(rnorm(sp$r * B), sp$r, B); zp0 <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr0 <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  P0 <- pieces(sp, contrasts(sp, delta0, Zk0, zp0, ssr0))
  cW <- apply(QW, 2, max); cL <- apply(QL, 2, max)
  pW <- interp_c(QW, P0$dh1); pL <- interp_c(QL, P0$dh1ml)
  rows <- list()
  for (k in seq_along(AL)) {
    add <- function(m, thr, sz, pw) rows[[length(rows) + 1]] <<- data.frame(config = id, method = m, alpha = AL[k],
      size_ratio = sz / AL[k], power = pw, threshold = thr)
    add("wald_nominal", tq[k], tail_w(tq[k], 0, P0), tail_w(tq[k], tq[k], P0))
    add("lr_nominal", cq[k], tail_lr(cq[k], 0, P0$LRt), tail_lr(cq[k], tq[k], P0$LRt))
    add("wald_protected", cW[k], tail_w(cW[k], 0, P0), tail_w(cW[k], tq[k], P0))
    add("lr_protected", cL[k], tail_lr(cL[k], 0, P0$LRt), tail_lr(cL[k], tq[k], P0$LRt))
    add("wald_plugin", median(pW[, k]), tail_w(pW[, k], 0, P0), tail_w(pW[, k], tq[k], P0))
    add("lr_plugin", median(pL[, k]), tail_lr(pL[, k], 0, P0$LRt), tail_lr(pL[, k], tq[k], P0$LRt))
  }
  out <- do.call(rbind, rows)
  cat(sprintf("\n%s  nu_delta = %.1f  grid: pieces %.0f s, quantiles %.0f s, total %.0f s\n", id,
              indices(sp, delta0)[["nu_delta"]], tp, tq_, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  cat("argmax delta of c_alpha (Wald):", DG[apply(QW, 2, which.max)], "  (LR):", DG[apply(QL, 2, which.max)], "\n")
  print(round(cbind(delta = DG, QW, QL), 3))
  attr(out, "QW") <- QW; attr(out, "QL") <- QL
  out
}

CF <- list(
  list(id = "geu_ERAP2_96252589_h0.05", kind = "geu", i = 1L, pos = 96252589L, h2 = 0.05),
  list(id = "syn_n1000_b5_h0.05", kind = "syn", n = 1000L, nb = 5L, h2 = 0.05))
res <- list()
for (z in CF) {
  if (z$kind == "geu") {
    W <- load_geu(z$i); sp <- spectral_setup(W$Z, j = match(z$pos, W$pos))
    seed <- 9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973)
  } else {
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    seed <- 4242 + 1000 * z$n + z$nb + round(100 * z$h2)
  }
  out <- run_cfg(sp, z$h2 / (1 - z$h2), seed, z$id)
  ref <- readRDS(file.path("calibration/results/calib_enlarged", paste0(z$id, ".rds")))
  ref <- ref[ref$method %in% c("wald_t", "lrt_chi2"), c("method", "alpha", "size_ratio", "power")]
  ref$method <- ifelse(ref$method == "wald_t", "wald_nominal", "lr_nominal")
  names(ref)[3:4] <- c("size_ratio_36", "power_36")
  out <- merge(out, ref, all.x = TRUE, sort = FALSE)
  print(out[order(out$alpha, out$method), ], digits = 3, row.names = FALSE)
  res[[length(res) + 1]] <- out
}
write.csv(do.call(rbind, res), "calibration/results/protected_plugin_pilot.csv", row.names = FALSE)
