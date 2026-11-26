# ==============================================================================
# 62_protected_wald_mixture.R
# Sizes of the protected Wald test at its maximizing ratio for the two designs whose evaluation
# in 52_protected_design.R exceeded the level by more than two standard errors, recomputed with
# the mixture sampler of 61_importance_sampling_pilot.R (equal-weight mixture of the laws of the
# error contrasts at the ratios delta x {1, 1/3, 1/10, 0, 3, 10}, weights bounded by six) and with
# plain draws on fresh seeds. The protected critical values are those frozen by 52.
# Usage (from the repository root): Rscript calibration/62_protected_wald_mixture.R [cores=4] [draws=2000000]
# Output: results/protected_wald_mixture.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R"); a0 <- commandArgs(TRUE); commandArgs <- function(...) character(0)
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # GENES, load_geu, spectral_setup, pick_grid, GRID2
suppressPackageStartupMessages(library(parallel))
NC <- if (length(a0) >= 1) as.integer(a0[1]) else 4L; B <- if (length(a0) >= 2) as.integer(as.numeric(a0[2])) else 2000000L
E <- read.csv("calibration/results/protected_design_eval.csv", stringsAsFactors = FALSE)
X <- E[E$kind == "argmax" & E$method == "wald_protected" & E$size_ratio - 2 * E$size_se > 1, ]
stopifnot(nrow(X) == 2)
MULT <- c(1, 1/3, 1/10, 0, 3, 10)
wald_size <- function(sp, delta0, c0, mult, seed) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; dm <- delta0 * mult; M <- length(dm)
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2; CH <- 2000L
  v <- unlist(mclapply(seq_len(B %/% CH), function(q) { set.seed(seed + q)
    comp <- sample.int(M, CH, replace = TRUE)
    s <- matrix(rnorm(r * CH), r, CH) * sqrt(1 + outer(ell, dm[comp]))
    spp <- if (sp$nperp) rnorm(CH) else numeric(CH); ssr <- if (sp$nrest > 0) rchisq(CH, sp$nrest) else numeric(CH)
    w <- if (M == 1) rep(1, CH) else { lf <- sapply(dm, function(d) -0.5 * (sum(log1p(d * ell)) + colSums(s^2 / (1 + d * ell))))
      lf0 <- -0.5 * (sum(log1p(delta0 * ell)) + colSums(s^2 / d0)); mx <- apply(lf, 1, max); exp(lf0 - mx) / rowMeans(exp(lf - mx)) }
    b <- (colSums((xi / d0) * s) + xp * spp) / S0; Rk <- s - outer(xi, b); Rp <- spp - xp * b; Rk2 <- Rk * Rk
    A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
    dh1 <- pick_grid(-0.5 * sweep((nu - 1) * log(A_r - sweep(B_r^2, 2, S2g, "/")), 2, ld2 + log(S2g), "+"), GRID2)
    iD <- 1 / (1 + outer(ell, dh1)); A1 <- colSums(Rk2 * iD) + Rp^2 + ssr; B1 <- colSums(xi * Rk * iD) + xp * Rp; S21 <- colSums(xi^2 * iD) + xp^2
    sig <- sqrt((A1 - B1^2 / S21) / (nu - 1)); wa <- sqrt(S21 / S0) / sig; wg <- (B1 / S21) * sqrt(S21) / sig
    w * (pnorm(-(c0 - wg) / wa) + pnorm((-c0 - wg) / wa)) }, mc.cores = NC))
  c(size = mean(v), se = sd(v) / sqrt(length(v)), ess = sum(v)^2 / sum(v^2), draws = length(v))
}
O <- do.call(rbind, lapply(seq_len(nrow(X)), function(k) { x <- X[k, ]; i <- match(x$gene, GENES$symbol); W <- load_geu(i)
  sp <- spectral_setup(W$Z, j = match(x$pos, W$pos))
  do.call(rbind, lapply(c("plain", "mixture"), function(sn) { z <- wald_size(sp, x$ratio, x$threshold, if (sn == "plain") 1 else MULT, 620000L + 100000L * k + 50000L * (sn == "mixture"))
    data.frame(gene = x$gene, pos = x$pos, nu_delta = x$nu_delta, ratio = x$ratio, alpha = x$alpha, threshold = x$threshold,
               size_ratio_52 = x$size_ratio, se_52 = x$size_se, ess_52 = x$size_ess, draws_52 = x$draws,
               scheme = sn, size_ratio = z[["size"]] / x$alpha, se = z[["se"]] / x$alpha, ess = z[["ess"]], draws = z[["draws"]]) })) }))
O$exceeds <- O$size_ratio - 2 * O$se > 1
write.csv(O, "calibration/results/protected_wald_mixture.csv", row.names = FALSE)
options(width = 200); print(O, digits = 4, row.names = FALSE)
