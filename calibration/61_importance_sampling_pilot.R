# ==============================================================================
# 61_importance_sampling_pilot.R
# Importance sampling of the error contrasts in the conditional tail calculation, on four
# configurations of 36_calibration_enlarged.R. The proposal is an equal-weight mixture of the
# Gaussian laws of the contrasts at the ratios delta_0 x {1, 1/3, 1/10, 0, 3, 10}, and each draw is
# weighted by the density at delta_0 over the density of the mixture, a weight bounded by six.
# Sizes of the Wald and likelihood-ratio tests with the plain draws and with the mixture, at equal
# numbers of draws, compared with the evaluation draws of 52_protected_design.R.
# Usage (from the repository root): Rscript calibration/61_importance_sampling_pilot.R [cores=4] [draws Wald=20000] [draws LR=4000]
# Output: results/importance_sampling_pilot.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R"); a0 <- commandArgs(TRUE); commandArgs <- function(...) character(0)
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))
suppressPackageStartupMessages(library(parallel))
NC <- if (length(a0) >= 1) as.integer(a0[1]) else 4L
BW <- if (length(a0) >= 2) as.integer(a0[2]) else 20000L; BL <- if (length(a0) >= 3) as.integer(a0[3]) else 4000L
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)
ids <- vapply(cfgs, cfg_id, "")
## per-draw conditional rejection probabilities of the Wald and LR tests, given the coordinates s of the proposal draws
perdraw <- function(sp, delta0, s, spp, ssr, what) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- nu + 1L; B <- ncol(s)
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp; RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1 <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, GRID2)
  ml1 <- apply(-sweep(n * log(RSS_r), 2, ld2, "+"), 1, max)
  iD <- 1 / (1 + outer(ell, dh1)); A1 <- colSums(Rk2 * iD) + Rp^2 + ssr; B1 <- colSums(xi * Rk * iD) + xp * Rp; S21 <- colSums(xi^2 * iD) + xp^2
  sig_h <- sqrt((A1 - B1^2 / S21) / (nu - 1)); wa <- sqrt(S21 / S0) / sig_h; wg <- (B1 / S21) * sqrt(S21) / sig_h
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1); out <- list(dh1 = dh1)
  if ("wald" %in% what) out$wald <- sapply(seq_along(AL), function(k) pnorm(-(tq_[k] - wg) / wa) + pnorm((-tq_[k] - wg) / wa))
  if ("lrt" %in% what) {
    ev <- function(u, key) { bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
      A2 <- crossprod(sk2, W2) + (s2p^2 + ssr)
      if (key == "lrt") return(ml1 - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max))
      dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2); iD0 <- 1 / (1 + outer(ell, dh0))
      (colSums(xi * sk * iD0) + xp * s2p) }
    lo <- rep(-30, B); hi <- rep(30, B)
    for (it in 1:30) { mid <- (lo + hi) / 2; up <- ev(mid, "Q") < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    split <- (lo + hi) / 2
    out$lrt <- sapply(seq_along(AL), function(k) { c0 <- cq[k]; f <- function(u) ev(u, "lrt")
      always <- f(split) >= c0
      lo <- split; hi <- rep(120, B); nr <- f(hi) < c0
      for (it in 1:26) { mid <- (lo + hi) / 2; up <- f(mid) < c0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      ur <- (lo + hi) / 2; ur[nr] <- Inf
      lo <- rep(-120, B); hi <- split; nl <- f(lo) < c0
      for (it in 1:26) { mid <- (lo + hi) / 2; dn <- f(mid) < c0; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
      ul <- (lo + hi) / 2; ul[nl] <- -Inf; ur[always] <- split[always]; ul[always] <- split[always]
      pnorm(-ur) + pnorm(ul) })
  }
  out
}
## mult: multipliers of delta_0 in the mixture; kap: variance inflations along the direction v below (one here; larger
## values, tried with the mixture, increased the variance)
one <- function(id, B, what, mult, kap, seed) {
  z <- cfgs[[match(id, ids)]]
  if (z$kind == "syn") { set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
  } else { W <- load_geu(z$i); sp <- spectral_setup(W$Z, j = match(z$pos, W$pos)) }
  delta0 <- z$h2 / (1 - z$h2); ell <- sp$ell; r <- sp$r; xi <- sp$xi; xp <- sp$xp
  ## coordinates: r kernel directions and the direction of the candidate outside them (variance 1)
  ellx <- c(ell, 0); xx <- c(xi, xp); d0 <- 1 + delta0 * ellx
  a <- xx / sum(xx^2) - (xx / d0) / sum(xx^2 / d0)          # OLS estimate minus GLS estimate at delta0, as a'c
  v <- sqrt(d0) * a; v <- v / sqrt(sum(v^2))                 # direction in the coordinates whitened at delta0
  G <- expand.grid(m = seq_along(mult), k = seq_along(kap)); dm <- delta0 * mult; NG <- nrow(G)
  CH <- 1000L; nch <- B %/% CH
  R <- mclapply(seq_len(nch), function(q) { set.seed(seed + q)
    comp <- sample.int(NG, CH, replace = TRUE); dd <- dm[G$m[comp]]; kk <- kap[G$k[comp]]
    zz <- matrix(rnorm((r + 1L) * CH), r + 1L, CH); tt <- drop(crossprod(v, zz))
    zz <- zz + outer(v, (kk - 1) * tt)                       # variance kappa^2 along v
    cc <- zz * sqrt(1 + outer(ellx, dd))                     # coordinates of the dataset
    ldens <- function(d, k) { zm <- cc / sqrt(1 + d * ellx); tm <- drop(crossprod(v, zm))
      -0.5 * sum(log1p(d * ellx)) - 0.5 * colSums(zm^2) - log(k) + 0.5 * tm^2 * (1 - 1 / k^2) }
    lf <- sapply(seq_len(NG), function(g) ldens(dm[G$m[g]], kap[G$k[g]]))
    lf0 <- ldens(delta0, 1); mx <- apply(lf, 1, max); w <- exp(lf0 - mx) / rowMeans(exp(lf - mx))
    ssr <- if (sp$nrest > 0) rchisq(CH, sp$nrest) else numeric(CH)
    p <- perdraw(sp, delta0, cc[seq_len(r), , drop = FALSE], if (sp$nperp) cc[r + 1L, ] else numeric(CH), ssr, what)
    list(w = w, dh1 = p$dh1, wald = p$wald, lrt = p$lrt) }, mc.cores = NC)
  w <- unlist(lapply(R, `[[`, "w")); dh1 <- unlist(lapply(R, `[[`, "dh1"))
  res <- do.call(rbind, lapply(what, function(key) { P <- do.call(rbind, lapply(R, `[[`, key))
    do.call(rbind, lapply(seq_along(AL), function(k) { vv <- w * P[, k]; m <- mean(vv)
      data.frame(config = id, method = key, alpha = AL[k], B = length(vv), size_ratio = m / AL[k], rel_se = sd(vv) / sqrt(length(vv)) / m,
                 ess = sum(vv)^2 / sum(vv^2), top1 = max(vv) / sum(vv), share_zero = sum(vv[dh1 < 1e-4]) / sum(vv), wmax = max(w), wmean = mean(w)) })) }))
  res
}
SCH <- list(plain = list(m = 1, k = 1), mixture = list(m = c(1, 1/3, 1/10, 0, 3, 10), k = 1))
WANT <- c("geu_SNHG5_86504582_h0.05", "geu_ERAP2_96252589_h0.05", "geu_ZNF266_9544276_h0.20", "geu_ERAP2_95943006_h0.05")
O <- do.call(rbind, lapply(WANT, function(id) do.call(rbind, lapply(names(SCH), function(sn) {
  x <- rbind(one(id, BW, "wald", SCH[[sn]]$m, SCH[[sn]]$k, 71000), one(id, BL, "lrt", SCH[[sn]]$m, SCH[[sn]]$k, 72000)); x$scheme <- sn; x }))))
E <- read.csv("calibration/results/protected_design_eval.csv")
E <- E[E$kind == "config" & E$method %in% c("wald_t", "lrt_chi2"), ]; E$method <- ifelse(E$method == "wald_t", "wald", "lrt")
M <- merge(O, E[, c("config", "method", "alpha", "nu_delta", "size_ratio", "size_se", "draws")], by = c("config", "method", "alpha"), suffixes = c("", "_ref"))
M$z <- (M$size_ratio - M$size_ratio_ref) / sqrt((M$rel_se * M$size_ratio)^2 + M$size_se^2)
P <- M[M$scheme == "plain", c("config", "method", "alpha", "rel_se")]; names(P)[4] <- "rel_se_plain"
M <- merge(M, P); M$variance_gain <- (M$rel_se_plain / M$rel_se)^2
M <- M[order(M$nu_delta, M$method, -M$alpha, M$scheme), ]
write.csv(M, "calibration/results/importance_sampling_pilot.csv", row.names = FALSE)
options(width = 200); print(M[M$scheme == "mixture", c("config", "nu_delta", "method", "alpha", "size_ratio", "size_ratio_ref", "z", "ess", "wmax", "variance_gain")], digits = 3, row.names = FALSE)
