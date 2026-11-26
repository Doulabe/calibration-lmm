# ==============================================================================
# 43_nu_plugin_boundary.R
# Two checks on the configurations of the calibration study (36_calibration_enlarged.R),
# with B null datasets drawn in the reduced model at the true ratio delta0:
#  (1) the index nu_j, defined at delta0, evaluated at the estimates used in practice, delta_hat (REML
#      with the candidate) and delta_hat_0 (REML of the null model): share of the draws in which it
#      falls in the same regime as nu_delta (below 20, 20 to 100, 100 and above), median ratio to
#      nu_delta, and the minimum of nu_j over a plausible range of ratios, h2 in [0.05, 0.5];
#  (2) the boundary of the domain: share of the null draws with delta_hat = 0 and their share of the
#      exact tail of the Wald test at 5e-8 (conditional probabilities as in 36).
# Usage (from the repository root): nice -n 10 Rscript calibration/43_nu_plugin_boundary.R [cores=3]
# Output: results/nu_plugin_boundary.csv (one row per configuration)
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # GENES, cfgs, load_geu, gen_blocks, GRID2, pick_grid
a0 <- commandArgs(TRUE)
NC <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 3L
B43 <- 1000L
H2R <- 10^seq(log10(0.05 / 0.95), 0, length.out = 25)            # delta for h2 from 0.05 to 0.5
regime <- function(v) cut(v, c(0, 20, 100, Inf), right = FALSE, labels = c("<20", "20-100", ">=100"))
nu_at <- function(sp, d) vapply(d, function(x) indices(sp, max(x, 1e-6))[["nu_delta"]], 0)
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else
  sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)

one_cfg <- function(k) {
  z <- cfgs[[k]]; id <- cfg_id(z)
  if (z$kind == "syn") {                                       # data as in 36
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    sp <- if (z$unlinked) { set.seed(77 + z$n); spectral_setup(Z, x = rnorm(z$n)) } else spectral_setup(Z, j = 25L)
  } else { w <- load_geu(z$i); sp <- spectral_setup(w$Z, j = match(z$pos, w$pos)) }
  set.seed(43000 + k); B <- B43
  d0 <- z$h2 / (1 - z$h2); ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  dd <- 1 + d0 * ell; S0 <- sum(xi^2 / dd) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(dd); spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / dd) * s) + xp * spp) / S0; Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  A_r <- crossprod(Rk^2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  dh1 <- pick_grid(-0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"), GRID2)   # REML with the candidate
  u0 <- rnorm(B); sk <- Rk + outer(xi, u0 / sqrt(S0)); s2p <- Rp + xp * u0 / sqrt(S0)  # one null dataset per draw
  A2 <- crossprod(sk^2, W2) + (s2p^2 + ssr)
  dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)                    # REML of the null model
  iD <- 1 / (1 + outer(ell, dh1)); A <- colSums(Rk^2 * iD) + Rp^2 + ssr; Bx <- colSums(xi * Rk * iD) + xp * Rp
  S2 <- colSums(xi^2 * iD) + xp^2; sig <- sqrt((A - Bx^2 / S2) / (nu - 1))
  wa <- sqrt(S2 / S0) / sig; wg <- (Bx / S2) * sqrt(S2) / sig; cc <- qt(1 - 2.5e-8, nu - 1)
  pr <- pnorm(-(cc - wg) / wa) + pnorm((-cc - wg) / wa)                                # Wald tail given r, 5e-8
  at0 <- dh1 < 1e-8
  nd <- indices(sp, d0)[["nu_delta"]]; n1 <- nu_at(sp, dh1); n0 <- nu_at(sp, dh0); nmin <- min(nu_at(sp, H2R))
  data.frame(config = id, kind = z$kind, h2 = z$h2, nu_delta = nd, regime = as.character(regime(nd)),
             agree_dh1 = mean(regime(n1) == regime(nd)), agree_dh0 = mean(regime(n0) == regime(nd)),
             ratio_dh1 = median(n1 / nd), ratio_dh0 = median(n0 / nd),
             nu_min = nmin, min_same_regime = regime(nmin) == regime(nd), min_le_nu_delta = nmin <= nd * (1 + 1e-6),
             p_dhat0 = mean(at0), tail_share_dhat0 = if (sum(pr) > 0) sum(pr[at0]) / sum(pr) else NA_real_,
             size_ratio_5e8 = mean(pr) / 5e-8)
}
suppressPackageStartupMessages(library(parallel))
O <- do.call(rbind, mclapply(seq_along(cfgs), function(k) tryCatch(one_cfg(k), error = function(e) { message("FAIL ", k, " ", conditionMessage(e)); NULL }),
                             mc.cores = NC, mc.preschedule = FALSE))
write.csv(O, "calibration/results/nu_plugin_boundary.csv", row.names = FALSE)
O$regime <- factor(O$regime, levels = c("<20", "20-100", ">=100"))
options(width = 200)
cat("configurations:", nrow(O), "\n")
cat("== nu_j at the estimates: regime agreement with nu_delta and median ratio (median and minimum over configurations) ==\n")
print(aggregate(cbind(agree_dh1, agree_dh0) ~ regime, O, function(v) round(c(med = median(v), min = min(v)), 3)))
print(aggregate(cbind(ratio_dh1, ratio_dh0) ~ regime, O, function(v) round(median(v), 3)))
cat("== minimum of nu_j over h2 in [0.05, 0.5] ==\n")
print(aggregate(cbind(min_same_regime, min_le_nu_delta) ~ regime, O, mean))
cat("== boundary delta_hat = 0: share of null draws and share of the Wald tail at 5e-8 (median and range) ==\n")
print(aggregate(cbind(p_dhat0, tail_share_dhat0) ~ regime, O, function(v) round(c(med = median(v), min = min(v), max = max(v)), 3)))
