# ==============================================================================
# 55_direct_simulation_check.R
# Direct simulation of the Wald test at alpha = 1e-4 on two configurations of the calibration
# study (36_calibration_enlarged.R). Each null dataset is drawn in full, delta_hat maximizes the
# restricted profile of the model with the candidate on that dataset, and the rejection indicator
# of |t(delta_hat)| is recorded, without the decomposition (b, r). The rejection rate is compared
# with the Rao-Blackwellized estimate on independent draws of r.
# Usage (from the repository root): nice -n 10 Rscript calibration/55_direct_simulation_check.R [cores=4]
# Output: results/direct_simulation_check.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # GENES, cfgs, load_geu, gen_blocks, GRID2, pick_grid
a0 <- commandArgs(TRUE); NC <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 4L
suppressPackageStartupMessages(library(parallel))
ALPHA <- 1e-4; NDIR <- 2e6L; NRB <- 2e5L; CH <- 5000L
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else
  sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)
ids <- vapply(cfgs, cfg_id, ""); WANT <- c(ids[1], "syn_n200_b5_h0.05")
one_cfg <- function(id) {
  k <- match(id, ids); z <- cfgs[[k]]
  if (z$kind == "syn") { set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
  } else { w <- load_geu(z$i); sp <- spectral_setup(w$Z, j = match(z$pos, w$pos)) }
  d0 <- z$h2 / (1 - z$h2); ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  dd <- 1 + d0 * ell; S0 <- sum(xi^2 / dd) + xp^2; cc <- qt(1 - ALPHA / 2, nu - 1)
  wald_at <- function(Sk, Sp, ssr) {        # delta_hat and Wald statistic computed from the coordinates of one dataset
    A <- crossprod(Sk^2, W2) + (Sp^2 + ssr); Bm <- crossprod(xi * Sk, W2) + xp * Sp
    dh <- pick_grid(-0.5 * sweep((nu - 1) * log(A - sweep(Bm^2, 2, S2g, "/")), 2, ld2 + log(S2g), "+"), GRID2)
    iD <- 1 / (1 + outer(ell, dh)); A1 <- colSums(Sk^2 * iD) + Sp^2 + ssr; B1 <- colSums(xi * Sk * iD) + xp * Sp
    S2 <- colSums(xi^2 * iD) + xp^2; sig <- sqrt((A1 - B1^2 / S2) / (nu - 1))
    list(t = (B1 / S2) * sqrt(S2) / sig, S2 = S2, sig = sig, B1 = B1)
  }
  draw <- function(B) list(s = matrix(rnorm(r * B), r, B) * sqrt(dd), spp = if (sp$nperp) rnorm(B) else numeric(B),
                           ssr = if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B))
  nch <- NDIR %/% CH
  rej <- unlist(mclapply(seq_len(nch), function(q) { set.seed(55000 + 1000 * k + q); D <- draw(CH)
    sum(abs(wald_at(D$s, D$spp, D$ssr)$t) > cc) }, mc.cores = NC))
  pr <- unlist(mclapply(seq_len(NRB %/% CH), function(q) { set.seed(55500000 + 1000 * k + q); D <- draw(CH)   # Rao-Blackwellized, independent draws
    b <- (colSums((xi / dd) * D$s) + xp * D$spp) / S0; Rk <- D$s - outer(xi, b); Rp <- D$spp - xp * b
    w <- wald_at(Rk, Rp, D$ssr); wa <- sqrt(w$S2 / S0) / w$sig; pnorm((w$t - cc) / wa) + pnorm((-cc - w$t) / wa) }, mc.cores = NC))
  pd <- sum(rej) / NDIR
  data.frame(config = id, nu_delta = indices(sp, d0)[["nu_delta"]], alpha = ALPHA, n_direct = NDIR, rejections = sum(rej),
             size_ratio_direct = pd / ALPHA, se_direct = sqrt(pd * (1 - pd) / NDIR) / ALPHA,
             n_rb = NRB, size_ratio_rb = mean(pr) / ALPHA, se_rb = sd(pr) / sqrt(NRB) / ALPHA)
}
O <- do.call(rbind, lapply(WANT, one_cfg))
O$z <- (O$size_ratio_direct - O$size_ratio_rb) / sqrt(O$se_direct^2 + O$se_rb^2)
write.csv(O, "calibration/results/direct_simulation_check.csv", row.names = FALSE)
options(width = 200); print(O, digits = 4)
