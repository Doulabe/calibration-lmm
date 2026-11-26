# ==============================================================================
# 60_skovgaard_rstar.R
# Modified signed likelihood root r* with Skovgaard's (1996) approximation, for the effect of
# the candidate in the mixed model, on the configurations and draws of
# 36_calibration_enlarged.R (same seeds, same order of the random numbers).
# The likelihood is that of the error contrasts s = H'y, which is free of the fixed effects F:
#   s ~ N(xi beta, tau_e I + tau_g L),  theta = (beta, tau_e, tau_g),  delta = tau_g / tau_e,
# with interest parameter beta and nuisance (tau_e, tau_g). With theta_hat the maximizer, theta_til
# the maximizer under beta = 0, r the signed root of the likelihood ratio and
#   S = Cov_hat{l'(theta_hat), l'(theta_til)},   q = Cov_hat{l'(theta_hat), l(theta_hat) - l(theta_til)},
# the covariances being taken under theta_hat,
#   u = [S^-1 q]_beta |j_hat|^{1/2} |i_hat|^{-1} |S| |j_til,nuisance|^{-1/2},   r* = r + log(u / r) / r,
# and r* is referred to N(0, 1). For this Gaussian model, with W = Sigma^-1 (diagonal in the
# spectral coordinates), d = xi beta_hat and c_e = 1, c_g = l_i the derivatives of Sigma,
#   S_bb = xi' W_til xi,  S_bk = xi' W_til C_k W_til d,  S_kb = 0,  S_jk = tr(C_j W_til C_k W_til) / 2,
#   q_b = xi' W_til d,    q_j = {tr(C_j W_til) - tr(C_j W_hat)} / 2.
# Sizes relative to alpha and powers at lambda = t_{nu-1, 1-alpha/2}, by the conditional
# calculation on a grid of the standardized effect u (as in 56). Rows:
#   lrc_chi2   likelihood ratio of the error contrasts, r^2 referred to chi2_1
#   rstar      r* referred to N(0, 1); r is used where u/r is not positive and for |r| < 1
# Mode check: direct simulation of complete null datasets at alpha = 1e-3 on two configurations,
# compared with the conditional calculation.
# Usage (from the repository root):
#   Rscript calibration/60_skovgaard_rstar.R [cores=1] [--smoke | --check]
# Output: results/skovgaard_rstar/<config>.rds, results/skovgaard_rstar.csv, results/skovgaard_rstar_check.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
a0 <- commandArgs(TRUE)
commandArgs <- function(...) character(0)
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # functions, GENES, load_geu, cfgs, GRID2, AL
NCORES <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 1L
SMOKE <- "--smoke" %in% a0; CHECK <- "--check" %in% a0             # set after 36, which defines its own
B_DRAW <- 1000L
UG <- c(seq(-40, -14, by = 2), seq(-12, 12, by = 0.2), seq(14, 40, by = 2))
R_MIN <- 1

## r and r* for the datasets in the columns of Sk (kernel coordinates), with the coordinate Sp of the
## candidate outside the kernel directions and the sum of squares Ssr of the remaining directions
rstar_stats <- function(sp, Sk, Sp, Ssr, pre) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n0 <- sp$nperp + sp$nrest; B <- ncol(Sk)
  Sk2 <- Sk * Sk
  A2 <- crossprod(Sk2, pre$W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, pre$W2) + xp * Sp
  lA <- log(A2); lR <- log(A2 - sweep(B2^2, 2, pre$S2g, "/"))
  dt <- pick_grid(-0.5 * sweep(nu * lA, 2, pre$ld2, "+"), GRID2)       # delta under beta = 0
  dh <- pick_grid(-0.5 * sweep(nu * lR, 2, pre$ld2, "+"), GRID2)       # delta in the model with the candidate
  Dh <- 1 + outer(ell, dh); Dt <- 1 + outer(ell, dt)                   # r x B
  Sh <- colSums(xi^2 / Dh) + xp^2; Bh <- colSums(xi * Sk / Dh) + xp * Sp; Ah <- colSums(Sk2 / Dh) + Sp^2 + Ssr
  At <- colSums(Sk2 / Dt) + Sp^2 + Ssr
  bh <- Bh / Sh; RSSh <- Ah - Bh^2 / Sh
  te <- RSSh / nu; tt <- At / nu                                       # tau_e at the two maximizers
  r2 <- colSums(log(Dt)) - colSums(log(Dh)) + nu * (log(At) - log(RSSh)); r2 <- pmax(r2, 0)
  r <- sign(bh) * sqrt(r2)
  wh <- sweep(1 / Dh, 2, te, "/"); wt <- sweep(1 / Dt, 2, tt, "/")     # weights W_hat, W_til on the kernel directions
  wh0 <- 1 / te; wt0 <- 1 / tt                                         # weights on the n0 directions with l = 0
  eh <- Sk - outer(xi, bh); ehp <- Sp - xp * bh                        # residuals at theta_hat
  cs <- function(M) colSums(M)
  ## S and q
  Sbb <- cs(xi^2 * wt) + xp^2 * wt0
  Sbe <- bh * (cs(xi^2 * wt^2) + xp^2 * wt0^2); Sbg <- bh * cs(xi^2 * ell * wt^2)
  See <- 0.5 * (cs(wt^2) + n0 * wt0^2); Seg <- 0.5 * cs(ell * wt^2); Sgg <- 0.5 * cs(ell^2 * wt^2)
  qb <- bh * Sbb
  qe <- 0.5 * (cs(wt) + n0 * wt0 - cs(wh) - n0 * wh0); qg <- 0.5 * (cs(ell * wt) - cs(ell * wh))
  detS <- See * Sgg - Seg^2
  num <- qb - (Sbe * (Sgg * qe - Seg * qg) + Sbg * (See * qg - Seg * qe)) / detS      # q_b - S_b,nuis S_nuis^-1 q_nuis
  ## expected and observed information at theta_hat
  ibb <- Sh / te
  iee <- 0.5 * (cs(wh^2) + n0 * wh0^2); ieg <- 0.5 * cs(ell * wh^2); igg <- 0.5 * cs(ell^2 * wh^2)
  deti <- ibb * (iee * igg - ieg^2)
  eh2 <- eh * eh; wh3 <- wh^3
  jbe <- cs(xi * eh * wh^2) + xp * ehp * wh0^2; jbg <- cs(xi * ell * eh * wh^2)
  jee <- -iee + cs(eh2 * wh3) + (ehp^2 + Ssr) * wh0^3; jeg <- -ieg + cs(ell * eh2 * wh3); jgg <- -igg + cs(ell^2 * eh2 * wh3)
  detj <- ibb * (jee * jgg - jeg^2) - jbe * (jbe * jgg - jeg * jbg) + jbg * (jbe * jeg - jee * jbg)
  ## observed information of the nuisance parameters at theta_til
  wt3 <- wt^3
  tee <- -See + cs(Sk2 * wt3) + (Sp^2 + Ssr) * wt0^3; teg <- -Seg + cs(ell * Sk2 * wt3); tgg <- -Sgg + cs(ell^2 * Sk2 * wt3)
  dett <- tee * tgg - teg^2
  u <- suppressWarnings(num * detS * sqrt(detj) / (deti * sqrt(dett)))
  ok <- is.finite(u) & r != 0 & u / r > 0
  ## r* is numerically unstable near r = 0, where log(u / r) / r divides two small numbers; the correction is
  ## applied for |r| >= R_MIN only, which leaves every rejection region considered here unchanged
  use <- ok & abs(r) >= R_MIN
  rs <- r; rs[use] <- r[use] + log(u[use] / r[use]) / r[use]
  list(r = r, rstar = rs, ok = ok | abs(r) < R_MIN, ratio = u / r, dh = dh, dt = dt, bh = bh)
}
precompute <- function(sp) { W2 <- 1 / (1 + outer(sp$ell, GRID2))
  list(W2 = W2, ld2 = colSums(log(1 + outer(sp$ell, GRID2))), S2g = colSums(sp$xi^2 * W2) + sp$xp^2) }
draw_r <- function(sp, delta0, B) {                          # draws of the error contrasts, as calib_config of 36
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  list(Rk = s - outer(xi, b), Rp = spp - xp * b, ssr = ssr, S0 = S0)
}
calib_rstar <- function(sp, delta0, tag, B = B_DRAW) {
  xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; pre <- precompute(sp); D <- draw_r(sp, delta0, B)
  ev <- function(u) { bv <- u / sqrt(D$S0); rstar_stats(sp, D$Rk + outer(xi, bv), D$Rp + xp * bv, D$ssr, pre) }
  lo <- rep(-30, B); hi <- rep(30, B)                       # zero of the effect estimate
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- ev(mid)$bh < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  nU <- length(UG); Fm <- list(lrc = matrix(0, B, nU), rstar = matrix(0, B, nU)); okm <- matrix(TRUE, B, nU)
  for (k in seq_len(nU)) { e <- ev(rep(UG[k], B)); Fm$lrc[, k] <- e$r^2; Fm$rstar[, k] <- e$rstar^2; okm[, k] <- e$ok }
  right <- outer(split, UG, "<"); kR0 <- rowSums(!right) + 1L
  cross <- function(key, c0) {                              # boundaries of {statistic >= c0} in u, zero statistic at the split
    Fk <- Fm[[key]]
    M <- (Fk >= c0) & right; any_r <- rowSums(M) > 0; k <- max.col(M, ties.method = "first"); k[!any_r] <- nU
    first <- k <= kR0; ka <- pmax(k - 1L, 1L)
    ua <- ifelse(first, split, UG[ka]); fa <- ifelse(first, 0, Fk[cbind(seq_len(B), ka)]); ub <- UG[k]; fb <- Fk[cbind(seq_len(B), k)]
    up <- pmin(pmax(ua + (c0 - fa) / pmax(fb - fa, 1e-300) * (ub - ua), ua), ub); up[!any_r] <- Inf
    M <- (Fk >= c0) & !right; any_l <- rowSums(M) > 0; k <- max.col(M, ties.method = "last"); k[!any_l] <- 1L
    last <- k >= kR0 - 1L; ka <- pmin(k + 1L, nU)
    ua <- ifelse(last, split, UG[ka]); fa <- ifelse(last, 0, Fk[cbind(seq_len(B), ka)]); ub <- UG[k]; fb <- Fk[cbind(seq_len(B), k)]
    dn <- pmax(pmin(ua + (c0 - fa) / pmax(fb - fa, 1e-300) * (ub - ua), ua), ub); dn[!any_l] <- -Inf
    list(up = up, dn = dn)
  }
  rej <- function(rt, lam) mean(pnorm(lam - rt$up) + pnorm(rt$dn - lam))
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  inner <- abs(UG) <= 12
  out <- do.call(rbind, lapply(seq_along(AL), function(k) do.call(rbind, lapply(c("lrc", "rstar"), function(key) {
    rt <- cross(key, cq[k])
    data.frame(tag, method = c(lrc = "lrc_chi2", rstar = "rstar")[[key]], alpha = AL[k], size_ratio = rej(rt, 0) / AL[k],
               power = rej(rt, tq_[k]), lambda = tq_[k]) }))))
  out$fallback <- mean(!okm[, inner])                       # share of evaluations where u/r is not positive
  ## monotonicity of r*^2 away from the zero of the effect, on the inner grid
  mono <- vapply(seq_len(B), function(b) { f <- Fm$rstar[b, ]; rr <- which(right[b, ] & inner); ll <- which(!right[b, ] & inner)
    all(diff(f[rr]) > -1e-6) && all(diff(f[ll]) < 1e-6) }, TRUE)
  out$monotone <- mean(mono)
  out
}
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else
  sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)
design <- function(z) {                                      # design, tag and seed of run_one of 36
  if (z$kind == "syn") {
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = cfg_id(z), design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    seed <- 4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L
  } else {
    g <- GENES[z$i, ]; W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j)); sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = cfg_id(z), design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    seed <- 9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973)
  }
  list(sp = sp, tag = tag, seed = seed)
}
suppressPackageStartupMessages(library(parallel))
ids <- vapply(cfgs, cfg_id, "")

if (CHECK) {                                                 # direct simulation against the conditional calculation
  WANT <- c("geu_ERAP2_96252589_h0.05", "syn_n200_b5_h0.05"); NDIR <- 400000L; CH <- 2000L; a <- 1e-3
  O <- do.call(rbind, lapply(WANT, function(id) { z <- cfgs[[match(id, ids)]]; d <- design(z); sp <- d$sp; pre <- precompute(sp)
    d0 <- z$h2 / (1 - z$h2); set.seed(d$seed); cond <- calib_rstar(sp, d0, d$tag)
    R <- Reduce(`+`, mclapply(seq_len(NDIR %/% CH), function(q) { set.seed(600000 + q); D <- draw_r(sp, d0, CH)
      bv <- rnorm(CH) / sqrt(D$S0); e <- rstar_stats(sp, D$Rk + outer(sp$xi, bv), D$Rp + sp$xp * bv, D$ssr, pre)
      c(sum(e$r^2 > qchisq(1 - a, 1)), sum(e$rstar^2 > qchisq(1 - a, 1)), sum(!e$ok)) }, mc.cores = NCORES))
    p <- R[1:2] / NDIR
    data.frame(config = id, method = c("lrc_chi2", "rstar"), alpha = a, direct = p / a, se = sqrt(p * (1 - p) / NDIR) / a,
               conditional = cond$size_ratio[cond$alpha == a], fallback_direct = R[3] / NDIR) }))
  write.csv(O, "calibration/results/skovgaard_rstar_check.csv", row.names = FALSE); print(O, digits = 4); quit(save = "no")
}
OUT <- "calibration/results/skovgaard_rstar"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_one <- function(z) { ck <- file.path(OUT, paste0(cfg_id(z), ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
  d <- design(z); set.seed(d$seed); out <- calib_rstar(d$sp, z$h2 / (1 - z$h2), d$tag); if (SMOKE) return(out); saveRDS(out, ck); ck }
if (SMOKE) {
  for (id in c("geu_ERAP2_96252589_h0.05", "syn_n200_b20_h0.05_unlinked", "geu_ERAP2_95943006_h0.05")) { t0 <- Sys.time(); out <- run_one(cfgs[[match(id, ids)]])
    cat(sprintf("\n%s  nu_delta = %.1f  (%.0f s)\n", id, out$nu_delta[1], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(out[, c("method", "alpha", "size_ratio", "power", "fallback", "monotone")], row.names = FALSE, digits = 4) }
} else {
  cat("configurations:", length(cfgs), "\n")
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) paste("ERROR", conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT, "rds$", full.names = TRUE)
  write.csv(do.call(rbind, lapply(fs, readRDS)), "calibration/results/skovgaard_rstar.csv", row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
