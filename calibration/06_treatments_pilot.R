# ==============================================================================
# 06_treatments_pilot.R
#
# Which treatment of the variance ratio delta gives calibrated evidence for a
# candidate in LD with the regional background, and at what power cost?
# Same spectral reduction and configurations as 02_score_null_exact.R.
#
# Every statistic is a deterministic function of (u, r), where r are the error
# contrasts (law free of beta) and u = b sqrt(S(delta0)) ~ N(lambda, 1) is
# independent of r. For each draw of r, the rejection region in u is found
# (closed form or bisection), so size (lambda = 0) and power (lambda > 0) are
# exact Rao-Blackwellised averages, P(reject | r) = Phi(lambda - u+) + Phi(u- - lambda).
#
# Treatments (nominal level alpha):
#   oracle_t   t(delta0), delta known                         (benchmark, exact)
#   wald_h1    t(delta_hat), REML with the candidate (GEMMA)  cut-off t_{nu-1}
#   score_h0   Q(delta_hat0), REML of the null (EMMAX, MLMM)  cut-off t_{nu-1}
#   lrt_ml     profile-ML LRT, delta re-estimated under each model,
#              cut-off = exact LRT quantile when delta is known
#   signflip   spectral sign-flip test (exact given |s|, any delta), statistic
#              = score numerator at delta_hat0, Lugannani-Rice/r* tail
#   bb_wald    Berger-Boos: sup of the exact Wald p-value over a REML-profile
#              confidence set for delta (gamma = alpha/10), tabulated exact law
# Bayes factors (slab tau2 = 0.04), thresholds matched to the
# alpha-quantile of the oracle BF (delta known):
#   bf_h1      BF at delta_hat (REML with the candidate)
#   bf_h0      BF at delta_hat0 (REML arm)
#   mbf        BF marginalised over delta under both hypotheses
#              (half-Cauchy prior on sqrt(delta), log-delta in [1e-4,1e3])
# eBIC decision: LRT_ml > log n + 2 log m, compared with the same
# rule applied to the oracle LRT (delta known).
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/06_treatments_pilot.R [cores=2] [--smoke]
# Output: calibration/results/treatments_pilot.csv (long format)
# ==============================================================================

args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
SMOKE  <- "--smoke" %in% args
B_DRAW <- if (SMOKE) 200L else 1000L
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))

GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))       # REML / ML maximisation grid
QG    <- 10^seq(-4, 3, length.out = 61)             # MBF quadrature range
TAU2  <- 0.04
AL    <- c(1e-3, 1e-5, 5e-8)
NBIS  <- 28L
CGRID <- seq(0, 80, by = 0.02)

pick_grid <- function(O, grid) {
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
rowLSE <- function(M) { mx <- apply(M, 1, max); mx + log(rowSums(exp(M - mx))) }

## log two-sided p-value of the sign-flip test: T* = sum a_i eps_i, eps Rademacher.
## Lugannani-Rice through the r* form of Barndorff-Nielsen (always a valid probability).
log_signflip_p <- function(A, x) {                  # A: k x B nonneg weights, x: |T_obs|
  K2_0 <- colSums(A^2); tot <- colSums(A)
  t <- pmax(x / K2_0, 1e-12)
  for (it in 1:60) {
    th <- tanh(sweep(A, 2, t, "*"))
    K1 <- colSums(A * th); K2 <- colSums(A^2 * (1 - th^2))
    tn <- t - (K1 - x) / pmax(K2, 1e-300)
    tn <- ifelse(tn <= 0 | !is.finite(tn), t * 0.5, tn)
    tn <- ifelse(K1 < x & tn < t, t * 2, tn)       # keep moving right when below target
    conv <- max(abs(tn - t) / pmax(t, 1e-12)); t <- tn
    if (conv < 1e-10) break
  }
  z  <- sweep(A, 2, t, "*"); az <- abs(z)
  Kt <- colSums(az + log1p(exp(-2 * az)) - log(2))
  th <- tanh(z); K2 <- colSums(A^2 * (1 - th^2))
  w  <- sqrt(pmax(2 * (t * x - Kt), 0)); uu <- t * sqrt(K2)
  rs <- w + log(uu / w) / w
  lp <- ifelse(w < 1e-3, log(0.5), pnorm(-rs, log.p = TRUE))
  lp[x >= tot * (1 - 1e-12)] <- -colSums(A > 0)[x >= tot * (1 - 1e-12)] * log(2)
  pmin(log(2) + lp, 0)
}

## exact Wald law at delta (same algebra as 01/03): per-draw (m, v)
draws_mv <- function(sp, delta, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu0 <- sp$nu0; r <- sp$r
  W  <- 1 / (1 + outer(ell, GRID2)); ld <- colSums(log(1 + outer(ell, GRID2)))
  Sg <- colSums(xi^2 * W) + xp^2
  d0 <- 1 + delta * ell
  s  <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  sp_ <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  s2 <- s * s
  Ag <- crossprod(s2, W) + (sp_^2 + ssr); Bg <- crossprod(xi * s, W) + xp * sp_
  O  <- -0.5 * (sweep((nu0 - 1) * log(Ag - sweep(Bg^2, 2, Sg, "/")), 2, ld + log(Sg), "+"))
  dh <- pick_grid(O, GRID2)
  at <- function(dl) {
    iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(s2 * iD) + sp_^2 + ssr; Bx <- colSums(xi * s * iD) + xp * sp_
    S2 <- colSums(xi^2 * iD) + xp^2
    list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
  }
  Q <- at(dh); Q0 <- at(rep(delta, B)); sig2 <- Q$RSS / (nu0 - 1)
  list(m = (Q$Bx / Q$S2 - Q0$Bx / Q0$S2) * sqrt(Q$S2 / sig2), v = Q$S2 / (Q0$S2 * sig2))
}
crit_exact <- function(sp, delta, alphas, B = 1500L) {
  D <- draws_mv(sp, delta, B); sd_ <- sqrt(D$v)
  G <- vapply(CGRID, function(c0) mean(pnorm((D$m - c0) / sd_) + pnorm((-c0 - D$m) / sd_)), 0)
  lg <- log(pmax(G, 1e-300))
  vapply(log(alphas), function(la) approx(rev(lg), rev(CGRID), xout = la, ties = mean)$y, 0)
}

treat_config <- function(sp, delta0, m_mark, tag, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  WQ <- 1 / (1 + outer(ell, QG)); ldQ <- colSums(log(1 + outer(ell, QG)))
  SQg <- colSums(xi^2 * WQ) + xp^2; vQ <- TAU2 * SQg
  tq <- log(QG); dq <- tq[2] - tq[1]
  lwq <- log(c(dq / 2, rep(dq, length(tq) - 2), dq / 2)) +
         0.5 * tq - log(pi) - log1p(QG)                  # trapezoid x half-Cauchy(sqrt delta) in t
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  ## ---- error contrasts ---------------------------------------------------------
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1  <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+")
  dh1 <- pick_grid(l1, GRID2)
  ml1 <- apply(-sweep(n * log(RSS_r), 2, ld2, "+"), 1, max)
  at_r <- function(dl) {
    iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Rk2 * iD) + Rp^2 + ssr; Bx <- colSums(xi * Rk * iD) + xp * Rp
    S2 <- colSums(xi^2 * iD) + xp^2
    list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
  }
  Wd <- at_r(dh1); sig_h <- sqrt(Wd$RSS / (nu - 1))
  wa <- sqrt(Wd$S2 / S0) / sig_h                 # t(delta_hat) = wa * u + wg
  wg <- (Wd$Bx / Wd$S2) * sqrt(Wd$S2) / sig_h
  O0 <- at_r(rep(delta0, B)); RSS0 <- O0$RSS; sig0 <- sqrt(RSS0 / (nu - 1))
  ## ---- statistics as functions of u (vectorised over draws) --------------------
  ev <- function(u, what) {
    bv <- u / sqrt(S0)
    sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    out <- list()
    if (any(what %in% c("Q", "bf0", "sf", "lrt"))) {
      A2 <- crossprod(sk2, W2) + (s2p^2 + ssr)
      if ("lrt" %in% what) out$lrt <- ml1 - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max)
      if (any(what %in% c("Q", "bf0", "sf"))) {
        dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
        iD <- 1 / (1 + outer(ell, dh0))
        Ah <- colSums(sk2 * iD) + s2p^2 + ssr
        Bh <- colSums(xi * sk * iD) + xp * s2p
        Sh <- colSums(xi^2 * iD) + xp^2
        out$Q <- Bh / sqrt(Sh * Ah / nu)
        if ("bf0" %in% what) {
          v <- TAU2 * Sh
          out$bf0 <- -0.5 * log1p(v) - (nu / 2) * log1p(-v * out$Q^2 / (nu * (1 + v)))
        }
        if ("sf" %in% what) {
          Aw <- rbind(abs(xi * iD * sk), abs(xp * s2p))
          out$sf <- -log_signflip_p(Aw, abs(Bh))       # -log p (U-shaped in u)
        }
      }
    }
    if ("mbf" %in% what) {
      AQ <- crossprod(sk2, WQ) + (s2p^2 + ssr); BQ <- crossprod(xi * sk, WQ) + xp * s2p
      l0 <- -0.5 * sweep(nu * log(AQ), 2, ldQ, "+")
      lbf <- sweep(-(nu / 2) * log1p(-sweep(BQ^2 / AQ, 2, TAU2 / (1 + vQ), "*")), 2,
                   0.5 * log1p(vQ), "-")
      lw <- sweep(l0, 2, lwq, "+")
      out$mbf <- rowLSE(lw + lbf) - rowLSE(lw)
    }
    out
  }
  ## split point: Q(delta_hat0(u)) = 0 (Q increasing in u)
  lo <- rep(-30, B); hi <- rep(30, B)
  for (it in 1:40) { mid <- (lo + hi) / 2; q <- ev(mid, "Q")$Q; up <- q < 0
    lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  rootsU <- function(key, c0) {                   # U-shaped statistic, level c0
    f <- function(u) ev(u, key)[[key]]
    always <- f(split) >= c0
    lo <- split; hi <- rep(120, B); never_r <- f(hi) < c0
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- f(mid) < c0
      lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    ur <- (lo + hi) / 2; ur[never_r] <- Inf
    lo <- rep(-120, B); hi <- split; never_l <- f(lo) < c0
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; dn <- f(mid) < c0
      hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
    ul <- (lo + hi) / 2; ul[never_l] <- -Inf
    ur[always] <- split[always]; ul[always] <- split[always]
    list(up = ur, dn = ul)
  }
  rootsQ <- function(c0) {                          # monotone score, crossings of +-c0
    one <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
      for (it in seq_len(NBIS + 6)) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < level
        lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      (lo + hi) / 2 }
    list(up = one(c0), dn = one(-c0))
  }
  lin <- function(a, g, c0) list(up = (c0 - g) / a, dn = (-c0 - g) / a)  # |a u + g| > c0
  rej <- function(rt, lam) mean(pnorm(lam - rt$up) + pnorm(rt$dn - lam))

  ## ---- thresholds ------------------------------------------------------------------
  tq_ <- qt(1 - AL / 2, df = nu - 1)
  qa  <- tq_ * sqrt(nu / (nu - 1 + tq_^2))                     # oracle |Q| quantile
  v0  <- TAU2 * S0
  cBF <- -0.5 * log1p(v0) - (nu / 2) * log1p(-v0 * qa^2 / (nu * (1 + v0)))
  cLR <- n * log1p(tq_^2 / (nu - 1))
  pen <- log(n) + 2 * log(m_mark)
  ## Berger-Boos: exact Wald critical values over delta, at level 0.9 alpha
  TG  <- c(0, 10^seq(-3, 3.5, length.out = 27))
  CT  <- t(vapply(TG, function(dl) crit_exact(sp, dl, 0.9 * AL), numeric(length(AL))))
  cG  <- apply(CT, 2, function(y) {
    out <- numeric(length(GRID2)); pos <- GRID2 > 0
    out[!pos] <- y[1]
    out[pos] <- approx(log(TG[-1]), y[-1], xout = log(pmin(pmax(GRID2[pos], TG[2]), max(TG))),
                       rule = 2)$y
    out })
  lmax <- apply(l1, 1, max)
  res <- list()
  add <- function(method, k, rt) {
    lam <- tq_[k]
    res[[length(res) + 1]] <<- data.frame(tag, method = method, alpha = AL[k],
      size_ratio = rej(rt, 0) / AL[k], power = rej(rt, lam), lambda = lam)
  }
  for (k in seq_along(AL)) {
    add("oracle_t", k, list(up = tq_[k] * sig0, dn = -tq_[k] * sig0))
    add("wald_h1", k, lin(wa, wg, tq_[k]))
    add("score_h0", k, rootsQ(tq_[k]))
    add("lrt_ml", k, rootsU("lrt", cLR[k]))
    add("signflip", k, rootsU("sf", -log(AL[k])))
    inCI <- 2 * (lmax - l1) <= qchisq(1 - 0.1 * AL[k], 1)
    cbb <- apply(inCI, 1, function(z) max(cG[z, k]))
    add("bb_wald", k, lin(wa, wg, cbb))
    ## Bayes factors at the oracle-matched threshold
    Q2c <- (nu / (nu - qa[k]^2))                                    # helper for oracle BF
    add("oracle_bf", k, list(up = qa[k] * sqrt(RSS0 / (nu - qa[k]^2)),
                             dn = -qa[k] * sqrt(RSS0 / (nu - qa[k]^2))))
    vh <- TAU2 * Wd$S2
    Qc2 <- nu * (1 + vh) / vh * (1 - exp(-(2 * cBF[k] + log1p(vh)) / nu))
    tc  <- ifelse(Qc2 < nu, sqrt((nu - 1) * Qc2 / (nu - Qc2)), Inf)
    add("bf_h1", k, lin(wa, wg, tc))
    add("bf_h0", k, rootsU("bf0", cBF[k]))
    add("mbf", k, rootsU("mbf", cBF[k]))
  }
  ## eBIC decision (stopping rule of the selection paths), with the oracle LRT at the same penalty
  te <- sqrt((nu - 1) * expm1(pen / n))
  ro <- list(up = te * sig0, dn = -te * sig0); rl <- rootsU("lrt", pen)
  for (lam in c(0, tq_[1], tq_[3])) {
    res[[length(res) + 1]] <- data.frame(tag, method = c("ebic_oracle", "ebic_lrt_ml"),
      alpha = NA, size_ratio = NA, power = c(rej(ro, lam), rej(rl, lam)), lambda = lam)
  }
  do.call(rbind, res)
}

## ---- configurations ----------------------------------------------------------------
cfgs <- list()
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2)
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
for (i in seq_len(nrow(GENES))) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "geu", i = i, h2 = h2)
if (SMOKE) cfgs <- cfgs[c(which(vapply(cfgs, function(z) z$kind == "geu" && z$i == 1 &&
                                          z$h2 == 0.05, TRUE)))]

gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t",
                fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, j = match(g$lead_pos, pos))
}

OUT <- "calibration/results/treatments_pilot"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", n = z$n, blocks = z$nb,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2))
    out <- treat_config(sp, z$h2 / (1 - z$h2), 50L * z$nb, tag)
  } else {
    id <- sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", n = nrow(L$Z), blocks = NA,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(7000 + z$i + round(100 * z$h2))
    out <- treat_config(sp, z$h2 / (1 - z$h2), ncol(L$Z), tag)
  }
  if (!SMOKE) saveRDS(out, ck)
  cat(sprintf("[%7.1f s] %s done\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), id))
  if (SMOKE) out else ck
}
if (SMOKE) {
  out <- run_one(cfgs[[1]])
  op <- options(width = 200); print(out[, c("method", "alpha", "size_ratio", "power", "lambda")],
                                    digits = 3); options(op)
} else {
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e)
    paste("ERROR", conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT, "rds$", full.names = TRUE)
  all <- do.call(rbind, lapply(fs, readRDS))
  write.csv(all, "calibration/results/treatments_pilot.csv", row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
cat("done in", format(Sys.time() - t0), "\n")
