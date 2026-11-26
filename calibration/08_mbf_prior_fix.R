# ==============================================================================
# 08_mbf_prior_fix.R
#
# Correction of the non-monotonicity of the marginal Bayes factor (MBF).
# With a fixed slab variance tau2, the H1 marginal likelihood of a very
# large effect decays like the Gaussian slab (sigma2 has to inflate), while the
# null model absorbs the effect into a large delta at a polynomial cost, so the
# MBF peaks and falls back towards 1 (information paradox of fixed g-priors,
# Liang et al. 2008). A heavy-tailed prior on the effect scale removes it.
#
# Versions compared (delta integrated under a half-Cauchy prior on sqrt(delta),
# log-delta in [1e-4, 1e3]):
#   mbf_fixed  tau2 = 0.04 fixed (results taken from 06b)
#   mbf_hc     tau ~ half-Cauchy(0, 0.2), i.e. the fixed value as prior scale,
#              v = tau2 * S(delta)  (same parameterisation as mbf_fixed)
#   mbf_hg     hyper-g prior, pi(g) = (a-2)/2 (1+g)^(-a/2), a = 3, with the
#              Zellner scaling v = g in the whitened metric (Liang et al. 2008)
# Each version is thresholded at the (1 - alpha) quantile of its own oracle
# (delta known), computed at |Q(delta0)| = q_alpha. Rejection regions are found
# as in 06b (side maximum by golden section, then two bisections), on the SAME
# error-contrast draws as 06/06b (same seeds and draw order).
# A diagnostic prints log MBF against u on four ERAP2 draws for the three versions.
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/08_mbf_prior_fix.R [cores=4] [--diag_only]
# Output: results/mbf_prior_fix/<config>.rds, results/mbf_prior_fix.csv,
#         results/mbf_prior_diag.csv
# ==============================================================================

args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 4L
DIAG_ONLY <- "--diag_only" %in% args
src <- readLines("calibration/06_treatments_pilot.R")
cut <- grep("^## ---- configurations", src)
eval(parse(text = src[1:(cut - 1)]))           # helpers and grids of 06
B_DRAW <- 1000L
GMAX <- 25; GFALL <- 30; NGOLD <- 32L; NB2 <- 26L
TAU0 <- 0.2; A_HG <- 3
T2  <- 10^seq(-6, 4, length.out = 41)          # tau2 nodes (half-Cauchy on tau)
GG  <- 10^seq(-4, 6, length.out = 41)          # g nodes (hyper-g)
lw_t2 <- { t <- log(T2); dt <- t[2] - t[1]
  w <- log(c(dt / 2, rep(dt, length(t) - 2), dt / 2)) - log(pi * TAU0) - log1p(T2 / TAU0^2) + 0.5 * t
  w - (max(w) + log(sum(exp(w - max(w))))) }  # normalised on the grid
lw_g <- { t <- log(GG); dt <- t[2] - t[1]
  w <- log(c(dt / 2, rep(dt, length(t) - 2), dt / 2)) + log((A_HG - 2) / 2) - (A_HG / 2) * log1p(GG) + t
  w - (max(w) + log(sum(exp(w - max(w))))) }

## log of the H1 marginal likelihood integrated over delta AND the effect-scale prior,
## vectorised over all (delta node, scale node) pairs in a single log-sum-exp.
## R = B^2/A (draws x delta nodes), l0w = l0 + log prior weight of delta (draws x delta nodes),
## S = S(delta) at the nodes, nu = residual df. Returns a vector over draws.
node_pairs <- function(S, version) {
  nq <- length(S); K <- if (version == "hc") length(T2) else length(GG)
  iq <- rep(seq_len(nq), each = K); ik <- rep(seq_len(K), times = nq)
  if (version == "hc") {
    v <- T2[ik] * S[iq]; list(iq = iq, lv = log1p(v), cc = T2[ik] / (1 + v), lw = lw_t2[ik])
  } else {
    list(iq = iq, lv = log1p(GG[ik]), cc = (GG[ik] / (1 + GG[ik])) / S[iq], lw = lw_g[ik])
  }
}
logp1_mix <- function(R, l0w, np, nu, chunk = 200L) {   # chunked over draws (memory)
  B <- nrow(R); J <- length(np$iq); out <- numeric(B)
  for (st in seq(1L, B, by = chunk)) {
    ii <- st:min(B, st + chunk - 1L); b <- length(ii)
    M <- -0.5 * matrix(rep(np$lv, each = b), b, J) -
         (nu / 2) * log1p(-R[ii, np$iq, drop = FALSE] * rep(np$cc, each = b)) +
         l0w[ii, np$iq, drop = FALSE] + rep(np$lw, each = b)
    out[ii] <- rowLSE(M)
  }
  out
}
## oracle log BF of a version at delta0 (single delta, scale prior integrated)
lbf_oracle <- function(Rval, S0, nu, version) {
  np <- node_pairs(S0, version)
  L <- -0.5 * np$lv - (nu / 2) * log1p(-Rval * np$cc) + np$lw
  max(L) + log(sum(exp(L - max(L))))
}

mbf_config <- function(sp, delta0, tag, B = B_DRAW, versions = c("hc", "hg")) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  WQ <- 1 / (1 + outer(ell, QG)); ldQ <- colSums(log(1 + outer(ell, QG)))
  SQg <- colSums(xi^2 * WQ) + xp^2
  tq <- log(QG); dq <- tq[2] - tq[1]
  lwq <- log(c(dq / 2, rep(dq, length(tq) - 2), dq / 2)) + 0.5 * tq - log(pi) - log1p(QG)
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  NP <- list(hc = node_pairs(SQg, "hc"), hg = node_pairs(SQg, "hg"))
  mbf_u <- function(u, version) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv
    AQ <- crossprod(sk * sk, WQ) + (s2p^2 + ssr); BQ <- crossprod(xi * sk, WQ) + xp * s2p
    l0w <- sweep(-0.5 * sweep(nu * log(AQ), 2, ldQ, "+"), 2, lwq, "+")
    logp1_mix(BQ^2 / AQ, l0w, NP[[version]], nu) - rowLSE(l0w)
  }
  numQ <- function(u) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv
    A2 <- crossprod(sk * sk, W2) + (s2p^2 + ssr)
    dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
    colSums(xi * sk / (1 + outer(ell, dh0))) + xp * s2p
  }
  lo <- rep(-30, B); hi <- rep(30, B)
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- numQ(mid) < 0
    lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  gr <- (sqrt(5) - 1) / 2
  side <- function(version, dir, cvec) {
    f <- function(g) mbf_u(split + dir * g, version)
    a <- rep(0, B); z <- rep(GMAX, B)
    x1 <- z - gr * (z - a); x2 <- a + gr * (z - a); f1 <- f(x1); f2 <- f(x2)
    for (it in seq_len(NGOLD)) {
      mv <- f1 < f2
      an <- ifelse(mv, x1, a); zn <- ifelse(mv, z, x2)
      x1n <- ifelse(mv, x2, zn - gr * (zn - an)); x2n <- ifelse(mv, an + gr * (zn - an), x1)
      fe <- f(ifelse(mv, x2n, x1n))
      f1n <- ifelse(mv, f2, fe); f2n <- ifelse(mv, fe, f1)
      a <- an; z <- zn; x1 <- x1n; x2 <- x2n; f1 <- f1n; f2 <- f2n
    }
    gpk <- (a + z) / 2; fpk <- f(gpk)
    f0 <- f(rep(0, B)); fend <- f(rep(GMAX, B))
    use0 <- f0 > fpk; gpk[use0] <- 0; fpk[use0] <- f0[use0]
    useE <- fend > fpk; gpk[useE] <- GMAX; fpk[useE] <- fend[useE]
    ffall <- f(rep(GFALL, B))
    interior <- gpk > 0.05 & gpk < GMAX - 0.05
    lapply(cvec, function(c0) {
      has <- fpk >= c0
      lo <- rep(0, B); hi <- gpk; start_in <- f0 >= c0
      for (it in seq_len(NB2)) { mid <- (lo + hi) / 2; up <- f(mid) < c0
        lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      g1 <- ifelse(start_in, 0, (lo + hi) / 2)
      lo <- gpk; hi <- rep(GFALL, B); stays <- ffall >= c0
      for (it in seq_len(NB2)) { mid <- (lo + hi) / 2; above <- f(mid) >= c0
        lo[above] <- mid[above]; hi[!above] <- mid[!above] }
      g2 <- ifelse(stays, Inf, (lo + hi) / 2)
      g1[!has] <- NA; g2[!has] <- NA
      list(g1 = g1, g2 = g2, interior = interior)
    })
  }
  prob <- function(R, L, lam) {
    pr <- ifelse(is.na(R$g1), 0, pnorm(split + R$g2 - lam) - pnorm(split + R$g1 - lam))
    pl <- ifelse(is.na(L$g1), 0, pnorm(split - L$g1 - lam) - pnorm(split - L$g2 - lam))
    mean(pr + pl)
  }
  tq_ <- qt(1 - AL / 2, df = nu - 1); qa <- tq_ * sqrt(nu / (nu - 1 + tq_^2))
  res <- list()
  for (vv in versions) {
    ## oracle threshold of this version: delta0 known, |Q(delta0)| = q_alpha, B^2/A = q^2 S0/nu
    cth <- vapply(qa, function(q) lbf_oracle(q^2 * S0 / nu, S0, nu, vv), 0)
    Rs <- side(vv, +1, cth); Ls <- side(vv, -1, cth)
    for (k in seq_along(AL)) {
      res[[length(res) + 1]] <- data.frame(tag, method = paste0("mbf_", vv), alpha = AL[k],
        size_ratio = prob(Rs[[k]], Ls[[k]], 0) / AL[k], power = prob(Rs[[k]], Ls[[k]], tq_[k]),
        lambda = tq_[k], nonmonotone_share = mean(c(Rs[[k]]$interior, Ls[[k]]$interior)),
        threshold = cth[k])
    }
  }
  list(res = do.call(rbind, res), mbf_u = mbf_u)
}

## ---- same configurations and seeds as 06 / 06b -----------------------------------
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

## ---- diagnostic: log MBF against u on four ERAP2 draws, three versions --------------
diag_erap2 <- function() {
  L <- load_geu(1L); sp <- spectral_setup(L$Z, j = L$j)
  out <- list()
  for (h2 in c(0.05, 0.20)) {
    set.seed(5)
    ob <- mbf_config(sp, h2 / (1 - h2), data.frame(config = "diag"), B = 4L, versions = character(0))
    ## fixed-tau2 MBF as in 06 (tau2 = TAU2), for comparison
    fx <- function(u) {
      e <- environment(ob$mbf_u)
      bv <- u / sqrt(e$S0); sk <- e$Rk + outer(e$xi, bv); s2p <- e$Rp + e$xp * bv
      AQ <- crossprod(sk * sk, e$WQ) + (s2p^2 + e$ssr); BQ <- crossprod(e$xi * sk, e$WQ) + e$xp * s2p
      l0 <- -0.5 * sweep(e$nu * log(AQ), 2, e$ldQ, "+"); vQ <- TAU2 * e$SQg
      lbf <- sweep(-(e$nu / 2) * log1p(-sweep(BQ^2 / AQ, 2, TAU2 / (1 + vQ), "*")), 2,
                   0.5 * log1p(vQ), "-")
      lw <- sweep(l0, 2, e$lwq, "+"); rowLSE(lw + lbf) - rowLSE(lw)
    }
    for (u0 in c(0, 5, 10, 20, 40, 80, 120, 300, 1000)) {
      uu <- rep(u0, 4)
      out[[length(out) + 1]] <- data.frame(h2 = h2, u = u0,
        fixed = median(fx(uu)), hc = median(ob$mbf_u(uu, "hc")), hg = median(ob$mbf_u(uu, "hg")))
    }
  }
  d <- do.call(rbind, out)
  write.csv(d, "calibration/results/mbf_prior_diag.csv", row.names = FALSE)
  print(d, digits = 3)
}
diag_erap2()
if (DIAG_ONLY) quit(save = "no")

OUT <- "calibration/results/mbf_prior_fix"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", n = z$n, blocks = z$nb,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2))
  } else {
    id <- sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", n = nrow(L$Z), blocks = NA,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(7000 + z$i + round(100 * z$h2))
  }
  out <- mbf_config(sp, z$h2 / (1 - z$h2), tag)$res
  saveRDS(out, ck)
  cat(sprintf("[%7.1f s] %s done\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), id))
  ck
}
st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e)
  paste("ERROR", conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st)[grepl("ERROR", unlist(st))])
fx <- do.call(rbind, lapply(list.files(OUT, "rds$", full.names = TRUE), readRDS))
write.csv(fx, "calibration/results/mbf_prior_fix.csv", row.names = FALSE)
cat("configs done:", length(unique(fx$config)), "of", length(cfgs), "\n")
cat("done in", format(Sys.time() - t0), "\n")
