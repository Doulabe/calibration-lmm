# ==============================================================================
# 02_score_null_exact.R
#
# Exact size of the SCORE-type conditional statistic when the variance ratio is
# estimated under the candidate-specific NULL model (EMMAX / MLMM
# style: delta_hat_0 = argmax of the restricted likelihood without the candidate
# as fixed effect), in the same spectral reduction as 01_exact_size_pilot.R.
#
# Here delta_hat_0 depends on the candidate direction, so t | r is not Gaussian.
# Exact computation: s = r + xi b with b = beta_hat(delta0) ~ N(0, 1/S2(delta0))
# independent of the error contrasts r. For each draw of r, T_s is a function of
# the scalar u = b sqrt(S2(delta0)) ~ N(0,1); the crossing points u+ and u- of
# the cut-off are found by vectorised bisection and
#   P(|T_s| > c | r) = Phi(-u+) + Phi(u-)      (monotonicity checked).
# Oracle check: the same machinery with delta fixed at delta0 must reproduce the
# analytic size of the score statistic, P(|t_{nu0-1}| > c sqrt((nu0-1)/(nu0-c^2))).
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/02_score_null_exact.R [B=2000]
# Output: calibration/results/score_null_exact.csv
# ==============================================================================

args   <- commandArgs(TRUE)
B_DRAW <- if (length(args) >= 1) as.integer(args[1]) else 2000L
source("calibration/00_functions.R")
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
NBIS  <- 26L

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

score_exact_size <- function(sp, delta0, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu0 <- sp$nu0; r <- sp$r
  cc  <- qt(1 - ALPHA / 2, df = nu0 - 1)
  Wg  <- 1 / (1 + outer(ell, GRID2)); ldg <- colSums(log(1 + outer(ell, GRID2)))
  d0  <- 1 + delta0 * ell
  S20 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S20
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b              # error contrasts
  Tfun <- function(u, fixed = FALSE) {                      # u: length-B vector
    bv <- u / sqrt(S20)
    sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv
    if (fixed) dh <- rep(delta0, B) else {
      A <- crossprod(sk * sk, Wg) + (s2p^2 + ssr)
      dh <- pick_grid(-0.5 * sweep(nu0 * log(A), 2, ldg, "+"), GRID2)
    }
    iD <- 1 / (1 + outer(ell, dh))
    Ah <- colSums(sk * sk * iD) + s2p^2 + ssr
    Bx <- colSums(xi * sk * iD) + xp * s2p
    S2 <- colSums(xi^2 * iD) + xp^2
    list(T = (Bx / sqrt(S2)) / sqrt(Ah / nu0), dh = dh)
  }
  root <- function(c0, side, fixed) {                       # side = +1 upper, -1 lower
    lo <- rep(-12, B); hi <- rep(80, B)
    if (side < 0) { lo <- rep(-80, B); hi <- rep(12, B) }
    f <- function(u) side * Tfun(u, fixed)$T - c0
    flo <- f(if (side > 0) lo else hi); fhi <- f(if (side > 0) hi else lo)
    # for side>0 we need f(lo)<0<f(hi); for side<0 we search u- with -T(u-)=c0
    a <- if (side > 0) lo else hi; z <- if (side > 0) hi else lo
    never <- fhi < 0; always <- flo > 0
    for (it in seq_len(NBIS)) {
      mid <- 0.5 * (a + z); fm <- f(mid)
      up <- fm < 0; a[up] <- mid[up]; z[!up] <- mid[!up]
    }
    u <- 0.5 * (a + z); u[never] <- side * Inf; u[always] <- -side * Inf
    u
  }
  prob <- function(fixed) {
    sapply(cc, function(c0) {
      up <- root(c0, +1, fixed); dn <- root(c0, -1, fixed)
      pnorm(-up) + pnorm(dn)
    })
  }
  pr  <- prob(FALSE); po <- prob(TRUE)
  ## monotonicity diagnostic on the first 60 draws
  ug  <- seq(-12, 12, by = 0.5); nb <- min(60L, B)
  Tm  <- sapply(ug, function(u0) Tfun(rep(u0, B))$T[seq_len(nb)])
  mono <- mean(apply(Tm, 1, function(v) all(diff(v) > 0)))
  dh0 <- Tfun(rep(0, B))$dh
  analytic <- 2 * pt(-cc * sqrt((nu0 - 1) / (nu0 - cc^2)), df = nu0 - 1)
  out <- c(monotone_frac = mono, p_dhat0_null = mean(dh0 < 1e-4))
  for (k in seq_along(ALPHA)) {
    nm <- format(ALPHA[k], scientific = TRUE)
    out[paste0("score_ratio_", nm)]  <- mean(pr[, k]) / ALPHA[k]
    out[paste0("score_rse_", nm)]    <- sd(pr[, k]) / sqrt(B) / mean(pr[, k])
    out[paste0("oracle_score_", nm)] <- mean(po[, k]) / ALPHA[k]
    out[paste0("analytic_score_", nm)] <- analytic[k] / ALPHA[k]
  }
  out
}

res <- list(); t0 <- Sys.time()
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) {
  set.seed(1000 * n + nb)
  Z  <- gen_blocks(n, nb, 50L, 0.95)
  sp <- spectral_setup(Z, j = 25L)
  for (h2 in c(0.05, 0.20)) {
    d0 <- h2 / (1 - h2)
    ex <- score_exact_size(sp, d0)
    res[[length(res) + 1]] <- c(list(design = "AR1_rho0.95", n = n, blocks = nb,
      h2 = h2, candidate = "in_LD", r = sp$r), as.list(ex))
    cat(sprintf("[%6.1f s] n=%4d blocks=%3d h2=%.2f done\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs")), n, nb, h2))
  }
}
## GEUVADIS leads
D <- "data/geuvadis"; LOC <- file.path(D, "loci")
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
for (i in seq_len(nrow(GENES))) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path(LOC, paste0(g$tag, ".gt"))), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  jl <- match(g$lead_pos, pos); if (is.na(jl)) next
  set.seed(77 + i)
  sp <- spectral_setup(Z, j = jl)
  for (h2 in c(0.05, 0.20)) {
    ex <- score_exact_size(sp, h2 / (1 - h2))
    res[[length(res) + 1]] <- c(list(design = "GEUVADIS_cis", n = nrow(Z), blocks = NA,
      h2 = h2, candidate = paste0("lead_", g$symbol), r = sp$r), as.list(ex))
  }
  cat(sprintf("[%6.1f s] GEUVADIS %s done\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), g$symbol))
}
out <- do.call(rbind, lapply(res, as.data.frame))
write.csv(out, "calibration/results/score_null_exact.csv", row.names = FALSE)
cat("done in", format(Sys.time() - t0), "\n")
