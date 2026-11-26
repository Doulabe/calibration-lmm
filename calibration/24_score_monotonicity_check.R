# ==============================================================================
# 24_score_monotonicity_check.R
# Numerical support for the proof.
# For each draw r of the error contrasts, the score statistic at the null-model
# REML estimate, Q(u; r) = Q(delta_hat0(u); r + xi sigma u / sqrt(S0)), is evaluated
# on a fine u-grid. Reported per configuration:
#   - share of draws where u -> Q is not increasing (a decrease larger than 1e-8),
#   - share of those decreases that coincide with a jump of delta_hat0(u),
#   - exact size at alpha = 1e-3, 1e-5, 5e-8 computed (i) with the two-crossing
#     formula (bisection, as in 02 and 06) and (ii) with the general
#     formula, integrating phi over the union of rejection intervals found on the grid.
# Usage (from the repository root): nice -n 15 Rscript calibration/24_score_monotonicity_check.R [cores=3]
# Output: results/score_monotonicity_check.csv
# ==============================================================================
args <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
B_DRAW <- 300L
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
AL <- c(1e-3, 1e-5, 5e-8); NBIS <- 34L
UG <- seq(-40, 40, by = 0.1)

pick_grid <- function(O, grid) {                          # as in 06
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

check_config <- function(sp, delta0, tag, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  Qf <- function(u) {                                     # u: length-B vector
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    A2 <- crossprod(sk2, W2) + (s2p^2 + ssr)
    dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
    iD <- 1 / (1 + outer(ell, dh0))
    Ah <- colSums(sk2 * iD) + s2p^2 + ssr; Bh <- colSums(xi * sk * iD) + xp * s2p
    Sh <- colSums(xi^2 * iD) + xp^2
    list(Q = Bh / sqrt(Sh * Ah / nu), dh0 = dh0)
  }
  QM <- matrix(0, B, length(UG)); DM <- QM
  for (g in seq_along(UG)) { z <- Qf(rep(UG[g], B)); QM[, g] <- z$Q; DM[, g] <- z$dh0 }
  dQ <- QM[, -1] - QM[, -ncol(QM)]
  nonmono <- apply(dQ < -1e-8, 1, any)
  jump <- abs(log(pmax(DM[, -1], 1e-6)) - log(pmax(DM[, -ncol(DM)], 1e-6))) > 0.5
  dec_with_jump <- if (any(dQ < -1e-8)) mean(jump[dQ < -1e-8]) else NA
  max_drop <- if (any(nonmono)) max(-dQ[dQ < 0]) else 0
  tq <- qt(1 - AL / 2, df = nu - 1)
  ## (i) two-crossing formula with bisection (under monotonicity)
  cross <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- Qf(mid)$Q < level
      lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    (lo + hi) / 2 }
  size_mono <- vapply(tq, function(c0) mean(pnorm(-cross(c0)) + pnorm(cross(-c0))), 0)
  ## (ii) general formula: union of rejection intervals on the grid (linear interpolation),
  ## tails beyond the grid added with the sign of Q at the grid ends
  size_gen <- vapply(tq, function(c0) {
    mean(vapply(seq_len(B), function(i) {
      f <- abs(QM[i, ]) - c0; rejg <- f > 0
      p <- 0
      if (rejg[1]) p <- p + pnorm(UG[1])
      if (rejg[length(UG)]) p <- p + pnorm(-UG[length(UG)])
      k <- which(diff(sign(f)) != 0)                       # crossings between UG[k] and UG[k+1]
      xs <- UG[k] - f[k] * (UG[k + 1] - UG[k]) / (f[k + 1] - f[k])
      pts <- sort(c(UG[1], xs, UG[length(UG)]))
      for (m in seq_len(length(pts) - 1)) {
        mid <- (pts[m] + pts[m + 1]) / 2
        if (abs(approx(UG, QM[i, ], xout = mid)$y) > c0) p <- p + pnorm(pts[m + 1]) - pnorm(pts[m])
      }
      p }, 0)) }, 0)
  data.frame(tag, B = B, share_nonmonotone = mean(nonmono), dec_with_dhat0_jump = dec_with_jump,
             max_drop = max_drop, alpha = AL, size_ratio_twocross = size_mono / AL,
             size_ratio_general = size_gen / AL, rel_diff = size_gen / size_mono - 1)
}

GENES <- data.frame(tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000174652"),
                    symbol = c("ERAP2", "RPS26", "ZNF266"),
                    lead_pos = c(96252589L, 56401085L, 9544276L), stringsAsFactors = FALSE)
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
cfgs <- c(lapply(seq_len(nrow(GENES)), function(i) list(kind = "geu", i = i, h2 = 0.05)),
          list(list(kind = "geu", i = 1L, h2 = 0.20)),
          list(list(kind = "syn", n = 200L, nb = 5L, h2 = 0.05)),
          list(list(kind = "syn", n = 1000L, nb = 20L, h2 = 0.20)))
run_one <- function(z) {
  if (z$kind == "syn") {
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2))
  } else {
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2))
  }
  tag$nu_delta <- indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]]
  set.seed(9100 + round(100 * z$h2) + nchar(tag$config))
  check_config(sp, z$h2 / (1 - z$h2), tag)
}
res <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) conditionMessage(e)),
                mc.cores = NCORES, mc.preschedule = FALSE)
bad <- vapply(res, is.character, TRUE); if (any(bad)) print(unlist(res[bad]))
out <- do.call(rbind, res[!bad])
write.csv(out, "calibration/results/score_monotonicity_check.csv", row.names = FALSE)
print(out, row.names = FALSE, digits = 3)
