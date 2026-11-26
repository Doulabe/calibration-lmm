# ==============================================================================
# 03_plugin_exact_calibration.R
#
# The exact law of the Wald statistic depends on the unknown delta0. In practice
# it is evaluated at delta_hat ("plug-in exact law"). Is the resulting test
# calibrated?  The event {p_exact(delta_hat) < alpha} equals
# {|t| > c_alpha(delta_hat)}, where c_alpha(delta) is the upper alpha/2 point of
# the exact law at delta. Since delta_hat is a function of the error contrasts r
# only, the conditional normal representation still applies:
#   P(reject | r) = Phi((m - c)/sqrt v) + Phi((-c - m)/sqrt v),  c = c_alpha(delta_hat(r)).
# c_alpha(delta) is tabulated on a grid of delta by inner Rao-Blackwellised Monte
# Carlo and interpolated on log(delta).
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/03_plugin_exact_calibration.R
# Output: calibration/results/plugin_exact_calibration.csv
# ==============================================================================

B_DRAW <- 4000L
source("calibration/00_functions.R")
B_IN  <- 4000L; B_OUT <- 20000L
CGRID <- seq(0, 40, by = 0.02)

## per-draw (m, v, delta_hat) under true delta0 (same algebra as exact_size)
draws <- function(sp, delta0, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu0 <- sp$nu0; r <- sp$r
  Wg  <- 1 / (1 + outer(ell, GRID)); ldg <- colSums(log(1 + outer(ell, GRID)))
  S2g <- colSums(xi^2 * Wg) + xp^2
  d0  <- 1 + delta0 * ell
  res <- list(); nch <- ceiling(B / CHUNK)
  for (ch in seq_len(nch)) {
    b <- CHUNK
    s <- matrix(rnorm(r * b), r, b) * sqrt(d0)
    sp_ <- if (sp$nperp) rnorm(b) else numeric(b)
    ssr <- if (sp$nrest > 0) rchisq(b, sp$nrest) else numeric(b)
    s2 <- s * s
    Ag <- crossprod(s2, Wg) + (sp_^2 + ssr)
    Bg <- crossprod(xi * s, Wg) + xp * sp_
    RSSg <- Ag - sweep(Bg^2, 2, S2g, "/")
    O <- -0.5 * (sweep((nu0 - 1) * log(RSSg), 2, ldg + log(S2g), "+"))
    g <- max.col(O, ties.method = "first"); dh <- GRID[g]
    int <- g > 2 & g < length(GRID)
    if (any(int)) {
      ii <- which(int); lg <- log(GRID)
      x1 <- lg[g[ii] - 1]; x2 <- lg[g[ii]]; x3 <- lg[g[ii] + 1]
      f1 <- O[cbind(ii, g[ii] - 1)]; f2 <- O[cbind(ii, g[ii])]; f3 <- O[cbind(ii, g[ii] + 1)]
      den <- (x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1)
      num <- (x2 - x1)^2 * (f2 - f3) - (x2 - x3)^2 * (f2 - f1)
      xs <- ifelse(abs(den) > 1e-14, x2 - 0.5 * num / den, x2)
      dh[ii] <- exp(pmin(pmax(xs, x1), x3))
    }
    at <- function(dl) {
      iD <- 1 / (1 + outer(ell, dl))
      A  <- colSums(s2 * iD) + sp_^2 + ssr
      Bx <- colSums(xi * s * iD) + xp * sp_
      S2 <- colSums(xi^2 * iD) + xp^2
      list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
    }
    Q <- at(dh); Q0 <- at(rep(delta0, b))
    bh <- Q$Bx / Q$S2; sig2 <- Q$RSS / (nu0 - 1)
    m <- (bh - Q0$Bx / Q0$S2) * sqrt(Q$S2 / sig2)
    v <- Q$S2 / (Q0$S2 * sig2)
    res[[ch]] <- data.frame(m = m, v = v, dh = dh)
  }
  do.call(rbind, res)
}

## exact upper points c_alpha(delta) for a vector of alphas
crit_points <- function(sp, delta, B = B_IN) {
  D <- draws(sp, delta, B)
  sd_ <- sqrt(D$v)
  G <- vapply(CGRID, function(c0)
    mean(pnorm((D$m - c0) / sd_) + pnorm((-c0 - D$m) / sd_)), 0)
  lg <- log(pmax(G, 1e-300))
  vapply(log(ALPHA), function(la) approx(rev(lg), rev(CGRID), xout = la, ties = mean)$y, 0)
}

calib <- function(sp, delta0, tag) {
  dgrid <- c(0, delta0 * exp(seq(log(1 / 200), log(200), length.out = 41)))
  CT <- t(vapply(dgrid, function(dl) crit_points(sp, dl), numeric(length(ALPHA))))
  D  <- draws(sp, delta0, B_OUT)
  ## interpolate c_alpha at delta_hat (linear in delta near 0, in log(delta) above)
  cfun <- function(k, dh) {
    y <- CT[, k]; lpos <- log(dgrid[-1])
    out <- numeric(length(dh))
    lo <- dh <= dgrid[2]
    out[lo] <- y[1] + (y[2] - y[1]) * dh[lo] / dgrid[2]
    out[!lo] <- approx(lpos, y[-1], xout = log(pmin(dh[!lo], max(dgrid))), rule = 2)$y
    out
  }
  cc_naive <- qt(1 - ALPHA / 2, df = sp$nu0 - 1)
  sd_ <- sqrt(D$v)
  out <- c(tag, list(delta0 = delta0))
  for (k in seq_along(ALPHA)) {
    nm <- format(ALPHA[k], scientific = TRUE)
    ca <- cfun(k, D$dh)
    pe <- pnorm((D$m - ca) / sd_) + pnorm((-ca - D$m) / sd_)
    pn <- pnorm((D$m - cc_naive[k]) / sd_) + pnorm((-cc_naive[k] - D$m) / sd_)
    out[[paste0("naive_ratio_", nm)]] <- mean(pn) / ALPHA[k]
    out[[paste0("plugin_exact_ratio_", nm)]] <- mean(pe) / ALPHA[k]
    out[[paste0("plugin_exact_rse_", nm)]] <- sd(pe) / sqrt(length(pe)) / mean(pe)
    out[[paste0("c_naive_", nm)]] <- cc_naive[k]
    out[[paste0("c_exact_at_delta0_", nm)]] <- approx(log(dgrid[-1]), CT[-1, k],
                                                      xout = log(delta0))$y
  }
  out
}

res <- list(); t0 <- Sys.time()
## synthetic: n = 1000, 5 blocks (rho 0.95), candidate inside block 1
set.seed(1000 * 1000 + 5)
Z <- gen_blocks(1000L, 5L, 50L, 0.95)
sp <- spectral_setup(Z, j = 25L)
for (h2 in c(0.05, 0.20)) {
  set.seed(31 + round(100 * h2))
  res[[length(res) + 1]] <- calib(sp, h2 / (1 - h2),
    list(design = "AR1_rho0.95", n = 1000L, blocks = 5L, h2 = h2, candidate = "in_LD"))
  cat(sprintf("[%6.1f s] synthetic h2=%.2f done\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), h2))
}
## GEUVADIS leads: ERAP2 and SNHG5
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
G2 <- data.frame(tag = c("ENSG00000164308", "ENSG00000203875"), symbol = c("ERAP2", "SNHG5"),
                 lead_pos = c(96252589L, 86387888L), stringsAsFactors = FALSE)
for (i in seq_len(nrow(G2))) {
  g <- G2[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t",
                fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  sp <- spectral_setup(Z, j = match(g$lead_pos, pos))
  for (h2 in c(0.05, 0.20)) {
    set.seed(57 + i + round(100 * h2))
    res[[length(res) + 1]] <- calib(sp, h2 / (1 - h2),
      list(design = "GEUVADIS_cis", n = nrow(Z), blocks = NA, h2 = h2,
           candidate = paste0("lead_", g$symbol)))
    cat(sprintf("[%6.1f s] GEUVADIS %s h2=%.2f done\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs")), g$symbol, h2))
  }
}
out <- do.call(rbind, lapply(res, as.data.frame))
write.csv(out, "calibration/results/plugin_exact_calibration.csv", row.names = FALSE)
cat("done in", format(Sys.time() - t0), "\n")
