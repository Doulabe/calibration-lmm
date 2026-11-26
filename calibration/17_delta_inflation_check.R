# ==============================================================================
# 17_delta_inflation_check.R
# Mechanism check for the selection-level demonstration: on the datasets of
# 14_selection_demo.R (regenerated from the same seeds), estimate delta in the
# model the path starts from (intercept only, kernel of the whole window) and in
# the model where the causal variants are fixed effects (kernel without them).
# Causal variants not yet selected are absorbed by the kernel and should inflate
# delta_hat towards (h2bg + Kc h2c) / sigma2_e instead of delta* = h2bg / sigma2_e.
# Usage (from the repository root): Rscript calibration/17_delta_inflation_check.R [ncores]
# Output: results/selection_delta_inflation.csv
# ==============================================================================
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(TRUE); NCORES <- if (length(args)) as.integer(args[1]) else 4L
B <- 25L
GR <- c(0, 10^seq(-4, 3, length.out = 281))

## ML and REML estimates of delta on the grid, V = sigma2 (I + delta K), fixed effects F
delta_hat <- function(y, F, ev, U) {
  n <- length(y); k <- ncol(F); Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F)
  ll <- vapply(GR, function(d) {
    w <- 1 / (1 + d * ev); M <- crossprod(UtF, w * UtF); p <- crossprod(UtF, w * Uty)
    R <- sum(w * Uty^2) - sum(p * solve(M, p)); ldV <- -sum(log(w))
    c(ml = -0.5 * ldV - 0.5 * n * log(R),
      reml = -0.5 * ldV - 0.5 * as.numeric(determinant(M)$modulus) - 0.5 * (n - k) * log(R)) }, c(0, 0))
  c(ml = GR[which.max(ll["ml", ])], reml = GR[which.max(ll["reml", ])])
}

GENES <- data.frame(tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000124587"),
                    symbol = c("ERAP2", "RPS26", "SLFN5", "PEX6-region"), stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_Z <- function(tag) {
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(tag, ".gt"))), "\t", fixed = TRUE)
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  Z <- scale(G[, keep, drop = FALSE]); Z[!is.finite(Z)] <- 0; Z
}
ZS <- lapply(GENES$tag, load_Z)
EG <- lapply(ZS, function(Z) eigen(tcrossprod(Z) / ncol(Z), symmetric = TRUE))
grid <- rbind(
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = c(1L, 3L), h2c = c(0.08, 0.15),
              b = seq_len(B), stringsAsFactors = FALSE))

one_rep <- function(i) {
  g <- grid[i, ]; Z <- ZS[[g$gene]]; n <- nrow(Z); m <- ncol(Z)
  set.seed(30000 + 1000 * g$gene + 100 * g$Kc + round(100 * g$h2bg) + 7919 * g$b + round(1000 * g$h2c))
  causal <- integer(0)
  while (length(causal) < g$Kc) {
    cand <- sample.int(m, 1)
    if (!length(causal) || max(cor(Z[, cand], Z[, causal, drop = FALSE])^2) < 0.5) causal <- c(causal, cand)
  }
  bvec <- rnorm(m); u <- as.numeric(Z %*% bvec); u <- u / sd(u) * sqrt(g$h2bg)
  sig <- as.numeric(Z[, causal, drop = FALSE] %*% rep(sqrt(g$h2c), g$Kc))
  s2e <- max(1 - g$h2bg - g$Kc * g$h2c, 0.05)
  y <- sig + u + rnorm(n) * sqrt(s2e)
  one <- matrix(1, n, 1)
  d0 <- delta_hat(y, one, pmax(EG[[g$gene]]$values, 0), EG[[g$gene]]$vectors)
  if (g$Kc > 0) {
    Zr <- Z[, -causal, drop = FALSE]; e <- eigen(tcrossprod(Zr) / ncol(Zr), symmetric = TRUE)
    dc <- delta_hat(y, cbind(one, Z[, causal, drop = FALSE]), pmax(e$values, 0), e$vectors)
  } else dc <- d0
  data.frame(gene = GENES$symbol[g$gene], h2bg = g$h2bg, Kc = g$Kc, h2c = g$h2c, b = g$b,
             delta_star = g$h2bg / s2e, delta_absorb = (g$h2bg + g$Kc * g$h2c) / s2e,
             d0_ml = d0[["ml"]], d0_reml = d0[["reml"]], dc_ml = dc[["ml"]], dc_reml = dc[["reml"]])
}

res <- do.call(rbind, mclapply(seq_len(nrow(grid)), one_rep, mc.cores = NCORES))
write.csv(res, "calibration/results/selection_delta_inflation.csv", row.names = FALSE)
med <- function(v) round(median(v), 2)
S <- aggregate(cbind(r0_reml = d0_reml / delta_star, r0_ml = d0_ml / delta_star,
                     rc_reml = dc_reml / delta_star, absorb = delta_absorb / delta_star) ~ h2bg + Kc + h2c,
               res, med)
print(S, row.names = FALSE)
