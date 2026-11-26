# ==============================================================================
# 20_nu_vs_n.R
# How does the effective df nu_j scale with n at a fixed cis window?
# Phased GEUVADIS haplotypes of two windows (RPS26, PEX6 region; MAF >= 5% in the 358
# EUR samples) are resampled to n in {87, 179, 358, 1000, 3000, 10000} individuals
# (pairs of haplotypes drawn with replacement), which keeps the LD of the window.
# For 30 fixed candidates per window, nu_j is computed at delta* = h2 / (1 - h2),
# h2 in {0.10, 0.25}, with the candidate-specific kernel K_-j = Z_-j Z_-j' / (m - 1)
# and F = intercept. The spectrum comes from the (m-1) x (m-1) Gram matrix, so large n
# costs one crossprod. Also reported: the Satterthwaite df nubar_j (1/nubar = 1/(nu-1)
# + 1/nu_j) and the Student tail ratio P(|t_nubar| > c) / P(|t_{nu-1}| > c) at the
# two-sided 5e-8 threshold, an order of magnitude of the plug-in distortion (note, Prop. 4).
# Usage (from the repository root): Rscript calibration/20_nu_vs_n.R [ncores]
# Output: results/nu_vs_n.csv
# ==============================================================================
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(TRUE); NCORES <- if (length(args)) as.integer(args[1]) else 4L
NS <- c(87L, 179L, 358L, 1000L, 3000L, 10000L); H2 <- c(0.10, 0.25); NCAND <- 30L
GENES <- c(RPS26 = "ENSG00000197728", `PEX6-region` = "ENSG00000124587")

load_haps <- function(tag) {                                # 2n x m haplotype matrix, MAF filter as in 14
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(tag, ".gt"))), "\t", fixed = TRUE)
  g <- lapply(L, function(r) r[-1L])
  h1 <- vapply(g, function(v) as.integer(substr(v, 1, 1)), integer(length(g[[1]])))
  h2 <- vapply(g, function(v) as.integer(substr(v, 3, 3)), integer(length(g[[1]])))
  stopifnot(all(substr(g[[1]], 2, 2) == "|"))
  G <- h1 + h2; maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  rbind(h1[, keep], h2[, keep])
}

nu_index <- function(ell, xi, xi0sq, nu, delta) {           # index nu_j with a zero block
  w <- 1 / (1 + delta * ell); a <- ell * w
  S <- sum(w * xi^2) + xi0sq; p <- w * xi^2 / S
  aU <- sum(a) / nu; VU <- sum(a^2) / nu - aU^2; ap <- sum(p * a)
  nu * VU / (ap - aU)^2
}

one_config <- function(gene, n) {
  H <- HAPS[[gene]]; m <- ncol(H)
  set.seed(1000 + n + match(gene, names(GENES)))
  X <- H[sample.int(nrow(H), n, TRUE), ] + H[sample.int(nrow(H), n, TRUE), ]
  ok <- apply(X, 2, sd) > 0
  Z <- scale(X); Z[, !ok] <- 0
  G <- crossprod(Z)                                         # m x m, Z'Z
  nu <- n - 1
  do.call(rbind, lapply(CANDS[[gene]], function(j) {
    if (!ok[j]) return(NULL)
    keep <- setdiff(which(ok), j); e <- eigen(G[keep, keep] / (length(keep)), symmetric = TRUE)
    pos <- e$values > 1e-8 * e$values[1]; ell <- e$values[pos]
    xi <- as.numeric(G[j, keep] %*% e$vectors[, pos]) / sqrt(length(keep) * ell)
    xi0sq <- max(G[j, j] - sum(xi^2), 0)
    do.call(rbind, lapply(H2, function(h2) {
      nuj <- nu_index(ell, xi, xi0sq, nu, h2 / (1 - h2))
      nubar <- 1 / (1 / (nu - 1) + 1 / nuj); cc <- qt(1 - 2.5e-8, nu - 1)
      data.frame(gene = gene, n = n, m = m, h2 = h2, cand = j, nu = nu, nu_j = nuj, nubar = nubar,
                 tail_ratio = pt(cc, nubar, lower.tail = FALSE) / pt(cc, nu - 1, lower.tail = FALSE))
    }))
  }))
}

HAPS <- lapply(GENES, load_haps)
CANDS <- lapply(HAPS, function(H) { set.seed(7); sort(sample.int(ncol(H), NCAND)) })
cfg <- expand.grid(gene = names(GENES), n = NS, stringsAsFactors = FALSE)
cfg <- cfg[order(-cfg$n), ]                                 # large n first (longest jobs)
res <- do.call(rbind, mclapply(seq_len(nrow(cfg)), function(i) one_config(cfg$gene[i], cfg$n[i]),
                               mc.cores = NCORES, mc.preschedule = FALSE))
write.csv(res, "calibration/results/nu_vs_n.csv", row.names = FALSE)
q <- function(v) sprintf("%.0f [%.0f-%.0f]", median(v), quantile(v, .25), quantile(v, .75))
S <- do.call(rbind, lapply(split(res, list(res$gene, res$h2, res$n), drop = TRUE), function(d)
  data.frame(gene = d$gene[1], h2 = d$h2[1], n = d$n[1], nu_j = q(d$nu_j), p_lt30 = mean(d$nu_j < 30),
             p_lt100 = mean(d$nu_j < 100), tail_med = median(d$tail_ratio), tail_max = max(d$tail_ratio))))
print(S[order(S$gene, S$h2, S$n), ], row.names = FALSE, digits = 3)
