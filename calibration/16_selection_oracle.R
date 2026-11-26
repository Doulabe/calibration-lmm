# ==============================================================================
# 16_selection_oracle.R
# Oracle arm of the selection-level demonstration (14_selection_demo.R): same
# datasets (regenerated from the same seeds, checked against the stored causal
# sets), same stepwise path and eBIC-equivalent level, but delta fixed at its
# design value delta* = h2bg / sigma2_e. Ranks and stops by the delta-known LRT
# n log(RSS0 / RSS1). Splits the recall gap between the mixed models and
# OLS / SuSiE into the cost of the kernel (delta known) and the cost of
# estimating delta.
# Usage (from the repository root): Rscript calibration/16_selection_oracle.R [ncores]
# Output: results/selection_oracle_replicates.csv
# ==============================================================================
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(TRUE); NCORES <- if (length(args) && !is.na(suppressWarnings(as.integer(args[1])))) as.integer(args[1]) else 4L
K_MAX <- 10L; B <- 25L

## engine algebra at one delta, identical to pieces_x() of 14_selection_demo.R
pieces_x <- function(Uty, UtF, UtX, s_rank, delta, a, n, k) {
  d_inv <- 1 / (1 + delta * s_rank); WU_y <- d_inv * Uty
  q0 <- sum(Uty * WU_y); p0 <- as.numeric(crossprod(UtF, WU_y))
  M0 <- crossprod(UtF, d_inv * UtF); L0 <- chol(M0 + 1e-12 * diag(k))
  w0 <- backsolve(L0, forwardsolve(t(L0), p0)); S00 <- sum(p0 * w0)
  alpha <- colSums(d_inv * UtX^2); beta <- as.numeric(crossprod(UtX, WU_y))
  FWX <- crossprod(UtF, d_inv * UtX); H <- backsolve(L0, forwardsolve(t(L0), FWX))
  phi <- colSums(FWX * H); psi <- as.numeric(crossprod(FWX, w0))
  t_j <- pmax(1 - a * alpha, 1e-12); g_j <- a / t_j; denom <- 1 + g_j * phi
  yZy <- q0 + g_j * beta^2; pMp <- S00 + 2 * g_j * beta * psi + g_j^2 * beta^2 * phi
  hp <- psi + g_j * beta * phi; S0 <- pMp - g_j * hp^2 / denom
  RSS0 <- pmax(yZy - S0, 1e-300)
  D <- pmax(alpha / t_j - phi / (denom * t_j^2), 1e-300); u <- (beta - hp / denom) / t_j
  list(D = D, u = u, RSS0 = RSS0)
}

oracle_path <- function(y, X, delta) {
  n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m); a_n <- matrix(1, n, 1)
  sel <- integer(0)
  for (k in seq_len(K_MAX)) {
    rem <- setdiff(seq_len(m), sel); m_k <- length(rem) - 1L
    eg <- eigen(tcrossprod(X[, rem, drop = FALSE]), symmetric = TRUE)
    F_k <- cbind(a_n, X[, sel, drop = FALSE]); U <- eg$vectors
    pc <- pieces_x(as.numeric(crossprod(U, y)), crossprod(U, F_k), crossprod(U, X[, rem, drop = FALSE]),
                   pmax(eg$values, 0) / m_k, delta, delta / m_k, n, ncol(F_k))
    lrt <- -n * log(pmax(1 - pc$u^2 / (pc$D * pc$RSS0), 1e-300))
    b <- which.max(lrt); if (!(lrt[b] > pen)) break
    sel <- c(sel, rem[b])
  }
  sel
}

metrics <- function(sel, causal, X) {                      # identical to 14_selection_demo.R
  kc <- length(causal); ks <- length(sel)
  if (ks == 0) return(c(K_hat = 0, tp = 0, fp08 = 0, fp025 = 0, recall = if (kc) 0 else NA,
                        precision = NA, f1 = if (kc) 0 else NA, any_off = 0, exact = 0))
  if (kc == 0) return(c(K_hat = ks, tp = 0, fp08 = ks, fp025 = ks, recall = NA, precision = 0,
                        f1 = NA, any_off = 1, exact = 0))
  r2 <- cor(X[, sel, drop = FALSE], X[, causal, drop = FALSE])^2
  best <- apply(r2, 1, max); tagged <- apply(r2 >= 0.8, 2, any)
  tp <- sum(best >= 0.8); prec <- tp / ks; rec <- mean(tagged)
  c(K_hat = ks, tp = tp, fp08 = sum(best < 0.8), fp025 = sum(best < 0.25), recall = rec,
    precision = prec, f1 = if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0,
    any_off = as.numeric(any(best < 0.25)), exact = sum(sel %in% causal))
}

## ---- data and design, identical to 14_selection_demo.R --------------------------------
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
grid <- rbind(
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = c(1L, 3L), h2c = c(0.08, 0.15),
              b = seq_len(B), stringsAsFactors = FALSE))
grid$cell <- sprintf("%s_bg%.2f_K%d_c%.2f", GENES$symbol[grid$gene], grid$h2bg, grid$Kc, grid$h2c)
IN <- "calibration/results/selection_demo"

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
  stored <- readRDS(file.path(IN, sprintf("%s_b%03d.rds", g$cell, g$b)))
  if (!identical(as.integer(stored$causal), as.integer(causal))) stop("causal set differs from 14")
  sel <- oracle_path(y, Z, g$h2bg / s2e)
  data.frame(cell = g$cell, gene = GENES$symbol[g$gene], h2bg = g$h2bg, Kc = g$Kc, h2c = g$h2c, b = g$b,
             method = "oracle_delta", t(metrics(sel, causal, Z)))
}

IDX <- if ("--smoke" %in% args) which(grid$b == 1 & grid$gene == 1 & grid$h2bg == 0.25) else seq_len(nrow(grid))
t0 <- Sys.time()
res <- mclapply(IDX, function(i) tryCatch(one_rep(i), error = function(e)
  paste("FAIL", grid$cell[i], grid$b[i], conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
bad <- vapply(res, is.character, TRUE)
cat("replicates:", sum(!bad), "ok,", sum(bad), "failed;", format(Sys.time() - t0), "\n")
if (any(bad)) print(head(unlist(res[bad])))
R <- do.call(rbind, res[!bad])
if (!"--smoke" %in% args) write.csv(R, "calibration/results/selection_oracle_replicates.csv", row.names = FALSE)
print(aggregate(cbind(K_hat, any_off, recall, precision) ~ h2bg + Kc + h2c, R,
                function(v) round(mean(v, na.rm = TRUE), 3), na.action = na.pass))
