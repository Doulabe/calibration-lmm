# ==============================================================================
# 14_selection_demo.R
#
# Selection-level demonstration: does the treatment of the variance ratio delta
# change the variants selected by stepwise fine-mapping with a regional background?
#
# Data: real genotypes of four GEUVADIS cis windows (n = 358); phenotypes simulated
# with a regional polygenic background (u = Z b, h2_bg in {0.10, 0.25}) and K_c in
# {0, 1, 3} causal variants (per-causal variance 0.08, pairwise r2 < 0.5).
#
# Procedures share the same stepwise path (candidate-specific kernel K_jk, exact
# rank-one pieces of the Bayes-factor stepwise engine) and the same nominal level per step,
# given by the eBIC penalty pen = log n + 2 log m. Only the treatment of delta changes:
#   wald_h1     rank and stop by t^2 at delta_hat under H1 (per-candidate REML, GEMMA-type)
#   score_h0    rank and stop by the score at delta_hat_0 under H0 (EMMAX / MLMM-type)
#   lrt         rank and stop by the profile-ML LRT, delta re-estimated under each model
#   mbf_fixed   rank by the marginal BF (tau2 = 0.04), stop by the LRT
#   mbf_hc      rank by the marginal BF with tau ~ half-Cauchy(0, 0.2), stop by the LRT
# t^2 and score thresholds are the t^2 at which the delta-known LRT equals pen.
# References: cbf_deployed (Bayes-factor stepwise engine: marginal BF + eBIC with the
# candidate in the null kernel), ols_ebic (no random effect), susie (top-PIP variant of
# each 95% credible set, L = 10).
# Metrics: LD-aware (tag = r2 >= 0.8 with a causal), off-target false selections
# (r2 < 0.25 with every causal), any selection under the global null.
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/14_selection_demo.R [cores=5] [--smoke]
# Output: results/selection_demo/<cell>_b<b>.rds, results/selection_demo_summary.csv
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 5L
SMOKE  <- "--smoke" %in% args
suppressPackageStartupMessages({ library(parallel); library(susieR) })
source("calibration/lib/bf_stepwise_exact.R")      # Bayes-factor stepwise engine
source("calibration/lib/stepwise_baseline.R")
stopifnot(exists("stepwise_baseline"))

B      <- if (SMOKE) 2L else 25L
K_MAX  <- 10L
GRIDD  <- c(0, 10^seq(-4, 3, length.out = 41))             # delta nodes (0 = no background)
TAU2   <- 0.04; TAU0 <- 0.2
T2N    <- 10^seq(-6, 4, length.out = 41)                    # tau2 nodes (half-Cauchy on tau)
lw_t2  <- { t <- log(T2N); dt <- t[2] - t[1]
  w <- log(c(dt / 2, rep(dt, length(t) - 2), dt / 2)) - log(pi * TAU0) - log1p(T2N / TAU0^2) + 0.5 * t
  w - (max(w) + log(sum(exp(w - max(w))))) }
lw_d   <- { t <- log(GRIDD[-1]); dt <- t[2] - t[1]          # delta prior weights (half-Cauchy sqrt delta)
  log(c(dt / 2, rep(dt, length(t) - 2), dt / 2)) + 0.5 * t - log(pi) - log1p(GRIDD[-1]) }
lse_rows <- function(M) { mx <- apply(M, 1, max); mx + log(rowSums(exp(M - mx))) }

## exact per-candidate pieces at one delta (engine algebra) + log-determinants
pieces_x <- function(Uty, UtF, UtX, s_rank, delta, a, n, k) {
  d_inv <- 1 / (1 + delta * s_rank); WU_y <- d_inv * Uty
  q0 <- sum(Uty * WU_y); p0 <- as.numeric(crossprod(UtF, WU_y))
  M0 <- crossprod(UtF, d_inv * UtF); L0 <- chol(M0 + 1e-12 * diag(k))
  w0 <- backsolve(L0, forwardsolve(t(L0), p0)); S00 <- sum(p0 * w0)
  logdetM0 <- 2 * sum(log(diag(L0))); logdetA <- sum(log1p(delta * s_rank))
  alpha <- colSums(d_inv * UtX^2); beta <- as.numeric(crossprod(UtX, WU_y))
  FWX <- crossprod(UtF, d_inv * UtX); H <- backsolve(L0, forwardsolve(t(L0), FWX))
  phi <- colSums(FWX * H); psi <- as.numeric(crossprod(FWX, w0))
  t_j <- pmax(1 - a * alpha, 1e-12); g_j <- a / t_j; denom <- 1 + g_j * phi
  yZy <- q0 + g_j * beta^2; pMp <- S00 + 2 * g_j * beta * psi + g_j^2 * beta^2 * phi
  hp <- psi + g_j * beta * phi; S0 <- pMp - g_j * hp^2 / denom
  RSS0 <- pmax(yZy - S0, 1e-300)
  D <- pmax(alpha / t_j - phi / (denom * t_j^2), 1e-300); u <- (beta - hp / denom) / t_j
  list(D = D, u = u, RSS0 = RSS0, ldZ = logdetA + log(t_j), ldM = logdetM0 + log(denom))
}

## all treatment statistics for the candidates of one step
step_stats <- function(y, Xc, F_k, U, s_raw, m_k) {
  n <- length(y); k <- ncol(F_k); nu <- n - k; mc <- ncol(Xc); G <- length(GRIDD)
  s_rank <- pmax(s_raw, 0) / m_k
  Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F_k); UtX <- crossprod(U, Xc)
  Dm <- um <- R0 <- lZ <- lM <- matrix(0, mc, G)
  for (g in seq_len(G)) {
    pc <- pieces_x(Uty, UtF, UtX, s_rank, GRIDD[g], GRIDD[g] / m_k, n, k)
    Dm[, g] <- pc$D; um[, g] <- pc$u; R0[, g] <- pc$RSS0; lZ[, g] <- pc$ldZ; lM[, g] <- pc$ldM
  }
  R1 <- pmax(R0 - um^2 / Dm, 1e-300)
  l0 <- -0.5 * (lZ + lM) - 0.5 * nu * log(R0)
  l1 <- -0.5 * (lZ + lM + log(Dm)) - 0.5 * (nu - 1) * log(R1)
  g1 <- max.col(l1, ties.method = "first"); g0 <- max.col(l0, ties.method = "first")
  ix1 <- cbind(seq_len(mc), g1); ix0 <- cbind(seq_len(mc), g0)
  t2w <- (nu - 1) * (um[ix1]^2 / Dm[ix1]) / R1[ix1]
  Q2 <- nu * um[ix0]^2 / (Dm[ix0] * R0[ix0]); t2s <- (nu - 1) * Q2 / pmax(nu - Q2, 1e-12)
  lrt <- apply(-lZ - n * log(R1), 1, max) - apply(-lZ - n * log(R0), 1, max)
  ## marginal BF, fixed tau2 and half-Cauchy tau (positive delta nodes only)
  pos <- -1L
  base <- -0.5 * (lZ + lM)[, pos]; v <- TAU2 * Dm[, pos]
  lp0 <- (l0[, pos]) + rep(lw_d, each = mc)
  lp1f <- base - 0.5 * nu * log(pmax(R0[, pos] - v * um[, pos]^2 / (Dm[, pos] * (1 + v)), 1e-300)) -
          0.5 * log1p(v) + rep(lw_d, each = mc)
  mbf_fixed <- lse_rows(lp1f) - lse_rows(lp0)
  lp1h <- matrix(-Inf, mc, G - 1)
  for (g in seq_len(G - 1)) {
    Dg <- Dm[, g + 1]; ug <- um[, g + 1]; Rg <- R0[, g + 1]
    vv <- outer(Dg, T2N)                                   # mc x T
    L <- -0.5 * log1p(vv) - 0.5 * nu * log(pmax(Rg - vv * (ug^2 / Dg) / (1 + vv), 1e-300))
    L <- sweep(L, 2, lw_t2, "+"); lp1h[, g] <- base[, g] + lse_rows(L)
  }
  mbf_hc <- lse_rows(lp1h + rep(lw_d, each = mc)) - lse_rows(lp0)
  list(t2w = t2w, t2s = t2s, lrt = lrt, mbf_fixed = mbf_fixed, mbf_hc = mbf_hc, nu = nu)
}

TREAT <- c("wald_h1", "score_h0", "lrt", "mbf_fixed", "mbf_hc")
run_treatments <- function(y, X) {
  n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m); a_n <- matrix(1, n, 1)
  cache <- new.env()
  get_step <- function(sel) {                              # memoised by selected set
    key <- paste0("k:", paste(sort(sel), collapse = ","))
    if (!is.null(cache[[key]])) return(cache[[key]])
    rem <- setdiff(seq_len(m), sel)
    raw <- tcrossprod(X[, rem, drop = FALSE]); eg <- eigen(raw, symmetric = TRUE)
    F_k <- cbind(a_n, X[, sel, drop = FALSE])
    st <- step_stats(y, X[, rem, drop = FALSE], F_k, eg$vectors, pmax(eg$values, 0), length(rem) - 1L)
    st$rem <- rem; cache[[key]] <- st; st
  }
  out <- list()
  for (tr in TREAT) {
    sel <- integer(0)
    for (k in seq_len(K_MAX)) {
      st <- get_step(sel); ct <- (st$nu - 1) * expm1(pen / n)
      rank <- switch(tr, wald_h1 = st$t2w, score_h0 = st$t2s, lrt = st$lrt,
                     mbf_fixed = st$mbf_fixed, mbf_hc = st$mbf_hc)
      b <- which.max(rank)
      ok <- switch(tr, wald_h1 = st$t2w[b] > ct, score_h0 = st$t2s[b] > ct, st$lrt[b] > pen)
      if (!isTRUE(ok)) break
      sel <- c(sel, st$rem[b])
    }
    out[[tr]] <- sel
  }
  out
}

metrics <- function(sel, causal, X) {
  kc <- length(causal); ks <- length(sel)
  if (ks == 0) return(c(K_hat = 0, tp = 0, fp08 = 0, fp025 = 0, recall = if (kc) 0 else NA,
                        precision = NA, f1 = if (kc) 0 else NA, any_off = 0, exact = 0))
  if (kc == 0) return(c(K_hat = ks, tp = 0, fp08 = ks, fp025 = ks, recall = NA, precision = 0,
                        f1 = NA, any_off = 1, exact = 0))
  r2 <- cor(X[, sel, drop = FALSE], X[, causal, drop = FALSE])^2    # ks x kc
  best <- apply(r2, 1, max)
  tagged <- apply(r2 >= 0.8, 2, any)
  tp <- sum(best >= 0.8); prec <- tp / ks; rec <- mean(tagged)
  c(K_hat = ks, tp = tp, fp08 = sum(best < 0.8), fp025 = sum(best < 0.25), recall = rec,
    precision = prec, f1 = if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0,
    any_off = as.numeric(any(best < 0.25)), exact = sum(sel %in% causal))
}

## ---- data ------------------------------------------------------------------------------
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
ZS <- lapply(GENES$tag, load_Z)                          # loaded once, shared by the forks
grid <- rbind(
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = c(1L, 3L), h2c = c(0.08, 0.15),
              b = seq_len(B), stringsAsFactors = FALSE))
grid$cell <- sprintf("%s_bg%.2f_K%d_c%.2f", GENES$symbol[grid$gene], grid$h2bg, grid$Kc, grid$h2c)
OUT <- "calibration/results/selection_demo"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(i) {
  g <- grid[i, ]; ck <- file.path(OUT, sprintf("%s_b%03d.rds", g$cell, g$b))
  if (file.exists(ck) && file.size(ck) > 0 && !SMOKE) return("cached")
  Z <- ZS[[g$gene]]; n <- nrow(Z); m <- ncol(Z)
  set.seed(30000 + 1000 * g$gene + 100 * g$Kc + round(100 * g$h2bg) + 7919 * g$b + round(1000 * g$h2c))
  causal <- integer(0)
  while (length(causal) < g$Kc) {                         # causals with pairwise r2 < 0.5
    cand <- sample.int(m, 1)
    if (!length(causal) || max(cor(Z[, cand], Z[, causal, drop = FALSE])^2) < 0.5) causal <- c(causal, cand)
  }
  bvec <- rnorm(m); u <- as.numeric(Z %*% bvec); u <- u / sd(u) * sqrt(g$h2bg)
  sig <- as.numeric(Z[, causal, drop = FALSE] %*% rep(sqrt(g$h2c), g$Kc))
  e <- rnorm(n) * sqrt(max(1 - g$h2bg - g$Kc * g$h2c, 0.05))
  y <- sig + u + e
  t0 <- Sys.time()
  sels <- run_treatments(y, Z)
  el_tr <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cbf <- tryCatch(bf_stepwise_exact(y, Z, tau2 = TAU2, K_max = K_MAX, n_nodes = 15L)$indices,
                  error = function(e) NA_integer_)
  ols <- tryCatch(stepwise_baseline(y, Z, kernel = "none", delta_mode = "once", stop = "ebic",
                                    K_max = K_MAX)$indices, error = function(e) NA_integer_)
  sus <- tryCatch({ f <- susie(Z, y, L = 10, verbose = FALSE)
    cs <- f$sets$cs; if (is.null(cs)) integer(0) else
      unique(vapply(cs, function(ix) ix[which.max(f$pip[ix])], 1L)) }, error = function(e) NA_integer_)
  sels <- c(sels, list(cbf_deployed = cbf, ols_ebic = ols, susie = sus))
  met <- do.call(rbind, lapply(names(sels), function(nm) {
    s <- sels[[nm]]; if (length(s) == 1 && is.na(s)) return(NULL)
    data.frame(method = nm, t(metrics(s, causal, Z)))
  }))
  res <- list(cell = g$cell, gene = GENES$symbol[g$gene], h2bg = g$h2bg, Kc = g$Kc, h2c = g$h2c, b = g$b,
              causal = causal, selections = sels, metrics = met, elapsed_treatments = el_tr)
  if (!SMOKE) saveRDS(res, ck)
  if (SMOKE) res else "done"
}

if (SMOKE) {
  for (i in c(which(grid$gene == 1 & grid$h2bg == 0.25 & grid$Kc == 3 & grid$b == 1),
              which(grid$gene == 1 & grid$h2bg == 0.25 & grid$Kc == 0 & grid$b == 1))) {
    t0 <- Sys.time(); r <- one_rep(i)
    cat(sprintf("\n%s b%d  causal: %s  (%.1f s total, %.1f s treatments)\n", r$cell, r$b,
                paste(r$causal, collapse = ","), as.numeric(difftime(Sys.time(), t0, units = "secs")),
                r$elapsed_treatments))
    for (nm in names(r$selections)) cat(sprintf("  %-13s %s\n", nm, paste(r$selections[[nm]], collapse = ",")))
    print(r$metrics, row.names = FALSE, digits = 3)
  }
} else {
  ord <- order(grid$b, grid$gene, grid$h2bg, grid$Kc, grid$h2c)
  st <- unlist(mclapply(ord, function(i) tryCatch(one_rep(i), error = function(e)
    paste("FAIL", grid$cell[i], grid$b[i], conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE))
  cat("replicates:", sum(st == "done"), "new,", sum(st == "cached"), "cached,", sum(grepl("FAIL", st)), "failed\n")
  print(head(st[grepl("FAIL", st)]))
  R <- do.call(rbind, lapply(list.files(OUT, "rds$", full.names = TRUE), function(f) {
    r <- readRDS(f); if (is.null(r$metrics)) return(NULL)
    cbind(data.frame(cell = r$cell, gene = r$gene, h2bg = r$h2bg, Kc = r$Kc, h2c = r$h2c, b = r$b), r$metrics) }))
  write.csv(R, "calibration/results/selection_demo_replicates.csv", row.names = FALSE)
  agg <- aggregate(cbind(K_hat, fp025, any_off, recall, precision, f1) ~ h2bg + Kc + h2c + method, R,
                   function(v) mean(v, na.rm = TRUE), na.action = na.pass)
  write.csv(agg, "calibration/results/selection_demo_summary.csv", row.names = FALSE)
  options(width = 200); print(agg, digits = 3, row.names = FALSE)
}
