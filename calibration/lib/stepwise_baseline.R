# ==============================================================================
#  stepwise_baseline.R
#  Stepwise baselines with a FIXED regional kernel (MLMM / StepLMM scheme,
#  EMMAX-type conditional scan, MLMe-type candidate exclusion) and conditional
#  OLS stepwise regression, all with a profile-ML eBIC stopping rule (or a
#  Bonferroni rule). The scripts 14, 14b, 14c, 27 and 27c use
#  stepwise_baseline() (conditional OLS stepwise regression).
#  Requires bf_stepwise_exact.R (exact_pieces_at_delta, neg2_profile_ml) to be
#  sourced first.
#
#  This file should NOT be invoked directly with Rscript; it defines functions
#  and exits.  No top-level side effects.
# ==============================================================================

# ── REML estimate of delta for a FIXED kernel (eigen-coordinates) ─────────────
reml_delta_fixed <- function(Uty, UtF, s_K, n) {
  p <- ncol(UtF)
  f <- function(ld) {
    delta <- exp(ld); d_inv <- 1 / (1 + delta * s_K)
    M  <- crossprod(UtF, d_inv * UtF)
    pv <- as.numeric(crossprod(UtF, d_inv * Uty))
    ch <- tryCatch(chol(M + 1e-12 * diag(p)), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    sol <- backsolve(ch, forwardsolve(t(ch), pv))
    rss <- sum(Uty^2 * d_inv) - sum(pv * sol)
    if (!is.finite(rss) || rss <= 0) return(1e10)
    sum(log1p(delta * s_K)) + 2 * sum(log(diag(ch))) + (n - p) * log(rss)
  }
  opt <- optimize(f, c(log(1e-6), log(1e6)))
  if (f(log(1e-6)) <= opt$objective) 1e-6 else exp(opt$minimum)
}

# ── conditional GLS statistics of every candidate at fixed delta, fixed kernel ─
cond_stats_fixed <- function(Uty, UtF, UtX, s_K, delta, n) {
  k <- ncol(UtF); d_inv <- 1 / (1 + delta * s_K)
  M0 <- crossprod(UtF, d_inv * UtF)
  p0 <- as.numeric(crossprod(UtF, d_inv * Uty))
  L0 <- chol(M0 + 1e-12 * diag(k))
  w0 <- backsolve(L0, forwardsolve(t(L0), p0))
  RSS0 <- sum(Uty^2 * d_inv) - sum(p0 * w0)
  alpha <- colSums(d_inv * UtX^2)
  beta  <- as.numeric(crossprod(UtX, d_inv * Uty))
  FWX   <- crossprod(UtF, d_inv * UtX)
  H     <- backsolve(L0, forwardsolve(t(L0), FWX))
  D     <- pmax(alpha - colSums(FWX * H), 1e-300)
  u     <- beta - as.numeric(crossprod(FWX, w0))
  list(Q2 = u^2 / (D * RSS0 / (n - k)), D = D, u = u, RSS0 = RSS0)
}

# ── precomputed rotation shared by the kernel arms of one dataset ─────────────
precompute_rotation <- function(X, y, eig = NULL) {
  if (is.null(eig)) eig <- eigen(tcrossprod(X), symmetric = TRUE)
  U <- eig$vectors
  list(U = U, s_raw = pmax(eig$values, 0),
       Uty = as.numeric(crossprod(U, y)), UtX = crossprod(U, X),
       UtA = crossprod(U, matrix(1, nrow(X), 1L)))
}

# ── generic forward path with eBIC (or Bonferroni) stopping ───────────────────
# kernel: "fixed" (K = XX'/m, all markers), "excl" (all markers but the
#         candidate, selected markers inside), "none" (OLS)
# delta_mode: "refit" (REML at every step), "once" (intercept-only REML, held)
# stop: "ebic" or "bonf"
stepwise_baseline <- function(y, X, kernel = c("fixed", "excl", "none"),
                              delta_mode = c("refit", "once"),
                              stop = c("ebic", "bonf"),
                              K_max = 10L, gamma = 1, alpha_bonf = 0.05,
                              eig = NULL, pre = NULL) {
  kernel <- match.arg(kernel); delta_mode <- match.arg(delta_mode)
  stop <- match.arg(stop)
  n <- nrow(X); m <- ncol(X); a_n <- matrix(1, n, 1L)
  pen_unit <- log(n) + 2 * gamma * log(m)
  t0 <- Sys.time()

  if (kernel == "none") {
    U <- NULL; Uty <- y; UtX <- X; UtA <- a_n
    s_full <- rep(0, n); s_raw <- rep(0, n)
  } else {
    if (is.null(pre)) pre <- precompute_rotation(X, y, eig)
    U <- pre$U; s_raw <- pre$s_raw; s_full <- s_raw / m
    Uty <- pre$Uty; UtX <- pre$UtX; UtA <- pre$UtA
  }
  rot <- function(Fm) if (is.null(U)) Fm else crossprod(U, Fm)

  selected <- integer(0); remaining <- seq_len(m)
  stat_path <- numeric(0); delta_path <- numeric(0); ebic_path <- numeric(0)

  delta0 <- if (kernel == "none") 0 else reml_delta_fixed(Uty, UtA, s_full, n)
  ebic_prev <- if (stop == "ebic") {
    if (kernel == "none") n * log(sum((y - mean(y))^2) / n)
    else neg2_profile_ml(y, a_n, U, s_full)
  } else Inf
  ebic_path <- c(ebic_path, ebic_prev)

  for (k_step in seq_len(K_max)) {
    if (length(remaining) <= 1L) break
    F_k <- if (length(selected) == 0L) a_n else cbind(a_n, X[, selected, drop = FALSE])
    UtF <- rot(F_k); k <- ncol(F_k)
    delta <- if (kernel == "none") 0
             else if (delta_mode == "once") delta0
             else reml_delta_fixed(Uty, UtF, s_full, n)
    if (kernel == "excl") {
      pc <- exact_pieces_at_delta(Uty, UtF, UtX[, remaining, drop = FALSE],
                                  s_rank = s_raw / (m - 1), delta = delta,
                                  a = delta / (m - 1), tau2 = 0.04, n = n, k = k)
      Q2 <- pc$u^2 / (pc$D * pc$RSS0 / (n - k)); D <- pc$D; u <- pc$u; RSS0 <- pc$RSS0
    } else {
      cs <- cond_stats_fixed(Uty, UtF, UtX[, remaining, drop = FALSE],
                             s_full, delta, n)
      Q2 <- cs$Q2; D <- cs$D; u <- cs$u; RSS0 <- cs$RSS0
    }
    best_idx <- which.max(Q2); best_j <- remaining[best_idx]
    if (!is.finite(Q2[best_idx])) break
    S_prop <- c(selected, best_j)

    if (stop == "ebic") {
      F_prop <- cbind(a_n, X[, S_prop, drop = FALSE])
      ll <- if (kernel == "none") {
        r <- lm.fit(F_prop, y)$residuals; n * log(sum(r^2) / n)
      } else neg2_profile_ml(y, F_prop, U, s_full)
      ebic_new <- ll + length(S_prop) * pen_unit
      if (ebic_new >= ebic_prev) break
    } else {
      RSS1 <- max(RSS0 - u[best_idx]^2 / D[best_idx], 1e-300)
      Fst  <- (RSS0 - RSS1) / (RSS1 / (n - k - 1))
      pval <- pf(Fst, 1, n - k - 1, lower.tail = FALSE)
      if (pval > alpha_bonf / m) break
      ebic_new <- NA_real_
    }
    selected <- S_prop; remaining <- setdiff(remaining, best_j)
    stat_path <- c(stat_path, Q2[best_idx]); delta_path <- c(delta_path, delta)
    ebic_path <- c(ebic_path, ebic_new); ebic_prev <- ebic_new
  }
  list(indices = selected, K_hat = length(selected), Q2_at_step = stat_path,
       delta_path = delta_path, ebic_path = ebic_path,
       elapsed = as.numeric(Sys.time() - t0, units = "secs"))
}

# ── the five baseline arms ────────────────────────────────────────────────────
ARMS <- list(
  FK_refit = list(kernel = "fixed", delta_mode = "refit", stop = "ebic"),
  FK_once  = list(kernel = "fixed", delta_mode = "once",  stop = "ebic"),
  FK_excl  = list(kernel = "excl",  delta_mode = "refit", stop = "ebic"),
  OLS_eBIC = list(kernel = "none",  delta_mode = "once",  stop = "ebic"),
  OLS_Bonf = list(kernel = "none",  delta_mode = "once",  stop = "bonf"))

run_baseline_arms <- function(y, X, truth, K_max = 10L, gamma = 1, eig = NULL) {
  pre <- precompute_rotation(X, y, eig)
  res <- lapply(names(ARMS), function(a) {
    p <- ARMS[[a]]
    fit <- stepwise_baseline(y, X, kernel = p$kernel, delta_mode = p$delta_mode,
                             stop = p$stop, K_max = K_max, gamma = gamma,
                             pre = if (p$kernel == "none") NULL else pre)
    c(list(arm = a, indices = fit$indices, elapsed = fit$elapsed,
           delta_path = fit$delta_path), as.list(metr(fit$indices, truth)))
  })
  names(res) <- names(ARMS); res
}

# ── metrics ───────────────────────────────────────────────────────────────────
metr <- function(idx, truth) {
  tp <- length(intersect(idx, truth)); K <- length(idx)
  prec <- if (K > 0) tp / K else 0; rec <- tp / length(truth)
  f1 <- if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0
  c(K_hat = K, tp = tp, precision = prec, recall = rec, f1 = f1)
}

# block-level metrics: a selection counts as a block hit when it lies in the LD
# block (block_size consecutive markers) of a causal variant; "proxy swaps" are
# selected non-causal variants inside a causal block
metr_block <- function(idx, truth, block_size = 10L) {
  blk <- function(j) (j - 1L) %/% block_size + 1L
  cb <- blk(truth); sb <- if (length(idx)) blk(idx) else integer(0)
  hits <- length(intersect(cb, sb))
  c(block_recall = hits / length(cb),
    block_precision = if (length(idx)) mean(sb %in% cb) else 0,
    proxy_swaps = sum(!(idx %in% truth) & sb %in% cb),
    off_block_fp = sum(!(sb %in% cb)))
}
