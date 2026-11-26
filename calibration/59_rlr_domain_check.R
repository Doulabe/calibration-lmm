# ==============================================================================
# 59_rlr_domain_check.R
# Dependence of the likelihood ratio (ML) and of the restricted likelihood ratio on a common
# residual space on the upper end delta_max of the domain Omega, in the GEUVADIS application
# with GEMMA's kernel and fixed effects (33e, functions of 33g).
#  obs   observed statistics and nominal chi2_1 p-values of the two likelihood ratios for the
#        38 variant-step cases of GEMMA's paths, for delta_max in {1e3, 1e4, 1e5, 1e6, 1e8};
#  prot  protected p-values of the two statistics for the second selection of GEMMA's
#        likelihood-ratio path in the Yoruba PEX6 region, for delta_max in {1e4, 1e6, 1e8}:
#        true ratios 0 and 1e-2, ..., 1e10 by half-decades with B1 draws, then tenth-decades
#        within 0.3 decade of each maximum with B2 draws. The two statistics use the same draws,
#        and the same seeds are used for every delta_max.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/59_rlr_domain_check.R obs
#   nice -n 10 Rscript calibration/59_rlr_domain_check.R prot [cores=1]
# Output: results/rlr_domain_obs.csv; results/rlr_domain/<cap>.rds, results/rlr_domain_protected.csv
# ==============================================================================
a0 <- commandArgs(TRUE); MODE <- a0[1]; stopifnot(MODE %in% c("obs", "prot"))
NC <- if (length(a0) >= 2 && grepl("^[0-9]+$", a0[2])) as.integer(a0[2]) else 1L
commandArgs <- function(...) character(0)
Sys.setenv(PC_SOURCE_ONLY = "1")
try(source("calibration/33g_protected_check.R"), silent = TRUE)      # data, spectral_fixed, pick_grid
stopifnot(exists("spectral_fixed"), exists("pick_grid"))
suppressPackageStartupMessages(library(parallel))
RES <- "calibration/results"

## the two likelihood ratios and the Wald statistic from the coordinates of datasets (columns)
lr_stats <- function(sp, Sk, Sp, Ssr, grid) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- sp$n
  W2 <- 1 / (1 + outer(ell, grid)); ld2 <- colSums(log(1 + outer(ell, grid))); S2g <- colSums(xi^2 * W2) + xp^2
  A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
  lA <- log(A2); lR <- log(A2 - sweep(B2^2, 2, S2g, "/")); ldm <- ld2 - sp$ldF
  l1 <- -0.5 * sweep((nu - 1) * lR, 2, ld2 + log(S2g), "+")                 # restricted profile, candidate model
  l01 <- -0.5 * sweep((nu - 1) * lA, 2, ld2 + log(S2g), "+")                # null model on the same residual space
  ml1 <- -0.5 * sweep(n * lR, 2, ldm, "+"); ml0 <- -0.5 * sweep(n * lA, 2, ldm, "+")
  G <- length(grid); g <- function(O) max.col(O, ties.method = "first")
  list(lrt = 2 * (apply(ml1, 1, max) - apply(ml0, 1, max)), rlr = 2 * (apply(l1, 1, max) - apply(l01, 1, max)),
       sgn = B2[cbind(seq_len(nrow(B2)), g(l1))],                             # sign of the effect estimate
       end = cbind(ml0 = g(ml0) == G, ml1 = g(ml1) == G, r1 = g(l1) == G, r01 = g(l01) == G))
}
case_setup <- function(v, grid) {
  g <- GENES[GENES$symbol == v$gene, ]
  R <- readRDS(file.path(COMP, sprintf("%s_%s.rds", v$pop, g$tag)))
  w <- R$variants[R$variants$pos == v$pos & R$variants$step == v$step, ]; stopifnot(nrow(w) == 1)
  before <- NULL
  for (tn in strsplit(w$tests, "+", fixed = TRUE)[[1]]) { idx <- R$gemma_paths[[tn]]$idx
    if (length(idx) >= w$step && idx[w$step] == w$j) { before <- idx[seq_len(w$step - 1L)]; break } }
  stopifnot(!is.null(before))
  W <- load_window(v$pop, g$tag); X <- scale(W$G); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); C <- covariates(v$pop, W$samples)
  n <- nrow(X); m <- ncol(X); K <- tcrossprod(X) / m; stopifnot(W$pos[w$j] == v$pos)
  pen <- log(n) + 2 * log(m); nu_lev <- n - (1L + ncol(C)) - v$step + 1L
  GRID2 <<- grid                                              # global read by spectral_fixed
  sp <- spectral_fixed(X, y, C, before, w$j, K)
  list(sp = sp, level = log10(2 * pt(-sqrt((nu_lev - 1) * expm1(pen / n)), nu_lev - 1)))
}
V <- read.csv(file.path(RES, "geuvadis_cov_gemma_variants.csv"), stringsAsFactors = FALSE)
l10 <- function(x) pchisq(x, 1, lower.tail = FALSE, log.p = TRUE) / log(10)

if (MODE == "obs") {
  CAPS <- c(1e3, 1e4, 1e5, 1e6, 1e8); out <- list()
  for (i in seq_len(nrow(V))) for (cap in CAPS) {
    grid <- c(0, 10^seq(-5, log10(cap), by = 0.1)); S <- case_setup(V[i, ], grid)
    o <- lr_stats(S$sp, matrix(S$sp$obs$sk), S$sp$obs$sp, S$sp$obs$ssr, grid)
    out[[length(out) + 1]] <- data.frame(V[i, c("pop", "gene", "pos", "step", "tests")], cap = cap, level = S$level,
      lrt = o$lrt, rlr = o$rlr, log10p_lrt = l10(o$lrt), log10p_rlr = l10(o$rlr), o$end)
  }
  A <- do.call(rbind, out); write.csv(A, file.path(RES, "rlr_domain_obs.csv"), row.names = FALSE)
  rg <- function(v) max(v) - min(v)
  S <- aggregate(cbind(log10p_lrt, log10p_rlr) ~ pop + gene + pos + step, A, rg)
  cat("largest change of the nominal log10 p over delta_max: LR", max(S$log10p_lrt), " RLR", max(S$log10p_rlr), "\n")
  cat("maximizers at the upper end (any case, any delta_max):\n"); print(colSums(A[, c("ml0", "ml1", "r1", "r01")]))
  B <- A[A$cap == 1e5, ]; cat("at 1e5: RLR minus LR nominal log10 p, range", range(B$log10p_rlr - B$log10p_lrt), "\n")
  cat("decisions at the level that differ between LR and RLR (1e5):\n")
  print(B[(B$log10p_lrt < B$level) != (B$log10p_rlr < B$level), c("pop", "gene", "pos", "step", "tests", "level", "log10p_lrt", "log10p_rlr")], digits = 4)
  quit(save = "no")
}

## ---- protected p-values of the two statistics, common draws ---------------------------------------
B1 <- 2000L; B2 <- 6000L; NB <- 32L
D1 <- c(0, 10^seq(-2, 10, by = 0.5))
v <- V[V$pop == "YRI" & V$gene == "PEX6-region" & V$pos == 43011766L & V$step == 2L, ]; stopifnot(nrow(v) == 1)
OUT <- file.path(RES, "rlr_domain"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
mc <- function(lp) { m <- max(lp); w <- exp(lp - m); p <- mean(w); se <- sd(w) / sqrt(length(w))
  c(log10p = (m + log(p)) / log(10), log10up = (m + log(p + 2 * se)) / log(10),
    log10lo = (m + log(max(p - 2 * se, 1e-300 * p))) / log(10), rel_se = se / p) }
exact_two <- function(sp, grid, delta0, obs, B, seed) {
  set.seed(seed)
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); lr_stats(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr, grid) }
  lo <- rep(-60, B); hi <- rep(60, B)                         # zero of the effect estimate
  for (it in 1:NB) { mid <- (lo + hi) / 2; up <- ev(mid)$sgn < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  lsum <- function(a1, a2) pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))
  res <- list()
  for (key in c("lrt", "rlr")) { level <- obs[[key]]; f <- function(u) ev(u)[[key]]
    lo <- split; hi <- rep(250, B)
    for (it in 1:NB) { mid <- (lo + hi) / 2; up <- f(mid) < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    ur <- (lo + hi) / 2
    lo <- rep(-250, B); hi <- split
    for (it in 1:NB) { mid <- (lo + hi) / 2; dn <- f(mid) < level; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
    ul <- (lo + hi) / 2
    res[[key]] <- mc(lsum(pnorm(-ur, log.p = TRUE), pnorm(ul, log.p = TRUE))) }
  e0 <- ev(rnorm(B))$end                                      # null datasets: maximizers at the upper end
  list(p = res, end = colMeans(e0))
}
run_cap <- function(cap) {
  ck <- file.path(OUT, sprintf("cap_%g.rds", cap)); if (file.exists(ck)) return(readRDS(ck))
  grid <- c(0, 10^seq(-5, log10(cap), by = 0.1)); S <- case_setup(v, grid); sp <- S$sp
  o <- lr_stats(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr, grid); obs <- list(lrt = o$lrt, rlr = o$rlr)
  row <- function(d0, B, stage, k) { e <- exact_two(sp, grid, d0, obs, B, 590000 + 1000 * stage + k)
    do.call(rbind, lapply(names(e$p), function(key) data.frame(cap = cap, stat = key, stage = stage, delta0 = d0, B = B,
      t(e$p[[key]]), t(e$end)))) }
  P1 <- do.call(rbind, lapply(seq_along(D1), function(k) row(D1[k], B1, 1L, k)))
  P2 <- do.call(rbind, lapply(c("lrt", "rlr"), function(key) { q <- P1[P1$stat == key, ]
    a <- log10(pmax(q$delta0[which.max(q$log10p)], 1e-3)); D2 <- setdiff(signif(10^seq(a - 0.3, a + 0.3, by = 0.1), 6), signif(D1, 6))
    z <- do.call(rbind, lapply(seq_along(D2), function(k) row(D2[k], B2, 2L, 100L * (key == "rlr") + k))); z[z$stat == key, ] }))
  out <- list(profile = rbind(P1, P2), level = S$level, obs = obs, nominal = c(lrt = l10(o$lrt), rlr = l10(o$rlr)))
  saveRDS(out, ck); out
}
CAPS <- c(1e4, 1e6, 1e8)
R <- mclapply(CAPS, function(cp) tryCatch(run_cap(cp), error = function(e) conditionMessage(e)), mc.cores = NC, mc.preschedule = FALSE)
bad <- vapply(R, is.character, TRUE); if (any(bad)) { print(R[bad]); quit(save = "no", status = 1) }
P <- do.call(rbind, lapply(R, `[[`, "profile")); write.csv(P, file.path(RES, "rlr_domain_protected.csv"), row.names = FALSE)
S <- do.call(rbind, lapply(R, function(r) do.call(rbind, lapply(c("lrt", "rlr"), function(key) { q <- r$profile[r$profile$stat == key, ]
  i <- which.max(q$log10p); up <- max(q$log10up); lo <- max(q$log10lo)
  data.frame(cap = q$cap[1], stat = key, level = r$level, nominal = r$nominal[[key]], max_p = q$log10p[i], argmax = q$delta0[i],
             upper = up, lower = lo, class = if (up < r$level) "supported" else if (lo > r$level) "unsupported" else "undetermined",
             end_ml0 = max(q$ml0), end_ml1 = max(q$ml1), end_r1 = max(q$r1), end_r01 = max(q$r01)) }))))
options(width = 200); print(S[order(S$stat, S$cap), ], digits = 4, row.names = FALSE)
write.csv(S, file.path(RES, "rlr_domain_summary.csv"), row.names = FALSE)
