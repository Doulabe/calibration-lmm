# ==============================================================================
# 33j_boundary_check.R
# Boundary solutions of the estimates of delta, which are maximized over a bounded domain
# D = [0, delta_max] because the maximum-likelihood profile is unbounded (33h).
#  (a) Observed data, paths of the article (27c): the decisions use the grid on [0, 1e3]
#      (GRIDD) and the reported p-values the grid on [0, 1e4] (GRID2); the statistics of the
#      selected variants are compared between the two.
#  (b) Draws of the exact distributions, decisive GEUVADIS cases of 33g: share of the draws in
#      which the REML estimate under H1 (Wald, Kenward-Roger), the ML estimate under H0 (score,
#      likelihood ratio) or the ML estimate under H1 (likelihood ratio) sits at the upper end of
#      D, at u = 0, 3 and 6, for
#        - D = [0, 1e4] (calibration study, coarse evaluation) at delta0 = 0.053 and 0.25, the two
#          ratios of the calibration study (h2 = 0.05 and 0.20);
#        - D = GEMMA's range [0, 1e-5 ... 1e5] (refined evaluation) at delta0 = 1, 1e2, 1e4, 1e6, 1e10.
# Usage (from the repository root): Rscript calibration/33j_boundary_check.R
# Output: results/boundary_check.csv (and printed summaries)
# ==============================================================================
Sys.setenv(PC_SOURCE_ONLY = "1")
try(source("calibration/33g_protected_check.R"), silent = TRUE)
stopifnot(exists("setup_case"), exists("gemma_stats"))
GRID_GEMMA <- GRID2; GRID_CALIB <- c(0, 10^seq(-4, 4, length.out = 81))

## ---- (a) observed statistics of the paths of the article, [0, 1e3] versus [0, 1e4] -------------
cmp <- list()
for (f in list.files("calibration/results/geuvadis_app_cov", "^(EUR|YRI)_ENSG[0-9]+[.]rds$", full.names = TRUE)) {
  r <- readRDS(f); if (is.null(r$variants)) next
  for (tr in names(r$rec)) { rc <- r$rec[[tr]]; rc <- rc[rc$accepted, , drop = FALSE]
    for (i in seq_len(nrow(rc))) { v <- r$variants[r$variants$j == rc$j[i] & r$variants$step == rc$step[i], ]
      if (nrow(v) != 1) next
      st <- switch(tr, wald_h1 = c(rc$t2w[i], v$t^2, rc$level_t2[i]), score_h0 = c(rc$t2s[i], v$t_null^2, rc$level_t2[i]),
                   lrt = c(rc$lrt[i], v$lrt, rc$pen[i]))
      cmp[[length(cmp) + 1]] <- data.frame(pop = r$pop, gene = r$gene, path = tr, step = rc$step[i],
        d_wald = rc$t2w[i] / v$t^2 - 1, d_null = rc$t2s[i] / v$t_null^2 - 1, d_lrt = rc$lrt[i] / v$lrt - 1,
        stat_decision = st[1], stat_reported = st[2], level = st[3], reported_passes = st[2] > st[3]) } }
}
CMP <- do.call(rbind, cmp)
cat("(a) selected variants of the paths of the article:", nrow(CMP), " max |relative difference| [0,1e3] vs [0,1e4]:",
    signif(max(abs(unlist(CMP[, c("d_wald", "d_null", "d_lrt")]))), 3), "\n")
print(signif(sapply(CMP[, c("d_wald", "d_null", "d_lrt")], function(x) max(abs(x))), 3))
cat("relative difference of the statistic of the path itself, by path:\n")
print(aggregate(cbind(d = abs(stat_decision / stat_reported - 1)) ~ path, CMP, max), digits = 3)
cat("selected variants whose reported statistic does not pass the level:", sum(!CMP$reported_passes), "\n")
print(CMP[!CMP$reported_passes, ], digits = 3, row.names = FALSE)

## ---- (b) share of the draws at the upper bound of D ---------------------------------------------
at_top <- function(sp, delta0, B = 4000L, us = c(0, 3, 6)) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r; n <- sp$n; top <- max(GRID2)
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0; Rk <- s - outer(xi, b); Rp <- spp - xp * b
  W2 <- 1 / (1 + outer(ell, GRID2)); ldm <- colSums(log(1 + outer(ell, GRID2))) - sp$ldF; S2g <- colSums(xi^2 * W2) + xp^2
  do.call(rbind, lapply(us, function(u) { bv <- rep(u / sqrt(S0), B); Sk <- Rk + outer(xi, bv); Sp <- Rp + xp * bv
    z <- gemma_stats(sp, Sk, Sp, ssr)
    A2 <- crossprod(Sk^2, W2) + (Sp^2 + ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
    RSS2 <- A2 - sweep(B2^2, 2, S2g, "/"); ml1 <- -0.5 * sweep(n * log(RSS2), 2, ldm, "+")
    data.frame(delta0 = delta0, u = u, reml_h1_top = mean(z$dh1 >= top * (1 - 1e-9)),
               ml_h0_top = mean(z$dm0 >= top * (1 - 1e-9)),
               ml_h1_top = mean(max.col(ml1, ties.method = "first") == length(GRID2))) }))
}
res <- list()
for (k in seq_len(nrow(CASES))) { cs <- CASES[k, ]
  for (dom in c("calib", "gemma")) {
    GRID2 <<- if (dom == "calib") GRID_CALIB else GRID_GEMMA
    set.seed(3500 + k); S <- setup_case(cs)
    d0s <- if (dom == "calib") c(0.053, 0.25) else c(1, 1e2, 1e4, 1e6, 1e10)
    for (d0 in d0s) res[[length(res) + 1]] <- cbind(case = cs$id, domain = dom, upper = max(GRID2), at_top(S$sp, d0)) } }
RES <- do.call(rbind, res)
write.csv(RES, "calibration/results/boundary_check.csv", row.names = FALSE)
cat("\n(b) share of draws with the estimate at the upper end of D\n")
print(aggregate(cbind(reml_h1_top, ml_h0_top, ml_h1_top) ~ domain + delta0, RES, max), digits = 3)
cat("\nby case, refined domain, delta0 >= 1e4:\n")
print(RES[RES$domain == "gemma" & RES$delta0 >= 1e4 & RES$u == 0, ], digits = 3, row.names = FALSE)
