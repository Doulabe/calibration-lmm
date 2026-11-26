# ==============================================================================
# 33h_ml_profile_check.R
# The maximum-likelihood profile of delta under the null model, used by GEMMA's LRT and
# score tests, keeps the fixed-effect directions in log|V|. With the centred regional
# kernel, whose null space contains the intercept, it grows like (1/2) log(delta) as
# delta -> infinity (proof sketch: sum log(1+delta l_i) ~ nu log delta, -(n/2) log A ~
# (n/2) log delta, log|F'V^-1F| ~ -(p-1) log delta, n - nu = p), so the global maximizer
# is at infinity and the ML estimate is a local maximum on a bounded range (GEMMA:
# lambda in [1e-5, 1e5]; 33e: grid on [0, 1e4]). This script prints the profile on the
# observed data of the cases of 33g, and the location of the interior maximum.
# Usage (from the repository root): Rscript calibration/33h_ml_profile_check.R
# ==============================================================================
Sys.setenv(PC_SOURCE_ONLY = "1")
try(source("calibration/33g_protected_check.R"), silent = TRUE)
GRID2 <<- c(0, 10^seq(-5, 12, by = 0.1))
for (k in c(1, 3, 4)) {
  cs <- CASES[k, ]; S <- setup_case(cs); sp <- S$sp
  W2 <- 1 / (1 + outer(sp$ell, GRID2)); ld2 <- colSums(log(1 + outer(sp$ell, GRID2)))
  A2 <- colSums(sp$obs$sk^2 * W2) + sp$obs$sp^2 + sp$obs$ssr
  ml0 <- -0.5 * (sp$n * log(A2) + ld2 - sp$ldF)
  int <- GRID2 <= 1e5; i0 <- which.max(ml0[int])
  at <- function(d) ml0[which.min(abs(log10(pmax(GRID2, 1e-12)) - log10(d)))]
  cat(sprintf("%-16s n=%d nu=%d  interior max at delta=%.3g (ML %.2f) | ML at 1e5 %.2f, 1e8 %.2f, 1e12 %.2f | slope per decade beyond 1e8: %.3f (theory %.3f)\n",
      cs$id, sp$n, sp$nu0, GRID2[int][i0], ml0[int][i0], at(1e5), at(1e8), at(1e12), (at(1e12) - at(1e8)) / 4, 0.5 * log(10)))
}
