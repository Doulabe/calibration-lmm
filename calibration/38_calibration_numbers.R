# ==============================================================================
# 38_calibration_numbers.R
# Numbers quoted in the calibration section of the article, from the
# enlarged exact-calibration study (36_calibration_enlarged.R, results/calib_enlarged/).
# Prints, by class of nu_delta and nominal level, the median and range of the size
# ratio and of the power of each treatment, the share of configurations within a
# factor of 3 of nominal, and the same by source (GEUVADIS candidates, AR(1) designs,
# candidates unrelated to the region).
# Usage (from the repository root): Rscript calibration/38_calibration_numbers.R
# ==============================================================================
RES <- "calibration/results"
fs <- list.files(file.path(RES, "calib_enlarged"), "rds$", full.names = TRUE)
C <- do.call(rbind, lapply(fs, readRDS))
fs2 <- list.files(file.path(RES, "calib_enlarged_score"), "rds$", full.names = TRUE)
if (length(fs2)) { E <- do.call(rbind, lapply(fs2, readRDS)); C <- rbind(C, E[E$method == "emmax_t", names(C)]) }
C$src <- ifelse(grepl("^geu_", C$config), "geuvadis", ifelse(grepl("unlinked", C$config), "unlinked", "ar1"))
C$cls <- cut(C$nu_delta, c(0, 20, 100, Inf), right = FALSE, labels = c("nu<20", "20-100", ">=100"))
cat("configurations:", length(unique(C$config)), "\n")
U <- unique(C[, c("config", "src", "nu_delta", "cls")])
print(table(U$src, U$cls))
cat("nu range by source:\n"); print(tapply(U$nu_delta, U$src, range))
f <- function(v) sprintf("%.3g [%.3g, %.3g]", median(v), min(v), max(v))
for (a in c(1e-3, 1e-5, 5e-8)) {
  cat("\n==== alpha =", a, "====\n")
  d <- C[abs(C$alpha - a) < a * 1e-6, ]
  for (m in c("oracle_t", "wald_t", "emmax_t", "score_t", "lrt_chi2", "kr", "satterthwaite", "wald_edgeworth", "lrt_edgeworth")) {
    z <- d[d$method == m, ]
    s <- tapply(z$size_ratio, z$cls, f); p <- tapply(z$power, z$cls, function(v) sprintf("%.2g", median(v)))
    w3 <- tapply(z$size_ratio, z$cls, function(v) sprintf("%.0f%%", 100 * mean(v > 1 / 3 & v < 3)))
    cat(sprintf("%-15s size %s | within x3 %s | power %s\n", m, paste(s, collapse = " ; "), paste(w3, collapse = " "), paste(p, collapse = " ")))
  }
}
cat("\n==== by source at 5e-8 (size ratio) ====\n")
d <- C[abs(C$alpha - 5e-8) < 1e-13, ]
for (m in c("wald_t", "emmax_t", "score_t", "lrt_chi2", "kr")) { z <- d[d$method == m, ]
  cat(m, ":\n"); print(tapply(z$size_ratio, z$src, function(v) signif(c(med = median(v), min = min(v), max = max(v)), 3))) }
cat("\n==== AR(1) and unlinked configurations at 5e-8 ====\n")
z <- d[d$src != "geuvadis" & d$method %in% c("wald_t", "emmax_t", "lrt_chi2", "kr"), c("config", "nu_delta", "method", "size_ratio")]
print(reshape(z, idvar = c("config", "nu_delta"), timevar = "method", direction = "wide"), digits = 3)
cat("\n==== LR: worst configurations ====\n")
z <- C[C$method == "lrt_chi2", ]; print(head(z[order(-abs(log10(z$size_ratio))), c("config", "nu_delta", "alpha", "size_ratio")], 8), digits = 3)
cat("\n==== GEUVADIS leads only (published lead, 16 configurations) ====\n")
L <- C[C$role == "lead" & abs(C$alpha - 5e-8) < 1e-13 & C$method %in% c("wald_t", "emmax_t", "score_t", "lrt_chi2", "kr"), ]
print(tapply(L$size_ratio, L$method, function(v) signif(c(med = median(v), min = min(v), max = max(v)), 3)))
cat("nu_delta of the leads:", signif(range(unique(L[, c("config", "nu_delta")])$nu_delta), 3), "\n")
cat("\n==== Kenward-Roger degrees of freedom (median and 10th percentile over draws), by class ====\n")
X <- do.call(rbind, lapply(fs, function(f) attr(readRDS(f), "extra")))
X$cls <- cut(X$nu_delta, c(0, 20, 100, Inf), right = FALSE, labels = c("nu<20", "20-100", ">=100"))
print(aggregate(cbind(kr_df_med, kr_df_q10, p_boundary) ~ cls, X, function(v) signif(median(v), 3)))
cat("\n==== spread within classes at 5e-8 (orders of magnitude between min and max) ====\n")
z <- d[d$method %in% c("wald_t", "kr", "satterthwaite", "wald_edgeworth"), ]
print(aggregate(size_ratio ~ cls + method, z, function(v) round(log10(max(v)) - log10(max(min(v), 1e-300)), 1)))
