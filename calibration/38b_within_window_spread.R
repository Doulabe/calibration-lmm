# ==============================================================================
# 38b_within_window_spread.R
# Spread, between the ten candidates of one GEUVADIS window at one heritability, of the
# exact size ratio of the Wald test (delta re-estimated) at 5e-8, from the enlarged
# calibration study (results/calib_enlarged/). Quoted in the calibration section and the
# abstract of the article: 2.9 to 4.2 orders of magnitude, 3.5 in median.
# Usage (from the repository root): Rscript calibration/38b_within_window_spread.R
# ==============================================================================
RES <- "calibration/results"
fs <- list.files(file.path(RES, "calib_enlarged"), "rds$", full.names = TRUE)
C <- do.call(rbind, lapply(fs, readRDS))
d <- C[grepl("^geu_", C$config) & abs(C$alpha - 5e-8) < 1e-13 & C$method %in% c("wald_t", "lrt_chi2", "emmax_t"), ]
w <- d[d$method == "wald_t", ]
s <- do.call(rbind, lapply(split(w, list(w$gene, w$h2), drop = TRUE), function(z)
  data.frame(gene = z$gene[1], h2 = z$h2[1], k = nrow(z), nu_min = min(z$nu_delta), nu_max = max(z$nu_delta),
             wald_min = min(z$size_ratio), wald_max = max(z$size_ratio),
             spread_log10 = log10(max(z$size_ratio)) - log10(min(z$size_ratio)))))
print(s, digits = 3, row.names = FALSE)
cat("median within-window spread (log10):", round(median(s$spread_log10), 2), " range:", round(range(s$spread_log10), 2), "\n")
