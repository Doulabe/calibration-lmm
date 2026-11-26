# ==============================================================================
# 52b_protected_window.R
# Position of the protected critical values of the likelihood ratio (52_protected_design.R)
# relative to the nominal chi2_1 quantiles. For each design and level, c* is the upper chi2_1
# quantile at some level alpha', reported as the ratio alpha/alpha'. Since c* lies between the
# quantiles at alpha and alpha/10, the protected and nominal decisions differ only for nominal
# p-values between alpha/10 and alpha (Discussion, Supplementary Section S2.3).
# Usage (from the repository root): Rscript calibration/52b_protected_window.R
# ==============================================================================
RES <- "calibration/results"
T <- read.csv(file.path(RES, "protected_design_thresholds.csv"))
E <- read.csv(file.path(RES, "protected_design_eval.csv"))
E <- E[E$kind == "config" & !is.na(E$config), ]
nu <- tapply(E$nu_delta, sub("_h[0-9.]+$", "", E$config), min)   # a design appears at one or two h2
T$nu <- nu[T$design]
T$ratio <- T$alpha / pchisq(T$lr, 1, lower.tail = FALSE)
cat(length(unique(T$design)), "designs\n")
for (a in sort(unique(T$alpha), decreasing = TRUE)) {
  d <- T[T$alpha == a, ]
  r <- tapply(d$ratio, cut(d$nu, c(0, 20, 100, Inf), right = FALSE,
                           labels = c("nu_delta < 20", "20 <= nu_delta < 100", "nu_delta >= 100")), max)
  cat(sprintf("alpha = %g: c* at the chi2_1 quantile of level alpha/%.2f to alpha/%.2f; largest ratio %s\n",
              a, max(d$ratio), min(d$ratio), paste(sprintf("%.2f (%s)", r, names(r)), collapse = ", ")))
}
cat(sprintf("all levels: alpha/%.2f to alpha/%.2f\n", max(T$ratio), min(T$ratio)))
