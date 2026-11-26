# ==============================================================================
# 42b_misspecified_table.R
# Supplementary table of the misspecified-covariance study (42_misspecified_covariance.R):
# exact size relative to the nominal level, median and range over the 24 candidates, for
# each scenario and treatment. Ratios below 1e-6 are censored as in Table 2 and Table S5.
# Usage (from the repository root):
#   ART=<article folder> Rscript calibration/42b_misspecified_table.R
# Output: calibration/output/tables/tabA_misspecified.tex
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/output")
M <- read.csv("calibration/results/misspecified.csv")
fm <- function(v) { if (!is.finite(v)) return("--"); if (v < 1e-6) return("$<10^{-6}$")
  if (v < 0.01 || v >= 1e4) { e <- floor(log10(v)); m <- v / 10^e; if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    return(if (abs(round(m, 1) - 1) < 1e-9) sprintf("$10^{%d}$", e) else sprintf("$%.1f\\times10^{%d}$", m, e)) }
  if (v >= 1000) return(sprintf("$%s$", formatC(round(v), format = "d", big.mark = "{,}")))
  if (v >= 10) sprintf("$%.0f$", v) else sprintf("$%.2f$", v) }
fr <- function(v) sprintf("%s {\\scriptsize [%s, %s]}", fm(median(v)), fm(min(v)), fm(max(v)))
SC <- c(correct = "Regional kernel (correct)", sparse = "Sparse, $20$ markers of the window",
        outside = "Half outside the window")
MT <- c(oracle_star = "Wald at the pseudo-true $\\delta^\\ast$", wald = "Wald, $\\widehat\\delta$ with the candidate",
        emmax = "Null estimate, $t(\\widehat\\delta_0)$", lrt = "Likelihood ratio, $\\chi^2_1$", kr = "Kenward--Roger")
body <- unlist(lapply(names(SC), function(s) c(if (s != "correct") "\\midrule",
  sprintf("\\multicolumn{4}{l}{\\emph{%s}} \\\\", SC[[s]]),
  vapply(names(MT), function(mt) { z <- M[M$scen == s & M$method == mt, ]
    sprintf("\\quad %s & %s \\\\", MT[[mt]], paste(vapply(c(1e-3, 1e-5, 5e-8), function(a) fr(z$size_ratio[abs(z$alpha - a) < a * 1e-6]), ""),
                                                  collapse = " & ")) }, ""))))
writeLines(c("\\begin{tabular}{p{5.2cm}ccc}", "\\toprule",
             "Background and test & $\\alpha=10^{-3}$ & $\\alpha=10^{-5}$ & $\\alpha=5\\cdot10^{-8}$ \\\\", "\\midrule",
             body, "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabA_misspecified.tex"))
cat("wrote tabA_misspecified.tex\n")
