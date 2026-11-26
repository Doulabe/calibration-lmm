# ==============================================================================
# 57b_nonnormal_table.R
# Supplementary table of the sensitivity to non-Gaussian errors (57_nonnormal_errors.R): size
# relative to the nominal level by direct simulation, for Gaussian, Student t5 and centred
# gamma errors, at 1e-3 and 1e-4, with the largest binomial standard error of each row.
# Usage (from the repository root): [ART=<article folder>] Rscript calibration/57b_nonnormal_table.R
# Output: <ART>/tables/tabS_nonnormal.tex
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/note")
D <- read.csv("calibration/results/nonnormal_errors.csv", stringsAsFactors = FALSE)
ERR <- c("normal", "t5", "gamma2"); LEV <- c(1e-3, 1e-4)
TEST <- c(known = "Known $\\delta$", wald = "Wald, $\\widehat\\delta$ with the candidate", null = "Null $\\widehat\\delta_0$, $t(\\widehat\\delta_0)$",
          kr = "Kenward--Roger", lrt = "Likelihood ratio", rlr = "Restricted likelihood ratio")
SYM <- c(ERAP2 = "\\textit{ERAP2}", SNHG5 = "\\textit{SNHG5}")
fm <- function(v) if (!is.finite(v)) "--" else if (v == 0) "$0$" else if (v >= 100) sprintf("$%.0f$", v) else
  if (v >= 10) sprintf("$%.1f$", v) else if (v >= 0.01) sprintf("$%.2f$", v) else {
    e <- floor(log10(v)); sprintf("$%.1f\\times10^{%d}$", v / 10^e, e) }
cf <- unique(D[order(D$nu_delta), c("config", "nu_delta", "n")])
lab <- function(id, nu, n) {
  if (grepl("^syn", id)) return(sprintf("AR(1), $n=%d$, five blocks, $\\nu_\\delta=%.1f$", n, nu))
  p <- strsplit(id, "_")[[1]]
  sprintf("%s, position %s, $n=%d$, $\\nu_\\delta=%.1f$", SYM[[p[2]]], gsub(",", "\\\\,", formatC(as.numeric(p[3]), format = "d", big.mark = ",")), n, nu)
}
body <- unlist(lapply(seq_len(nrow(cf)), function(q) {
  d <- D[D$config == cf$config[q], ]
  c(if (q > 1) "\\midrule", sprintf("\\multicolumn{8}{l}{\\emph{%s}} \\\\", lab(cf$config[q], cf$nu_delta[q], cf$n[q])),
    vapply(names(TEST), function(tt) {
      z <- d[d$test == tt, ]
      v <- unlist(lapply(LEV, function(a) vapply(ERR, function(e) { w <- z[z$error == e & abs(z$alpha - a) < a * 1e-6, ]
        if (nrow(w)) fm(w$size_ratio) else "--" }, "")))
      w <- z[z$alpha %in% LEV, ]; rse <- max(ifelse(w$rejections > 0, w$se / w$size_ratio, NA), na.rm = TRUE)
      sprintf("\\quad %s & %s & %s \\\\", TEST[[tt]], paste(v, collapse = " & "), if (is.finite(rse)) sprintf("$%.2f$", rse) else "--") }, "")) }))
dir.create(file.path(ART, "tables"), showWarnings = FALSE, recursive = TRUE)
writeLines(c("\\begin{tabular}{lccccccc}", "\\toprule",
             "& \\multicolumn{3}{c}{$\\alpha=10^{-3}$} & \\multicolumn{3}{c}{$\\alpha=10^{-4}$} & \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
             "Test & Gaussian & $t_5$ & Gamma & Gaussian & $t_5$ & Gamma & Rel. s.e. \\\\", "\\midrule",
             body, "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabS_nonnormal.tex"))
cat("wrote", file.path(ART, "tables", "tabS_nonnormal.tex"), "\n")
## summary for the text: ratio of the non-Gaussian to the Gaussian size, by test
G <- merge(D[D$error != "normal", ], D[D$error == "normal", c("config", "test", "alpha", "size_ratio")], by = c("config", "test", "alpha"), suffixes = c("", "_normal"))
G$rel <- G$size_ratio / G$size_ratio_normal
options(width = 200)
print(G[G$alpha <= 1e-3 & is.finite(G$rel), c("config", "error", "test", "alpha", "rejections", "size_ratio", "size_ratio_normal", "rel")], digits = 3, row.names = FALSE)
