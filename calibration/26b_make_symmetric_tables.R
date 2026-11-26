# ==============================================================================
# 26b_make_symmetric_tables.R
# Tables from 26_symmetric_check.R.
#   note/tab_laplace.tex     variation of 2 log MBF - RV along the effect u
#   note/tab_firstorder.tex  first-order predictions vs exact size ratios
# Usage (from the repository root): Rscript calibration/26b_make_symmetric_tables.R
# ==============================================================================
RES <- "calibration/results"; N <- "calibration/note"
L <- read.csv(file.path(RES, "symmetric_check_laplace.csv"))
lab <- function(cfg) {
  if (grepl("^geu_", cfg)) { p <- strsplit(sub("^geu_", "", cfg), "_h")[[1]]
    sprintf("GEUVADIS %s, $h^2=%s$", p[1], p[2]) } else {
    m <- regmatches(cfg, regexec("syn_n([0-9]+)_b([0-9]+)_h([0-9.]+)", cfg))[[1]]
    sprintf("AR(1), $n=%s$, %s blocs, $h^2=%s$", m[2], m[3], m[4]) } }
rows <- vapply(unique(L$config), function(cf) {
  d <- L[L$config == cf, ]; f <- d[grepl("fixe", d$prior), ]; h <- d[!grepl("fixe", d$prior), ]
  sprintf("%s & $%.1f$ & $%.0f$ & $%.1f$ & $%.0f$ & $%.1f$ & $%.0f$ \\\\", lab(cf), d$nu_delta[1],
          d$lr_range_u8_med[1], f$range_u8_med, f$range_u30_med, h$range_u8_med, h$range_u30_med) }, "")
writeLines(c("\\begin{tabular}{lrrrrrr}", "\\toprule",
  "& & RV & \\multicolumn{2}{c}{$\\tau^2$ fixé} & \\multicolumn{2}{c}{$\\tau$ demi-Cauchy} \\\\",
  "\\cmidrule(lr){3-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "Configuration & $\\nu_\\delta$ & $u\\le8$ & $u\\le8$ & $u\\le30$ & $u\\le8$ & $u\\le30$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_laplace.tex"))

F <- read.csv(file.path(RES, "symmetric_check_firstorder.csv"))
g <- function(v) { f1 <- function(x) if (x >= 100) sprintf("%.0f", x) else if (x >= 0.1) sprintf("%.2g", x) else {
  e <- floor(log10(x)); m <- round(x / 10^e); if (m == 10) { m <- 1; e <- e + 1 }
  if (m == 1) sprintf("10^{%d}", e) else sprintf("%d\\cdot10^{%d}", m, e) }
  sprintf("$%s$", f1(median(v))) }
rows <- vapply(sort(unique(F$alpha), decreasing = TRUE), function(a) { d <- F[abs(F$alpha - a) < a * 1e-6, ]
  sprintf("$%s$ & %s & %s & %s & %s & %s \\\\", if (a == 5e-8) "5\\cdot10^{-8}" else sprintf("10^{%d}", round(log10(a))),
          g(d$wald_pred), g(d$wald_exact), g(d$score_pred), g(d$score_exact), g(d$lrt_exact)) }, "")
writeLines(c("\\begin{tabular}{lccccc}", "\\toprule",
  "& \\multicolumn{2}{c}{Wald, $\\hat\\delta$} & \\multicolumn{2}{c}{Score, $\\hat\\delta_0$} & RV \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-6}",
  "$\\alpha$ & premier ordre & exact & premier ordre & exact & exact \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_firstorder.tex"))
cat(readLines(file.path(N, "tab_laplace.tex")), sep = "\n"); cat(readLines(file.path(N, "tab_firstorder.tex")), sep = "\n")
