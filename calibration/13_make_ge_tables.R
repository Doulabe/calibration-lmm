# ==============================================================================
# 13_make_ge_tables.R
# Table of the group-exchange test (12_group_exchange_test.R, reduced run) with the
# sign-flip CLR of 10_clr_test.R on the same eight GEUVADIS configurations.
# Usage (from the repository root): Rscript calibration/13_make_ge_tables.R
# Output: calibration/note/tab_ge.tex
# ==============================================================================
g <- read.csv("calibration/results/group_exchange_reduced.csv")
c <- read.csv("calibration/results/clr_test.csv")
c <- c[c$config %in% unique(g$config), ]
fm <- function(v) if (v < 0.01 && v > 0) { e <- floor(log10(v)); sprintf("%.1f\\times10^{%d}", v / 10^e, e) } else sprintf("%.3f", v)
cellm <- function(x) sprintf("$%s$ {\\scriptsize $[%s,\\ %s]$}", fm(median(x)), fm(min(x)), fm(max(x)))
row <- function(lab, d, m) {
  sz <- vapply(c(0.05, 0.01), function(a) cellm(d$reject[d$method == m & d$lambda == 0 & abs(d$alpha - a) < 1e-12]), "")
  pw <- vapply(c(1e-3, 1e-5, 5e-8), function(a) cellm(d$reject[d$method == m & d$lambda > 0 & abs(d$alpha - a) < 1e-12]), "")
  sprintf("%s & %s & %s \\\\", lab, paste(sz, collapse = " & "), paste(pw, collapse = " & "))
}
rows <- c(row("Oracle $t(\\delta_0)$, $\\delta$ connu", g, "oracle_t"),
          row("Échange, groupes de rapport 2", g, "ge_r2"),
          row("Échange, groupes de rapport 4", g, "ge_r4"),
          row("Rapport de vraisemblance ML", g, "lrt_ml"),
          row("Retournement de signes, CLR$^{*}$", c, "clr"))
writeLines(c("\\begin{tabular}{lccccc}", "\\toprule",
  "& \\multicolumn{2}{c}{Taille} & \\multicolumn{3}{c}{Puissance} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
  "Traitement & $0.05$ & $0.01$ & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "calibration/note/tab_ge.tex")
cat("tab_ge written\n")
