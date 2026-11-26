# ==============================================================================
# 21_make_nu_table.R
# Table of the effective df nu_j against n at fixed cis windows (20_nu_vs_n.R).
# Pools the two windows (RPS26, PEX6 region; 30 candidates each).
# Usage (from the repository root): Rscript calibration/21_make_nu_table.R
# Output: note/tab_nu_n.tex
# ==============================================================================
R <- read.csv("calibration/results/nu_vs_n.csv")
cell <- function(d) sprintf("$%.0f$ {\\scriptsize [%.0f, %.0f]} & $%.2f$ & $%s$",
  median(d$nu_j), quantile(d$nu_j, .25), quantile(d$nu_j, .75), mean(d$nu_j < 100),
  trimws(formatC(median(d$tail_ratio), format = "fg", digits = 2)))
rows <- vapply(sort(unique(R$n)), function(n) sprintf("$%s$ & %s & %s \\\\", if (n >= 1000) sprintf("%d\\,%03d", n %/% 1000, n %% 1000) else n,
  cell(R[R$n == n & R$h2 == 0.10, ]), cell(R[R$n == n & R$h2 == 0.25, ])), "")
writeLines(c("\\begin{tabular}{rcccccc}", "\\toprule",
  "& \\multicolumn{3}{c}{$h^2_{\\rm fond}=0.10$} & \\multicolumn{3}{c}{$h^2_{\\rm fond}=0.25$} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "$n$ & $\\nu_j$ médian [quartiles] & $P(\\nu_j<100)$ & queue & $\\nu_j$ médian [quartiles] & $P(\\nu_j<100)$ & queue \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "calibration/note/tab_nu_n.tex")
cat(readLines("calibration/note/tab_nu_n.tex"), sep = "\n")
