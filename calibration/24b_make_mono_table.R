# ==============================================================================
# 24b_make_mono_table.R
# Table of the score monotonicity check (24_score_monotonicity_check.R).
# Usage (from the repository root): Rscript calibration/24b_make_mono_table.R
# Output: note/tab_mono.tex
# ==============================================================================
x <- read.csv("calibration/results/score_monotonicity_check.csv")
lab <- function(cfg) {
  if (grepl("^geu_", cfg)) { p <- strsplit(sub("^geu_", "", cfg), "_h")[[1]]
    sprintf("GEUVADIS %s, $h^2=%s$", p[1], p[2]) } else {
    m <- regmatches(cfg, regexec("syn_n([0-9]+)_b([0-9]+)_h([0-9.]+)", cfg))[[1]]
    sprintf("AR(1), $n=%s$, %s blocs, $h^2=%s$", m[2], m[3], m[4]) } }
pct <- function(v) if (is.na(v)) "--" else sprintf("$%.1f$", 100 * v)
rows <- vapply(unique(x$config), function(cf) {
  d <- x[x$config == cf, ]; d <- d[order(-d$alpha), ]
  sprintf("%s & $%.1f$ & %s & %s & %s \\\\", lab(cf), d$nu_delta[1], pct(d$share_nonmonotone[1]),
          pct(d$dec_with_dhat0_jump[1]),
          paste(sprintf("$%+.1f$", 100 * d$rel_diff), collapse = " & ")) }, "")
writeLines(c("\\begin{tabular}{lrrrrrr}", "\\toprule",
  "& & \\multicolumn{2}{c}{Tirages (\\%)} & \\multicolumn{3}{c}{Écart de taille (\\%)} \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-7}",
  "Configuration & $\\nu_\\delta$ & non monotones & dont saut de $\\hat\\delta_0$ & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "calibration/note/tab_mono.tex")
cat(readLines("calibration/note/tab_mono.tex"), sep = "\n")
