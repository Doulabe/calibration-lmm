# ==============================================================================
# 23_make_edgeworth_tables.R
# Table of the Edgeworth check (22_edgeworth_check.R): size ratios of the Wald and LR
# tests with and without the Rothenberg / Zhou-Stephens correction, by class of nu_delta.
# Also cross-checks wald_t and lrt_exactref against wald_h1 and lrt_ml of
# 06_treatments_pilot.R (same draws and seeds).
# Usage (from the repository root): Rscript calibration/23_make_edgeworth_tables.R
# Output: note/tab_edgeworth.tex, results/edgeworth_summary.txt
# ==============================================================================
E <- read.csv("calibration/results/edgeworth_check.csv")
T6 <- read.csv("calibration/results/treatments_pilot.csv")
out <- character(0); say <- function(...) out <<- c(out, sprintf(...))
chk <- merge(E[E$method %in% c("wald_t", "lrt_exactref"), c("config", "method", "alpha", "size_ratio")],
             transform(T6[T6$method %in% c("wald_h1", "lrt_ml"), c("config", "method", "alpha", "size_ratio")],
                       method = ifelse(method == "wald_h1", "wald_t", "lrt_exactref")),
             by = c("config", "method", "alpha"), suffixes = c("", "_06"))
say("cross-check with 06 (same draws): max |log ratio| = %.2g over %d rows",
    max(abs(log(chk$size_ratio / chk$size_ratio_06))), nrow(chk))

cfg <- unique(E[, c("config", "nu_delta")])
brk <- c(0, 20, 100, Inf); lbl <- c("$\\nu_\\delta<20$", "$20\\le\\nu_\\delta<100$", "$\\nu_\\delta\\ge100$")
E$cls <- cut(E$nu_delta, brk, labels = FALSE, right = FALSE)
meth <- c(wald_t = "Wald, Student", wald_edgeworth = "Wald, Edgeworth",
          lrt_chi2 = "RV, $\\chi^2_1$", lrt_edgeworth = "RV, Edgeworth")
fmt <- function(v) { f1 <- function(x) if (x >= 100) sprintf("%.0f", x) else if (x >= 0.1) sprintf("%.2g", x) else
  { e <- floor(log10(x)); m <- round(x / 10^e); if (m == 10) { m <- 1; e <- e + 1 }
  if (m == 1) sprintf("10^{%d}", e) else sprintf("%d\\cdot10^{%d}", m, e) }
  sprintf("$%s$ {\\scriptsize $[%s, %s]$}", f1(median(v)), f1(min(v)), f1(max(v))) }
rows <- character(0)
for (k in seq_along(lbl)) {
  d <- E[E$cls == k, ]; if (!nrow(d)) next
  ncf <- length(unique(d$config))
  rows <- c(rows, sprintf("\\multicolumn{4}{l}{%s, %d configurations} \\\\", lbl[k], ncf))
  for (m in names(meth)) {
    cells <- vapply(c(1e-3, 1e-5, 5e-8), function(a) fmt(d$size_ratio[d$method == m & abs(d$alpha - a) < a * 1e-6]), "")
    rows <- c(rows, sprintf("\\quad %s & %s \\\\", meth[m], paste(cells, collapse = " & ")))
  }
  if (k < length(lbl)) rows <- c(rows, "\\midrule")
}
writeLines(c("\\begin{tabular}{lccc}", "\\toprule",
  "Test & $\\alpha=10^{-3}$ & $\\alpha=10^{-5}$ & $\\alpha=5\\cdot10^{-8}$ \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}"), "calibration/note/tab_edgeworth.tex")

for (a in c(1e-3, 1e-5, 5e-8)) for (m in c("oracle_t", "wald_t", "wald_chi2", "wald_edgeworth", "lrt_chi2",
                                          "lrt_exactref", "lrt_edgeworth", "score_t", "score_edgeworth")) {
  d <- E[E$method == m & abs(E$alpha - a) < a * 1e-6, ]
  say("alpha %-6g %-16s size ratio median %.3g [%.3g, %.3g]; median power %.3f", a, m, median(d$size_ratio),
      min(d$size_ratio), max(d$size_ratio), median(d$power))
}
say("configs by nu_delta class: %s", paste(table(cut(cfg$nu_delta, brk, right = FALSE)), collapse = " / "))
say("nu_delta range: %.1f to %.1f", min(cfg$nu_delta), max(cfg$nu_delta))
## where does the corrected Wald cross from anti-conservative to conservative?
w <- E[E$method == "wald_edgeworth" & abs(E$alpha - 5e-8) < 1e-12, c("config", "nu_delta", "size_ratio")]
say("Wald Edgeworth at 5e-8 by config (nu_delta: ratio): %s",
    paste(sprintf("%.0f: %.2g", w$nu_delta[order(w$nu_delta)], w$size_ratio[order(w$nu_delta)]), collapse = "; "))
writeLines(out, "calibration/results/edgeworth_summary.txt"); cat(out, sep = "\n")
