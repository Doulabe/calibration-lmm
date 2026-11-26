# ==============================================================================
# 07_make_extension_tables.R
# LaTeX tables for the pilot extension (06_treatments_pilot.R): size ratio and
# power by treatment of delta, medians and ranges over configurations.
# Usage (from the repository root): Rscript calibration/07_make_extension_tables.R
# Output: calibration/note/tab_treat_*.tex, tab_ebic.tex
# ==============================================================================
f_final <- "calibration/results/treatments_pilot_final.csv"
d <- read.csv(if (file.exists(f_final)) f_final else "calibration/results/treatments_pilot.csv")
N <- "calibration/note"
fm <- function(v) {
  if (is.na(v)) return("--")
  if (v == 0) return("$0$")
  if (v >= 1000 || v < 0.01) {
    e <- floor(log10(v)); m <- v / 10^e
    return(sprintf("$%.1f\\times10^{%d}$", m, e))
  }
  if (v >= 10) return(sprintf("$%.0f$", v))
  sprintf("$%.2f$", v)
}
lab <- c(oracle_t = "Oracle $t(\\delta_0)$, $\\delta$ connu",
         wald_h1 = "Wald, $\\hat\\delta$ sous $H_1$ (GEMMA)",
         score_h0 = "Score, $\\hat\\delta_0$ sous $H_0$ (EMMAX, MLMM)",
         lrt_ml = "Rapport de vraisemblance ML",
         signflip = "Retournement de signes",
         bb_wald = "Berger--Boos (Wald exact)",
         oracle_bf = "Oracle BF($\\delta_0$), $\\delta$ connu",
         bf_h1 = "BF en $\\hat\\delta$ sous $H_1$",
         bf_h0 = "BF en $\\hat\\delta_0$ (bras REML)",
         mbf = "BF marginal (principal)")
meth <- names(lab)
AL <- c(1e-3, 1e-5, 5e-8)
summ <- function(sub, col) {
  vapply(meth, function(mm) vapply(AL, function(a) {
    x <- sub[[col]][sub$method == mm & abs(sub$alpha - a) < 1e-12]
    if (!length(x)) return(c(NA, NA, NA))
    c(median(x), min(x), max(x))
  }, numeric(3)), matrix(0, 3, 3))
}
cellm <- function(v) sprintf("%s {\\scriptsize $[%s,\\ %s]$}", fm(v[1]),
                             gsub("\\$", "", fm(v[2])), gsub("\\$", "", fm(v[3])))
## benchmark of attainable power: Student with the effective df,
## 1/nubar = 1/(nu-1) + 1/nu_delta, at the same noncentrality lambda = t_{nu-1, alpha}
bench <- unique(d[!is.na(d$alpha) & d$method == "oracle_t",
                  c("config", "design", "n", "nu_delta", "alpha", "lambda")])
nub <- 1 / (1 / (bench$n - 2) + 1 / bench$nu_delta)
cb  <- qt(1 - bench$alpha / 2, nub)
bench$power <- pt(-cb, nub, ncp = bench$lambda) + pt(cb, nub, ncp = bench$lambda, lower.tail = FALSE)
bench$method <- "bench_tnu"
make_tab <- function(col, file, design) {
  sub <- d[d$design == design, ]
  S <- summ(sub, col)                      # 3 stats x 3 alpha x methods
  body <- vapply(seq_along(meth), function(i) {
    cells <- vapply(1:3, function(k) cellm(S[, k, i]), "")
    paste(lab[i], "&", paste(cells, collapse = " & "), "\\\\")
  }, "")
  body <- append(body, "\\midrule", after = 6)
  if (col == "power") {
    bz <- bench[bench$design == design, ]
    cells <- vapply(AL, function(a) { x <- bz$power[abs(bz$alpha - a) < 1e-12]
      cellm(c(median(x), min(x), max(x))) }, "")
    body <- c(body, "\\midrule",
              paste("Référence $t_{\\bar\\nu_j}$ (ordre de grandeur) &",
                    paste(cells, collapse = " & "), "\\\\"))
  }
  writeLines(c("\\begin{tabular}{lccc}", "\\toprule",
               "Traitement de $\\delta$ & $\\alpha=10^{-3}$ & $\\alpha=10^{-5}$ & $\\alpha=5\\cdot10^{-8}$ \\\\",
               "\\midrule", body, "\\bottomrule", "\\end{tabular}"), file.path(N, file))
}
make_tab("size_ratio", "tab_treat_size_geu.tex", "GEUVADIS_cis")
make_tab("power", "tab_treat_power_geu.tex", "GEUVADIS_cis")
make_tab("size_ratio", "tab_treat_size_syn.tex", "AR1_rho0.95")
make_tab("power", "tab_treat_power_syn.tex", "AR1_rho0.95")

## eBIC decision: oracle LRT vs LRT with delta re-estimated
e <- d[d$method %in% c("ebic_oracle", "ebic_lrt_ml"), ]
rows <- c()
for (des in c("GEUVADIS_cis", "AR1_rho0.95")) for (h in c(0.05, 0.20)) {
  z <- e[e$design == des & abs(e$h2 - h) < 1e-9, ]
  g <- function(mm, lam0) {
    x <- z$power[z$method == mm & (if (lam0 == 0) z$lambda == 0 else z$lambda > lam0)]
    median(x) }
  lam_hi <- 5
  rows <- c(rows, sprintf("%s & %.2f & %s & %s & %s & %s \\\\",
    ifelse(des == "GEUVADIS_cis", "GEUVADIS", "AR(1)"), h,
    fm(g("ebic_oracle", 0)), fm(g("ebic_lrt_ml", 0)),
    fm(g("ebic_oracle", lam_hi)), fm(g("ebic_lrt_ml", lam_hi))))
}
writeLines(c("\\begin{tabular}{lrcccc}", "\\toprule",
  "& & \\multicolumn{2}{c}{Sous $H_0$} & \\multicolumn{2}{c}{Puissance, $\\lambda=t_{5\\cdot10^{-8}}$} \\\\",
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "Plan & $h^2$ & oracle & ML, $\\delta$ ré-estimé & oracle & ML, $\\delta$ ré-estimé \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_ebic.tex"))

## ratio of power to the oracle, by nu_delta (for the text)
pw <- d[d$method %in% meth & !is.na(d$alpha), c("config", "design", "h2", "nu_delta", "method",
                                                   "alpha", "size_ratio", "power")]
or <- pw[pw$method == "oracle_t", c("config", "alpha", "power")]; names(or)[3] <- "p_or"
pw <- merge(pw, or, by = c("config", "alpha"))
pw$rel_power <- pw$power / pw$p_or
write.csv(pw, "calibration/results/treatments_relpower.csv", row.names = FALSE)
## share of draws whose statistic peaks inside |u - split| < 25 (non-monotone in the effect)
if ("nonmonotone_share" %in% names(d)) {
  nm <- unique(d[d$method %in% c("bf_h0", "mbf", "score_h0", "signflip") & !is.na(d$nonmonotone_share),
                 c("config", "design", "method", "nonmonotone_share")])
  rows <- c()
  for (mm in c("score_h0", "signflip", "bf_h0", "mbf")) {
    g <- function(des) { x <- nm$nonmonotone_share[nm$method == mm & nm$design == des]
      sprintf("$%.0f$ {\\scriptsize $[%.0f,\\ %.0f]$}", 100 * median(x), 100 * min(x), 100 * max(x)) }
    rows <- c(rows, sprintf("%s & %s & %s \\\\", lab[mm], g("GEUVADIS_cis"), g("AR1_rho0.95")))
  }
  writeLines(c("\\begin{tabular}{lcc}", "\\toprule",
               "Statistique & GEUVADIS (\\%) & AR(1) (\\%) \\\\", "\\midrule", rows,
               "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_nonmono.tex"))
}
cat("extension tables written\n")
