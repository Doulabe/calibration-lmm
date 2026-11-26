# ==============================================================================
# 11_make_clr_tables.R
# Tables for the conditional likelihood-ratio test (10_clr_test.R):
#   tab_clr_size.tex   rejection rate under H0 at 0.05 and 0.01 (exactness check)
#   tab_clr_power.tex  paired power at lambda = t_{nu-1, alpha}
#   tab_clr_nu.tex     power of CLR and LRT by class of nu_j (all configurations)
# Usage (from the repository root): Rscript calibration/11_make_clr_tables.R
# ==============================================================================
d <- read.csv("calibration/results/clr_test.csv")
N <- "calibration/note"
fm <- function(v) {
  if (is.na(v)) return("--")
  if (v == 0) return("0")
  if (v < 0.01) { e <- floor(log10(v)); return(sprintf("%.1f\\times10^{%d}", v / 10^e, e)) }
  sprintf("%.2f", v)
}
cellm <- function(x) sprintf("$%s$ {\\scriptsize $[%s,\\ %s]$}", fm(median(x)), fm(min(x)), fm(max(x)))
lab <- c(oracle_t = "Oracle $t(\\delta_0)$, $\\delta$ connu",
         clr = "Rapport de vraisemblance conditionnel (exact)",
         lrt_ml = "Rapport de vraisemblance ML",
         mbf_hc = "BF marginal, $\\tau$ demi-Cauchy",
         signflip = "Retournement de signes sur le score",
         wald_h1 = "Wald, $\\hat\\delta$ sous $H_1$ (GEMMA)")
des_lab <- c(GEUVADIS_cis = "GEUVADIS, 16 configurations", AR1_rho0.95 = "AR(1), 8 configurations")

## size check
rows <- c()
for (des in names(des_lab)) {
  rows <- c(rows, sprintf("\\multicolumn{3}{l}{\\emph{%s}} \\\\", des_lab[des]))
  for (mm in c("oracle_t", "clr", "lrt_ml", "signflip", "wald_h1")) {
    z <- d[d$design == des & d$method == mm & d$lambda == 0, ]
    rows <- c(rows, sprintf("%s & %s & %s \\\\", lab[mm], cellm(z$reject[z$alpha == 0.05]),
                            cellm(z$reject[z$alpha == 0.01])))
  }
}
writeLines(c("\\begin{tabular}{lcc}", "\\toprule",
             "Traitement & $\\alpha=0.05$ & $\\alpha=0.01$ \\\\", "\\midrule", rows,
             "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_clr_size.tex"))

## paired power
AL <- c(1e-3, 1e-5, 5e-8)
rows <- c()
for (des in names(des_lab)) {
  rows <- c(rows, sprintf("\\multicolumn{4}{l}{\\emph{%s}} \\\\", des_lab[des]))
  for (mm in c("oracle_t", "clr", "lrt_ml", "mbf_hc", "signflip", "wald_h1")) {
    z <- d[d$design == des & d$method == mm & d$lambda > 0, ]
    cells <- vapply(AL, function(a) cellm(z$reject[abs(z$alpha - a) < 1e-12]), "")
    rows <- c(rows, sprintf("%s & %s \\\\", lab[mm], paste(cells, collapse = " & ")))
  }
}
writeLines(c("\\begin{tabular}{lccc}", "\\toprule",
             "Traitement & $\\alpha=10^{-3}$ & $\\alpha=10^{-5}$ & $\\alpha=5\\cdot10^{-8}$ \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_clr_power.tex"))

## CLR vs LRT power by class of nu_j (all 24 configurations)
p <- d[d$lambda > 0 & d$method %in% c("clr", "lrt_ml", "mbf_hc", "oracle_t"), ]
cls <- cut(p$nu_delta, c(0, 20, 50, 150, Inf), right = FALSE,
           labels = c("$\\nu_j<20$", "$20\\le\\nu_j<50$", "$50\\le\\nu_j<150$", "$\\nu_j\\ge150$"))
p$cls <- cls
rows <- c()
for (cl in levels(cls)) {
  z <- p[p$cls == cl, ]
  ncf <- length(unique(z$config))
  if (!ncf) next
  cells <- c()
  for (a in AL) for (mm in c("clr", "lrt_ml", "mbf_hc")) {
    x <- z$reject[z$method == mm & abs(z$alpha - a) < 1e-12]
    cells <- c(cells, sprintf("$%s$", fm(median(x))))
  }
  rows <- c(rows, sprintf("%s & %d & %s \\\\", cl, ncf, paste(cells, collapse = " & ")))
}
writeLines(c("\\begin{tabular}{lrccccccccc}", "\\toprule",
  "& & \\multicolumn{3}{c}{$\\alpha=10^{-3}$} & \\multicolumn{3}{c}{$\\alpha=10^{-5}$} & \\multicolumn{3}{c}{$\\alpha=5\\cdot10^{-8}$} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
  "Classe & configs & CLR & LRT & BF & CLR & LRT & BF & CLR & LRT & BF \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_clr_nu.tex"))
cat("CLR tables written\n")
