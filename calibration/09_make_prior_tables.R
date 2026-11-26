# ==============================================================================
# 09_make_prior_tables.R
# Tables for the correction of the MBF non-monotonicity (08_mbf_prior_fix.R):
#   tab_prior_compare.tex  size ratio, power and non-monotone share of the three
#                          MBF versions (fixed tau2, half-Cauchy tau, hyper-g)
#   tab_prior_diag.tex     log MBF against the effect u, ERAP2, three versions
# Usage (from the repository root): Rscript calibration/09_make_prior_tables.R
# ==============================================================================
R <- "calibration/results"; N <- "calibration/note"
fin <- read.csv(file.path(R, "treatments_pilot_final.csv"))
fx  <- read.csv(file.path(R, "mbf_prior_fix.csv"))
old <- fin[fin$method == "mbf", ]; old$method <- "mbf_fixed"
keep <- c("config", "design", "h2", "method", "alpha", "size_ratio", "power", "nonmonotone_share")
d <- rbind(old[, keep], fx[, keep])
fm <- function(v) {
  if (is.na(v)) return("--")
  if (v == 0) return("0")
  if (v >= 1000 || v < 0.01) { e <- floor(log10(v)); return(sprintf("%.1f\\times10^{%d}", v / 10^e, e)) }
  if (v >= 10) return(sprintf("%.0f", v))
  sprintf("%.2f", v)
}
cellm <- function(x) sprintf("$%s$ {\\scriptsize $[%s,\\ %s]$}", fm(median(x)), fm(min(x)), fm(max(x)))
lab <- c(mbf_fixed = "$\\tau^2=0.04$ fixé",
         mbf_hc = "$\\tau\\sim$ demi-Cauchy$(0,\\,0.2)$",
         mbf_hg = "hyper-$g$, $a=3$")
AL <- c(1e-3, 1e-5, 5e-8)
rows <- c()
for (des in c("GEUVADIS_cis", "AR1_rho0.95")) {
  rows <- c(rows, sprintf("\\multicolumn{8}{l}{\\emph{%s}} \\\\",
                          ifelse(des == "GEUVADIS_cis", "GEUVADIS, 16 configurations",
                                 "AR(1), 8 configurations")))
  for (mm in names(lab)) {
    z <- d[d$design == des & d$method == mm, ]
    sz <- vapply(AL, function(a) cellm(z$size_ratio[abs(z$alpha - a) < 1e-12]), "")
    pw <- vapply(AL, function(a) cellm(z$power[abs(z$alpha - a) < 1e-12]), "")
    nm <- unique(z[, c("config", "nonmonotone_share")])$nonmonotone_share
    rows <- c(rows, sprintf("%s & %s & %s & %s \\\\", lab[mm],
      paste(sz, collapse = " & "), paste(pw, collapse = " & "), sprintf("$%.0f$", 100 * median(nm))))
  }
}
writeLines(c("\\begin{tabular}{lccccccc}", "\\toprule",
  "& \\multicolumn{3}{c}{Taille rapportée au nominal} & \\multicolumn{3}{c}{Puissance} & Non monotone \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "A priori sur l'effet & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ & (\\%) \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_prior_compare.tex"))

dg <- read.csv(file.path(R, "mbf_prior_diag.csv"))
us <- sort(unique(dg$u))
rows <- c()
for (h in c(0.05, 0.20)) {
  z <- dg[abs(dg$h2 - h) < 1e-9, ]; z <- z[order(z$u), ]
  rows <- c(rows, sprintf("\\multicolumn{%d}{l}{\\emph{$h^2=%.2f$}} \\\\", length(us) + 1, h))
  for (v in c("fixed", "hc", "hg"))
    rows <- c(rows, sprintf("\\quad %s & %s \\\\", c(fixed = "$\\tau^2$ fixé", hc = "demi-Cauchy",
                                                     hg = "hyper-$g$")[v],
                            paste(sprintf("%.1f", z[[v]]), collapse = " & ")))
}
writeLines(c(sprintf("\\begin{tabular}{l%s}", strrep("r", length(us))), "\\toprule",
  sprintf("Effet standardisé $u$ & %s \\\\", paste(sprintf("%.0f", us), collapse = " & ")),
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_prior_diag.tex"))
cat("prior tables written\n")
