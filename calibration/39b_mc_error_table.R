# ==============================================================================
# 39b_mc_error_table.R
# Supplementary table of the Monte Carlo error check (39_mc_error.R): size ratio with the
# stored 1,000 draws and with 10,000 fresh draws, and standard error of log10 of the size
# ratio for 1,000 draws (from the per-draw variance of the 10,000-draw run).
# Usage (from the repository root): ART=<article folder> Rscript calibration/39b_mc_error_table.R
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/output")
M <- read.csv("calibration/results/mc_error.csv")
fm <- function(v) { if (!is.finite(v)) return("--"); if (v < 0.01 || v >= 1e4) { e <- floor(log10(v)); m <- v / 10^e
    if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    return(if (abs(round(m, 1) - 1) < 1e-9) sprintf("$10^{%d}$", e) else sprintf("$%.1f\\times10^{%d}$", m, e)) }
  if (v >= 1000) return(sprintf("$%s$", formatC(round(v), format = "d", big.mark = "{,}")))
  if (v >= 10) sprintf("$%.0f$", v) else sprintf("$%.2f$", v) }
lab <- function(cf) { if (grepl("unlinked", cf)) return("AR(1), unrelated candidate")
  if (grepl("^syn", cf)) return(sub("^syn_n([0-9]+)_b([0-9]+)_h.*$", "AR(1), $n=\\1$, \\2 blocks", cf))
  sub("^geu_([^_]+)_.*$", "\\\\textit{\\1}", cf) }
tests <- c(wald_t = "Wald", emmax_t = "$t(\\widehat\\delta_0)$", lrt_chi2 = "LR", kr = "KR")
U <- unique(M[, c("config", "nu_delta")]); U <- U[order(U$nu_delta), ]
body <- unlist(lapply(seq_len(nrow(U)), function(i) { cf <- U$config[i]
  rows <- vapply(names(tests), function(te) { z <- M[M$config == cf & M$method == te, ]
    cells <- unlist(lapply(c(1e-5, 5e-8), function(a) { r <- z[abs(z$alpha - a) < a * 1e-6, ]
      c(fm(r$ratio_1k), fm(r$ratio_10k), sprintf("$%.2f$", r$se_log10_1k), sprintf("$%s$", formatC(round(r$ess), format = "d", big.mark = "{,}"))) }))
    paste(c(if (te == "wald_t") sprintf("%s, $\\nu_\\delta=%s$", lab(cf), formatC(signif(U$nu_delta[i], 3), format = "fg", big.mark = "{,}")) else "",
            tests[[te]], cells), collapse = " & ") }, "")
  c(if (i > 1) "\\addlinespace", paste0(rows, " \\\\")) }))
writeLines(c("\\begin{tabular}{llcccccccc}", "\\toprule",
             "& & \\multicolumn{4}{c}{$\\alpha=10^{-5}$} & \\multicolumn{4}{c}{$\\alpha=5\\cdot10^{-8}$} \\\\",
             "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
             "Configuration & Test & $B=10^3$ & $B=10^4$ & s.e. & eff. & $B=10^3$ & $B=10^4$ & s.e. & eff. \\\\", "\\midrule",
             body, "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabS_mc_error.tex"))
cat("wrote tabS_mc_error.tex\n")
