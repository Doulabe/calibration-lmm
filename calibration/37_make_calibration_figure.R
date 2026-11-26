# ==============================================================================
# 37_make_calibration_figure.R
# Figure 1 of the article: exact size of four tests relative to the
# nominal level as a function of nu_delta (36_calibration_enlarged.R), small multiples
# with a common axis, and the distribution of nu_j over all variants of the GEUVADIS
# windows on the same horizontal axis (results/nu_delta_per_snp_geuvadis.csv).
# Colours: slots 1-2 of the validated reference palette of the dataviz method
# (blue #2a78d6 for alpha = 1e-5, orange #eb6834 for alpha = 5e-8); text in neutral ink.
# Also writes the nu-stratified table of the main text (tab_corrections_main.tex) and the
# supplementary table of ranges (tabS_calib_enlarged.tex), with the plug-in rows of 50_protected_plugin_eval.R
# and the protected rows of 52_protected_design.R (50-51 when 52 has not run). The Wald and LR
# rows and points also come from the evaluation draws of 52. Each panel is titled.
# Usage (from the repository root):
#   ART=<article folder> Rscript calibration/37_make_calibration_figure.R
# ==============================================================================
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))  # the en dash of Kenward-Roger needs a UTF-8 locale
ART <- Sys.getenv("ART", "calibration/output")
RES <- "calibration/results"
fs <- list.files(file.path(RES, "calib_enlarged"), "rds$", full.names = TRUE)
C <- do.call(rbind, lapply(fs, readRDS))
## t(delta_hat_0) referred to t_{nu-1} (EMMAX, P3D, MLMM step), on the same draws (36b_score_reference.R)
fs2 <- list.files(file.path(RES, "calib_enlarged_score"), "rds$", full.names = TRUE)
E <- do.call(rbind, lapply(fs2, readRDS)); C <- rbind(C, E[E$method == "emmax_t", names(C)])
cat("configurations:", length(unique(C$config)), " with emmax_t:", length(unique(E$config)), "\n")
## Wald and LR rows from the evaluation draws of 52_protected_design.R (25,000 to 1,000,000
## independent draws per configuration, the draws of the protected rows) in place of the 1,000 draws of 36,
## in the figure and in the tables; the other treatments keep the 1,000 draws of 36
fs5 <- list.files(file.path(RES, "protected_design"), "rds$", full.names = TRUE)
E5 <- if (length(fs5)) do.call(rbind, lapply(fs5, function(f) { e <- readRDS(f)$eval; e[e$kind == "config", ] }))
if (!is.null(E5)) {
  N5 <- E5[E5$method %in% c("wald_t", "lrt_chi2"), ]
  C <- rbind(C[!paste(C$config, C$method) %in% paste(N5$config, N5$method), ], N5[, names(C)])
  cat("Wald and LR rows from 52:", nrow(N5), "rows,", length(unique(N5$config)), "configurations\n")
  if (length(unique(N5$config)) < length(unique(C$config))) warning("Wald and LR rows of 52 incomplete")
}
## plug-in and protected tests, thresholds from the exact distributions (49_protected_plugin_v2.R)
## and size and power at delta0 from 100,000 independent draws (50_protected_plugin_eval.R; 10,000 in 49
## when 50 has not run); rows of the tables only, not of the figure
fs3 <- list.files(file.path(RES, "protected_plugin_eval"), "rds$", full.names = TRUE)
if (!length(fs3)) fs3 <- list.files(file.path(RES, "protected_plugin_v2"), "rds$", full.names = TRUE)
if (length(fs3)) {
  V <- do.call(rbind, lapply(fs3, readRDS))
  ## protected thresholds recomputed with more common draws where 50 exceeded the level (51_protected_refine.R)
  fs4 <- list.files(file.path(RES, "protected_refine"), "rds$", full.names = TRUE)
  if (length(fs4)) {
    Rf <- do.call(rbind, lapply(fs4, readRDS))
    key <- function(z) paste(z$config, z$method, signif(z$alpha, 6))
    V <- rbind(V[!key(V) %in% key(Rf), names(C)], Rf[, names(C)])
    cat("protected rows refined:", nrow(Rf), "in", length(fs4), "configurations\n")
  }
  ## protected critical values computed once per design without delta0 and frozen
  ## (52_protected_design.R), evaluated at the true ratio of each configuration; they replace 50-51
  if (!is.null(E5)) {
    P5 <- E5[E5$method %in% c("wald_protected", "lr_protected"), ]
    V <- rbind(V[!V$method %in% c("wald_protected", "lr_protected"), names(C)], P5[, names(C)])
    cat("protected rows from design-level critical values:", nrow(P5), "rows,", length(unique(P5$config)), "configurations\n")
    if (length(unique(P5$config)) < length(unique(C$config))) warning("design-level protected rows incomplete")
  }
  C <- rbind(C, V[V$method %in% c("wald_plugin", "lr_plugin", "wald_protected", "lr_protected"), names(C)])
  cat("plug-in and protected tests:", length(fs3), "configurations from", dirname(fs3[1]), "\n")
  if (length(fs3) < length(unique(C$config))) warning("plug-in and protected rows incomplete")
}
## restricted likelihood ratio (chi2_1) and size-adjusted powers (56_restricted_lr_adjusted_power.R)
fs6 <- list.files(file.path(RES, "rlr_adjusted"), "rds$", full.names = TRUE)
ADJ <- NULL
if (length(fs6)) {
  R6 <- do.call(rbind, lapply(fs6, readRDS))
  C <- rbind(C, R6[R6$method == "rlr_chi2", names(C)])
  ADJ <- R6[grepl("_adj$", R6$method), ]
  ADJ$method <- c(wald_adj = "wald_t", score_adj = "emmax_t", lrt_adj = "lrt_chi2", rlr_adj = "rlr_chi2")[ADJ$method]
  cat("restricted LR and size-adjusted powers:", length(fs6), "configurations\n")
  if (length(fs6) < length(unique(C$config))) warning("restricted LR rows incomplete")
}
## modified signed likelihood root r* of Skovgaard (60_skovgaard_rstar.R), supplementary table only
fs7 <- list.files(file.path(RES, "skovgaard_rstar"), "rds$", full.names = TRUE)
if (length(fs7)) {
  R7 <- do.call(rbind, lapply(fs7, readRDS))
  C <- rbind(C, R7[R7$method == "rstar", names(C)])
  cat("Skovgaard r*:", length(fs7), "configurations\n")
  if (length(fs7) < length(unique(C$config))) warning("r* rows incomplete")
}
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; COL <- c("1e-05" = "#2a78d6", "5e-08" = "#eb6834")
## ratios below 1e-6 rest on a few draws (Table S4), so they are censored at 1e-6
## in the figure and in the tables (triangles, "<10^{-6}") instead of being plotted down to 1e-20.
FLOOR <- -6; XL <- c(0.5, 4.2); YT <- seq(-6, 4, 2)
## each panel carries the name of its test
panels <- c("wald_t", "emmax_t", "kr", "lrt_chi2", "rlr_chi2", "lr_protected")
PTITLE <- c("Wald, candidate estimate", "Null estimate", "Kenward\u2013Roger", "Likelihood ratio",
            "Restricted likelihood ratio", "Likelihood ratio, protected")
lx <- function(v) log10(pmin(v, 10^XL[2]))
FIG <- file.path(ART, "figures", "fig_calibration_nu.pdf")
pdf(sub("\\.pdf$", "_raw.pdf", FIG), width = 7.2, height = 8.3, pointsize = 9, encoding = "WinAnsi")
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7), 4, 2, byrow = TRUE), heights = c(1, 1, 1, 0.8))
par(oma = c(2.6, 3.6, 2.6, 0.4), mar = c(1.2, 1.6, 1.5, 0.6), mgp = c(2, 0.45, 0), tcl = -0.25, las = 1,
    col.axis = INK2, col.lab = INK, fg = INK2)
for (k in seq_along(panels)) {
  d <- C[C$method == panels[k] & C$alpha %in% c(1e-5, 5e-8), ]
  plot(NA, xlim = XL, ylim = c(FLOOR - 0.5, 5.2), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  abline(h = YT, col = GRID, lwd = 0.6); abline(v = log10(c(20, 100)), col = GRID, lwd = 0.8, lty = 2)
  rect(XL[1] - 1, -log10(3), XL[2] + 1, log10(3), col = adjustcolor("#9a9893", 0.18), border = NA)
  abline(h = 0, col = INK2, lwd = 0.9)
  mtext(PTITLE[k], side = 3, line = 0.25, adj = 0, cex = 0.85, font = 2, col = INK)
  axis(1, at = 1:4, labels = if (k >= 5) c("10", "100", "1,000", "10,000") else FALSE, col = INK2, lwd = 0.6)
  axis(2, at = YT, labels = if (k %in% c(1, 3, 5)) c("<-6", "-4", "-2", "0", "2", "4") else FALSE, col = INK2, lwd = 0.6)
  ## AR(1) designs and unlinked candidates drawn last, as larger outlined squares
  for (syn in c(FALSE, TRUE)) for (a in c("1e-05", "5e-08")) {
    z <- d[(format(d$alpha, scientific = TRUE) == a | abs(d$alpha - as.numeric(a)) < 1e-12) & grepl("^syn", d$config) == syn, ]
    y <- log10(pmax(z$size_ratio, 1e-300)); below <- y < FLOOR; y[below] <- FLOOR
    points(lx(z$nu_delta), y, pch = ifelse(below, 25, if (syn) 22 else 21), cex = if (syn) 1.55 else 1.35,
           bg = COL[a], col = if (syn) INK else "white", lwd = if (syn) 0.8 else 0.6) }
}
## distribution of nu_j over all variants of the windows
P <- read.csv(file.path(RES, "nu_delta_per_snp_geuvadis.csv"))
br <- seq(XL[1], XL[2], by = 0.1)
h1 <- hist(lx(P$nu_delta[abs(P$h2 - 0.05) < 1e-9]), breaks = c(-Inf, br[-1]), plot = FALSE)
h2 <- hist(lx(P$nu_delta[abs(P$h2 - 0.20) < 1e-9]), breaks = c(-Inf, br[-1]), plot = FALSE)
par(mar = c(1.2, 1.6, 1.0, 0.6))
ymax <- max(c(h1$counts / sum(h1$counts), h2$counts / sum(h2$counts)))
plot(NA, xlim = XL, ylim = c(0, ymax * 1.08), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
abline(v = log10(c(20, 100)), col = GRID, lwd = 0.8, lty = 2)
stp <- function(h, col, lwd) { x <- br; f <- h$counts / sum(h$counts); lines(rep(x, each = 2)[-c(1, 2 * length(x))], rep(f, each = 2), col = col, lwd = lwd) }
stp(h1, "#3d3c39", 1.6); stp(h2, "#9a9893", 1.6)
axis(1, at = 1:4, labels = c("10", "100", "1,000", "10,000"), col = INK2, lwd = 0.6)
axis(2, at = pretty(c(0, ymax)), labels = paste0(100 * pretty(c(0, ymax)), "%"), col = INK2, lwd = 0.6)
m1 <- which.max(h1$counts); m2 <- which.max(h2$counts)
text(br[m1] - 0.02, 1.04 * h1$counts[m1] / sum(h1$counts), expression(h^2 == "0.05"), pos = 2, cex = 0.85, col = "#3d3c39")
text(br[m2] + 0.12, 1.02 * h2$counts[m2] / sum(h2$counts), expression(h^2 == "0.20"), pos = 4, cex = 0.85, col = "#77756f")
mtext(expression("index " * nu[j] * " (" * nu[delta] * " at the true ratio), log scale"), side = 1, line = 1.2, outer = TRUE, cex = 0.85, col = INK)
mtext(expression(log[10] * " (estimated size / nominal level)"), side = 2, line = 2.2, outer = TRUE, cex = 0.85, col = INK, at = 0.61, las = 0)
mtext("share of variants", side = 2, line = 2.2, outer = TRUE, cex = 0.85, col = INK, at = 0.10, las = 0)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
legend("top", inset = c(0, 0.004), horiz = TRUE, bty = "n", pch = c(21, 21, 21, 22), pt.bg = c(COL, "#9a9893", "#9a9893"),
       col = c("white", "white", "white", INK), pt.lwd = c(0.6, 0.6, 0.6, 0.8), pt.cex = c(1.9, 1.9, 1.8, 1.9), cex = 1.1,
       legend = c(expression(alpha == 10^-5), expression(alpha == 5 %.% 10^-8), "GEUVADIS candidates", "AR(1) designs"),
       text.col = INK, x.intersp = 1.1, text.width = c(0.13, 0.17, 0.24, 0.16))
dev.off()
## embed the fonts (Ghostscript), as for the other figures of the article
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true",
                paste0("-sOutputFile=", FIG), sub("\\.pdf$", "_raw.pdf", FIG)))
file.remove(sub("\\.pdf$", "_raw.pdf", FIG))
cat("wrote figure\n")

## nu-stratified table for the main text: size and power, including the score and the reference
M <- C[C$alpha %in% c(1e-3, 1e-5, 5e-8), ]
M$cls <- cut(M$nu_delta, c(0, 20, 100, Inf), right = FALSE)
fm <- function(v) { if (!is.finite(v)) return("--"); if (v == 0) return("$0$")
  if (v < 1e-6) return("$<10^{-6}$")                            # censored, see the FLOOR comment
  if (v < 0.01 || v >= 1e4) { e <- floor(log10(v)); m <- v / 10^e; if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    return(if (abs(round(m, 1) - 1) < 1e-9) sprintf("$10^{%d}$", e) else sprintf("$%.1f\\times10^{%d}$", m, e)) }
  if (v >= 1000) return(sprintf("$%s$", formatC(round(v), format = "d", big.mark = "{,}")))
  if (v >= 10) sprintf("$%.0f$", v) else sprintf("$%.2f$", v) }
rows <- list(c("oracle_t", "Known $\\delta$ (reference)"), c("wald_t", "Wald, $\\widehat\\delta$ with the candidate"),
             c("emmax_t", "Null $\\widehat\\delta_0$, $t(\\widehat\\delta_0)$"), c("lrt_chi2", "Likelihood ratio, $\\chi^2_1$"),
             c("rlr_chi2", "Restricted likelihood ratio, $\\chi^2_1$"),
             c("kr", "Kenward--Roger"), c("satterthwaite", "Satterthwaite"), c("wald_edgeworth", "Wald, Edgeworth"),
             c("sep", ""), c("lr_protected", "Likelihood ratio, protected"), c("wald_protected", "Wald, protected"))
rows <- Filter(function(rr) rr[1] == "sep" || any(C$method == rr[1]), rows)
if (rows[[length(rows)]][1] == "sep") rows <- rows[-length(rows)]
lab <- c("$\\nu_\\delta<20$", "$20\\le\\nu_\\delta<100$", "$\\nu_\\delta\\ge100$")
body <- unlist(lapply(seq_along(levels(M$cls)), function(q) { d <- M[M$cls == levels(M$cls)[q], ]
  c(if (q > 1) "\\midrule", sprintf("\\multicolumn{10}{l}{%s, %d configurations} \\\\", lab[q], length(unique(d$config))),
    vapply(rows, function(rr) { if (rr[1] == "sep") return("\\addlinespace[3pt]"); z <- d[d$method == rr[1], ]
      sz <- vapply(c(1e-3, 1e-5, 5e-8), function(a) fm(median(z$size_ratio[abs(z$alpha - a) < a * 1e-6])), "")
      pw <- vapply(c(1e-3, 1e-5, 5e-8), function(a) fm(median(z$power[abs(z$alpha - a) < a * 1e-6])), "")
      za <- if (is.null(ADJ)) NULL else ADJ[ADJ$method == rr[1] & ADJ$config %in% d$config, ]
      pa <- if (rr[1] == "oracle_t") pw else if (is.null(za) || !nrow(za)) rep("--", 3) else
        vapply(c(1e-3, 1e-5, 5e-8), function(a) fm(median(za$power[abs(za$alpha - a) < a * 1e-6])), "")
      sprintf("\\quad %s & %s & %s & %s \\\\", rr[2], paste(sz, collapse = " & "), paste(pw, collapse = " & "), paste(pa, collapse = " & ")) }, "")) }))
writeLines(c("\\begin{tabular}{lccccccccc}", "\\toprule",
             "& \\multicolumn{3}{c}{Size relative to $\\alpha$} & \\multicolumn{3}{c}{Power} & \\multicolumn{3}{c}{Size-adjusted power} \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
             "Test & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ \\\\", "\\midrule",
             body, "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tab_corrections_main.tex"))
cat("wrote table\n")

## supplementary table: median and range of the size ratio by class of nu_delta, all treatments of the enlarged design
rows2 <- list(c("oracle_t", "Known $\\delta$, Student $t_{\\nu-1}$ (reference)"),
              c("wald_t", "Wald, $\\widehat\\delta$ with the candidate (GEMMA, MM4LMM)"),
              c("emmax_t", "Null $\\widehat\\delta_0$, $t(\\widehat\\delta_0)$ (EMMAX, P3D, MLMM step)"),
              c("score_t", "Null $\\widehat\\delta_0$, score $Q(\\widehat\\delta_0)$ (GEMMA score test)"),
              c("lrt_chi2", "Likelihood ratio, $\\delta$ under each model (GEMMA LRT)"),
              c("rlr_chi2", "Restricted likelihood ratio on a common residual space"),
              c("rstar", "Modified signed likelihood root $r^*$ (Skovgaard)"),
              c("kr", "Kenward--Roger (pbkrtest, MM4LMM, SAS)"),
              c("satterthwaite", "Satterthwaite degrees of freedom (lmerTest)"),
              c("wald_edgeworth", "Wald, Edgeworth threshold (GEMMA \\texttt{-crt})"),
              c("lrt_edgeworth", "Likelihood ratio, Edgeworth threshold"),
              c("wald_plugin", "Wald, plug-in at $\\widehat\\delta$"),
              c("lr_plugin", "Likelihood ratio, plug-in at $\\widehat\\delta$"),
              c("wald_protected", "Wald, protected"),
              c("lr_protected", "Likelihood ratio, protected"))
rows2 <- Filter(function(rr) any(C$method == rr[1]), rows2)
fr <- function(v) sprintf("%s {\\scriptsize [%s, %s]}", fm(median(v)), fm(min(v)), fm(max(v)))
body2 <- unlist(lapply(seq_along(levels(M$cls)), function(q) { d <- M[M$cls == levels(M$cls)[q], ]
  cf <- unique(d$config); ng <- sum(grepl("^geu_", cf))
  c(if (q > 1) "\\midrule", sprintf("\\multicolumn{4}{l}{\\emph{%s: %d configurations (%d GEUVADIS, %d synthetic)}} \\\\",
                                    lab[q], length(cf), ng, length(cf) - ng),
    vapply(rows2, function(rr) { z <- d[d$method == rr[1], ]
      sprintf("%s & %s \\\\", rr[2], paste(vapply(c(1e-3, 1e-5, 5e-8), function(a) fr(z$size_ratio[abs(z$alpha - a) < a * 1e-6]), ""), collapse = " & ")) }, "")) }))
writeLines(c("\\begin{tabular}{p{6.6cm}ccc}", "\\toprule",
             "Treatment of $\\delta$ & $\\alpha=10^{-3}$ & $\\alpha=10^{-5}$ & $\\alpha=5\\cdot10^{-8}$ \\\\", "\\midrule",
             body2, "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabS_calib_enlarged.tex"))
cat("wrote supplementary table\n")
## summaries for the text
S <- aggregate(size_ratio ~ cls + method + alpha, M, function(v) c(med = median(v), min = min(v), max = max(v)))
print(S[S$alpha == 5e-8, ], digits = 3)
cat("nu range:", range(unique(C[, c("config", "nu_delta")])$nu_delta), "\n")
print(table(cut(unique(C[, c("config", "nu_delta")])$nu_delta, c(0, 20, 100, Inf), right = FALSE)))
