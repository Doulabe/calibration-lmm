# ==============================================================================
# 35e_protected_profiles_figure.R
# Supplementary figure: profile of the exact-distribution p-value p(delta0) of the observed
# statistic over delta0, with its Monte Carlo band (p_hat +/- 2 se), for four decisive cases of
# the GEUVADIS application (33g and 27d results): ZNF266 and the PEX6-region lead in the Yoruba
# sample (GEMMA's Wald test), the second step of the likelihood-ratio path of this study in the
# PEX6 region (Yoruba), and TRA2A-AS (GEMMA's score, Yoruba). The coarse-grid maximum, the
# per-test level and the nominal p-value are drawn for reference.
# each panel carries the name of its case.
# Usage (from the repository root): [ART=<article folder>] Rscript calibration/35e_protected_profiles_figure.R
# Output: <ART>/figures/fig_protected_profiles.pdf
# ==============================================================================
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
ART <- Sys.getenv("ART", "calibration/output"); FIG <- file.path(ART, "figures", "fig_protected_profiles.pdf")
RES <- "calibration/results/protected_check"
cases <- list(
  list(f = "ZNF266_YRI_wald.rds", t = "ZNF266, Wald (GEMMA)"),                             # ZNF266, Yoruba, Wald test of GEMMA
  list(f = "PEX6_YRI_wald.rds", t = "PEX6 region, Wald (GEMMA)"),                               # PEX6 region, Yoruba, Wald test of GEMMA
  list(f = "PEX6_YRI_lrpath_step2.rds", t = "PEX6 region, LR path, step 2"),                       # PEX6 region, Yoruba, likelihood ratio, step 2
  list(f = "TRA2A_YRI_score.rds", t = "TRA2A-AS, score (GEMMA)"))                             # TRA2A-AS, Yoruba, score test of GEMMA
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; BAND <- adjustcolor("#2a78d6", 0.18)
X0 <- -3.7                                                     # abscissa of delta0 = 0
xpos <- function(d) ifelse(d == 0, X0, log10(d))
pdf(sub("\\.pdf$", "_raw.pdf", FIG), width = 7.2, height = 5.6, pointsize = 9, encoding = "WinAnsi")
par(mfrow = c(2, 2), mar = c(2.6, 3.2, 1.6, 0.6), oma = c(1.6, 1.4, 2.2, 0), mgp = c(2, 0.45, 0), tcl = -0.25,
    col.axis = INK2, col.lab = INK, fg = INK2)
for (k in seq_along(cases)) {
  r <- readRDS(file.path(RES, cases[[k]]$f)); inf <- r$info
  fb <- file.path(sub("protected_check$", "protected_check_bigB", RES), cases[[k]]$f)
  p4 <- r$profile[order(r$profile$delta0), ]                        # 4,000 / 10,000 draws (33g, 27d)
  p <- if (file.exists(fb)) { rb <- readRDS(fb); rb$profile[order(rb$profile$delta0), ] } else p4
  lo <- p$log10p + log10(pmax(1 - 2 * p$rel_se, 0.01)); up <- p$log10up
  shown <- p$delta0 == 0 | p$delta0 >= 1e-3
  yl <- range(c(p$log10p[p$delta0 >= 0.01], inf$level, inf$nominal, inf$stored_protected, up, p4$log10p[p4$delta0 >= 0.01]))
  yl <- c(floor(yl[1] - 0.3), ceiling(yl[2] + 0.3)); yl[1] <- max(yl[1], yl[2] - 12)
  plot(NA, xlim = c(X0 - 0.2, 10.2), ylim = yl, xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "l")
  abline(h = seq(ceiling(yl[1]), floor(yl[2]), 1), col = GRID, lwd = 0.6)
  abline(v = -3.35, col = GRID, lwd = 3)                       # break between 0 and 10^-3
  axis(1, at = c(X0, seq(-2, 10, 2)), labels = c("0", parse(text = paste0("10^", seq(-2, 10, 2)))), col = INK2, lwd = 0.6, cex.axis = 0.9)
  axis(2, at = seq(ceiling(yl[1]), floor(yl[2]), 1), las = 1, col = INK2, lwd = 0.6, cex.axis = 0.9)
  q <- p[shown & p$delta0 > 0, ]; ql <- lo[shown & p$delta0 > 0]; qu <- up[shown & p$delta0 > 0]
  if (file.exists(fb)) { g <- p4[p4$delta0 > 0 & p4$stage == 1, ]                # profile with 4,000 draws, in grey
    lines(log10(g$delta0), pmax(g$log10p, yl[1]), col = "#9a9893", lwd = 0.9) }
  polygon(c(log10(q$delta0), rev(log10(q$delta0))), c(pmax(ql, yl[1]), rev(pmin(qu, yl[2]))), col = BAND, border = NA)
  s1 <- q$stage == 1
  lines(log10(q$delta0[s1]), q$log10p[s1], col = BLUE, lwd = 1.2)
  points(log10(q$delta0[s1]), q$log10p[s1], pch = 16, cex = 0.45, col = BLUE)
  points(log10(q$delta0[!s1]), q$log10p[!s1], pch = 21, cex = 0.7, bg = "white", col = BLUE, lwd = 0.8)
  z <- p[p$delta0 == 0, ][1, ]; points(X0, z$log10p, pch = 16, cex = 0.6, col = BLUE)
  segments(X0, max(lo[p$delta0 == 0][1], yl[1]), X0, min(up[p$delta0 == 0][1], yl[2]), col = BLUE, lwd = 0.8)
  abline(h = inf$level, col = INK, lwd = 1, lty = 2)
  mtext(cases[[k]]$t, side = 3, line = 0.3, adj = 0, cex = 0.85, font = 2, col = INK)
  abline(h = inf$nominal, col = INK2, lwd = 1, lty = 3)
  points(xpos(inf$stored_argmax), inf$stored_protected, pch = 4, cex = 1.1, lwd = 1.6, col = ORANGE)
}
mtext(expression("true ratio " * delta[0] * " (log scale; " * delta[0] * " = 0 at the left)"), side = 1, line = 0.4, outer = TRUE, cex = 0.85, col = INK)
mtext(expression(log[10] * " p" * (delta[0]) * " of the observed statistic"), side = 2, line = 0.2, outer = TRUE, cex = 0.85, col = INK, las = 0)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", inset = c(0, 0.002), horiz = TRUE, bty = "n", cex = 0.85, text.col = INK, x.intersp = 0.6, seg.len = 1.6,
       legend = c("refined estimate", "near the maximum", expression(hat(p) %+-% 2 * " s.e."), "4,000 draws (Wald)", "coarse maximum", "level", "nominal"),
       col = c(BLUE, BLUE, adjustcolor(BLUE, 0.3), "#9a9893", ORANGE, INK, INK2), lty = c(1, NA, NA, 1, NA, 2, 3), lwd = c(1.2, NA, NA, 0.9, 1.6, 1, 1),
       pch = c(16, 21, 15, NA, 4, NA, NA), pt.bg = c(NA, "white", NA, NA, NA, NA, NA), pt.cex = c(0.6, 0.9, 1.8, NA, 1.1, NA, NA))
## the band colour in the legend
dev.off()
system2("gs", c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true",
                paste0("-sOutputFile=", FIG), sub("\\.pdf$", "_raw.pdf", FIG)))
file.remove(sub("\\.pdf$", "_raw.pdf", FIG))
cat("wrote", FIG, "\n")
