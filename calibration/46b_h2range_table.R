# ==============================================================================
# 46b_h2range_table.R
# Supplementary table of 46_calibration_h2_range.R: exact size relative to the nominal level at
# h2 = 0 and h2 = 0.5, as medians over the configurations of each class of nu_delta, at
# alpha = 5e-8 (1e-3 and the test with delta known are printed for the text), for the tests of Table 2 of the article; and paired comparisons of each
# candidate with itself at the design heritabilities of 36 (h2 = 0.05 for h2 = 0, 0.20 for 0.5).
# Usage (from the repository root): ART=<article folder> Rscript calibration/46b_h2range_table.R
# Output: calibration/output/tables/tabS_calib_h2range.tex
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/output")
RES <- "calibration/results"
rd <- function(d) do.call(rbind, lapply(list.files(file.path(RES, d), "rds$", full.names = TRUE), readRDS))
H <- rd("calib_h2range")
C0 <- rd("calib_enlarged"); E0 <- rd("calib_enlarged_score"); C0 <- rbind(C0, E0[E0$method == "emmax_t", names(C0)])
X <- read.csv(file.path(RES, "calib_h2range_extra.csv"))
fm <- function(v) { if (!is.finite(v)) return("--")              # as 37, ratios below 1e-6 censored
  if (v < 1e-6) return("$<10^{-6}$")
  if (v < 0.01 || v >= 1e4) { e <- floor(log10(v)); m <- v / 10^e; if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    return(if (abs(round(m, 1) - 1) < 1e-9) sprintf("$10^{%d}$", e) else sprintf("$%.1f\\times10^{%d}$", m, e)) }
  if (v >= 1000) return(sprintf("$%s$", formatC(round(v), format = "d", big.mark = "{,}")))
  if (v >= 10) sprintf("$%.0f$", v) else sprintf("$%.2f$", v) }
rows <- list(c("wald_t", "Wald, $\\widehat\\delta$ with the candidate"),             # the test with delta known is checked in the text
             c("emmax_t", "Null $\\widehat\\delta_0$, $t(\\widehat\\delta_0)$"), c("lrt_chi2", "Likelihood ratio, $\\chi^2_1$"),
             c("kr", "Kenward--Roger"), c("satterthwaite", "Satterthwaite"), c("wald_edgeworth", "Wald, Edgeworth"))
H$cls <- cut(H$nu_delta, c(0, 20, 100, Inf), right = FALSE)
cells <- expand.grid(cls = levels(H$cls), h2 = c(0, 0.5), stringsAsFactors = FALSE)[, c("h2", "cls")]
ncf <- vapply(seq_len(nrow(cells)), function(k) length(unique(H$config[H$h2 == cells$h2[k] & H$cls == cells$cls[k]])), 0L)
cells <- cells[ncf > 0, ]; ncf <- ncf[ncf > 0]                 # no configuration has nu_delta < 20 at h2 = 0.5
med <- function(m, a, k) { v <- H$size_ratio[H$method == m & abs(H$alpha - a) < a * 1e-6 & H$h2 == cells$h2[k] & H$cls == cells$cls[k]]
  if (!length(v)) "--" else fm(median(v)) }
blk <- function(a, lab) c(if (!is.null(lab)) sprintf("\\multicolumn{%d}{l}{\\emph{%s}} \\\\", nrow(cells) + 1L, lab),
  vapply(rows, function(rr) sprintf("%s & %s \\\\", rr[2], paste(vapply(seq_len(nrow(cells)), function(k) med(rr[1], a, k), ""), collapse = " & ")), ""))
cl <- c("$<20$", "$20$--$100$", "$\\ge100$")
n0 <- sum(cells$h2 == 0); n5 <- sum(cells$h2 == 0.5)
writeLines(c(sprintf("\\begin{tabular}{l%s}", strrep("c", nrow(cells))), "\\toprule",
             sprintf("& \\multicolumn{%d}{c}{$h^2=0$} & \\multicolumn{%d}{c}{$h^2=0.5$} \\\\", n0, n5),
             sprintf("\\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}", n0 + 1, n0 + 2, n0 + n5 + 1),
             sprintf("$\\nu_\\delta$ & %s \\\\", paste(cl[match(cells$cls, levels(H$cls))], collapse = " & ")),
             sprintf("Configurations & %s \\\\", paste(ncf, collapse = " & ")), "\\midrule",
             blk(5e-8, NULL),                                                    # 1e-3 given in the text
             "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabS_calib_h2range.tex"))
cat("wrote table; configurations per cell:", ncf, "\n")
## numbers for the text
options(width = 200)
print(aggregate(size_ratio ~ method + alpha + cls + h2, H[H$alpha %in% c(1e-3, 5e-8) & H$method %in% c("oracle_t", sapply(rows, `[`, 1)), ],
                function(v) signif(c(med = median(v), min = min(v), max = max(v)), 3)))
cat("\nboundary share P(delta_hat = 0) by h2 (median, range):\n"); print(aggregate(p_boundary ~ h2, X, function(v) round(c(med = median(v), min = min(v), max = max(v)), 3)))
cat("nu_delta range by h2:\n"); print(aggregate(nu_delta ~ h2, unique(H[, c("config", "h2", "nu_delta")]), range))
## paired comparison with the design heritabilities of 36: same candidate, h2 0.05 -> 0 and 0.20 -> 0.5
key <- function(d) sub("_h[0-9.]+$", "", d$config)
P <- merge(transform(H, k = key(H))[, c("k", "h2", "method", "alpha", "size_ratio", "nu_delta")],
           transform(C0, k = key(C0))[, c("k", "h2", "method", "alpha", "size_ratio", "nu_delta")],
           by = c("k", "method", "alpha"), suffixes = c("", ".d"))
P <- P[(P$h2 == 0 & abs(P$h2.d - 0.05) < 1e-9) | (P$h2 == 0.5 & abs(P$h2.d - 0.20) < 1e-9), ]
cat("\npaired: median ratio of nu_delta (new / design) by h2:\n"); print(aggregate(I(nu_delta / nu_delta.d) ~ h2, unique(P[, c("k", "h2", "nu_delta", "nu_delta.d")]), median))
cat("rank correlation of log10 size ratio with log nu_delta at 5e-8, by h2 and method:\n")
for (h in c(0, 0.5)) for (m in c("wald_t", "emmax_t", "lrt_chi2", "kr")) { d <- H[H$h2 == h & H$method == m & H$alpha == 5e-8, ]
  cat(sprintf("  h2 = %.1f %-10s %.2f\n", h, m, cor(log10(pmax(d$size_ratio, 1e-6)), log(d$nu_delta), method = "spearman"))) }
