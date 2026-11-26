# ==============================================================================
# 35f_merge_supp_tables.R
# Merges supplementary tables of the article to shorten the supplement, from the tables already
# written in the article's tables/ folder by 35_make_article_tables.R (parts calib and appendix):
#   tabA_calib_power.tex    size (tab_calibration) and power (tab_power) side by side, first design
#   tabA_pilot.tex          published GEUVADIS leads (tabA_geuvadis) and AR(1) pilot (tabA_synthetic)
#   tabA_clr_sizepower.tex  size (tabA_clr_size) and paired power (tabA_clr_power) of the exact tests
#   tabA_prior_u.tex        log marginal Bayes factor along u (tabA_prior_diag), six values of u
# The numbers are copied, not recomputed. Run after 35 (calib, appendix) and 35b.
# Usage (from the repository root): ART=<article folder> Rscript calibration/35f_merge_supp_tables.R
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/output")
TAB <- file.path(ART, "tables")
rd <- function(f) readLines(file.path(TAB, f))
wr <- function(lines, f) { writeLines(lines, file.path(TAB, f)); cat("wrote", f, "\n") }
is_row <- function(x) grepl("&", x, fixed = TRUE) & grepl("\\\\\\\\\\s*$", x)       # data or header row
is_block <- function(x) grepl("^\\\\multicolumn", x)                                   # design title row
cells <- function(x) trimws(strsplit(sub("\\\\\\\\\\s*$", "", x), "&", fixed = TRUE)[[1]])
row <- function(v) paste(paste(v, collapse = " & "), "\\\\")
block <- function(x, k) sub("^\\\\multicolumn\\{[0-9]+\\}", sprintf("\\\\multicolumn{%d}", k), x)
## data rows of a table, grouped by design block, as a list of label -> cells
by_block <- function(x) {
  out <- list(); cur <- "none"
  for (l in x) {
    if (is_block(l)) { cur <- l; out[[cur]] <- list(); next }
    if (is_row(l) && cur != "none") { v <- cells(l); out[[cur]][[v[1]]] <- v[-1] }
  }
  out
}

## ---- size and power of the first design, side by side ------------------------------------------
S <- rd("tab_calibration.tex"); P <- rd("tab_power.tex")
PB <- by_block(P); cur <- NULL; out <- c()
for (l in S) {
  if (grepl("^\\\\begin\\{tabular\\}", l)) { out <- c(out, "\\begin{tabular}{p{6.2cm}cccccc}"); next }
  if (grepl("^Treatment of", l)) {
    out <- c(out, "& \\multicolumn{3}{c}{Size relative to the nominal level} & \\multicolumn{3}{c}{Power at $\\lambda=t_{\\nu-1,1-\\alpha/2}$} \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
             row(c("Treatment of $\\delta$", rep(c("$\\alpha=10^{-3}$", "$\\alpha=10^{-5}$", "$\\alpha=5\\cdot10^{-8}$"), 2))))
    next }
  if (is_block(l)) { cur <- l; out <- c(out, block(l, 7)); next }
  if (is_row(l)) { v <- cells(l); p <- PB[[cur]][[v[1]]]; stopifnot(length(p) == 3)
    out <- c(out, row(c(v, p))); next }
  out <- c(out, l)
}
wr(out, "tabA_calib_power.tex")

## ---- pilot ratios: published GEUVADIS leads and AR(1) blocks -----------------------------------
G <- rd("tabA_geuvadis.tex"); A <- rd("tabA_synthetic.tex")
G <- G[is_row(G)][-1]; A <- A[is_row(A)][-1]                    # data rows (header dropped)
gr <- vapply(G, function(l) { v <- cells(l); row(c(v, "--")) }, "")            # Gene, h2, nu, W.05 ... S5e-8, unlinked
ar <- vapply(A, function(l) { v <- cells(l)                                     # n, blocks, h2, nu, W1e-3 ... unlinked
  row(c(sprintf("$n=%s$, %s blocks", v[1], v[2]), v[3], v[4], "--", v[5:10])) }, "")
wr(c("\\begin{tabular}{lrrrrrrrrr}", "\\toprule",
     row(c("Configuration", "$h^2$", "$\\nu_j$", "Wald $0.05$", "Wald $10^{-3}$", "Wald $5\\cdot10^{-8}$",
           "Satt.\\ $5\\cdot10^{-8}$", "Score $10^{-3}$", "Score $5\\cdot10^{-8}$", "Unlinked $5\\cdot10^{-8}$")),
     "\\midrule", "\\multicolumn{10}{l}{\\emph{Published lead of the eight GEUVADIS windows, $n=358$}} \\\\", unname(gr),
     "\\midrule", "\\multicolumn{10}{l}{\\emph{AR(1) blocks with $\\rho=0.95$, candidate at the centre of the first block}} \\\\",
     unname(ar), "\\bottomrule", "\\end{tabular}"), "tabA_pilot.tex")

## ---- exact tests: size and paired power side by side --------------------------------------------
Z <- by_block(rd("tabA_clr_size.tex")); W <- rd("tabA_clr_power.tex"); out <- c(); cur <- NULL
for (l in W) {
  if (grepl("^\\\\begin\\{tabular\\}", l)) { out <- c(out, "\\begin{tabular}{lccccc}"); next }
  if (grepl("^Treatment", l)) {
    out <- c(out, "& \\multicolumn{2}{c}{Size under $H_0$} & \\multicolumn{3}{c}{Paired power at $\\lambda=t_{\\nu-1,1-\\alpha/2}$} \\\\",
             "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}",
             row(c("Treatment", "$\\alpha=0.05$", "$\\alpha=0.01$", "$\\alpha=10^{-3}$", "$\\alpha=10^{-5}$", "$\\alpha=5\\cdot10^{-8}$")))
    next }
  if (is_block(l)) { cur <- sub("\\{4\\}", "{3}", l); out <- c(out, block(l, 6)); next }
  if (is_row(l)) { v <- cells(l); z <- Z[[cur]][[v[1]]]; if (is.null(z)) z <- c("--", "--")   # no size run for this treatment
    out <- c(out, row(c(v[1], z, v[-1]))); next }
  out <- c(out, l)
}
wr(out, "tabA_clr_sizepower.tex")

## ---- log marginal Bayes factor along u, six values ----------------------------------------------
D <- rd("tabA_prior_diag.tex"); h <- cells(D[grepl("^Standardized", D)]); keep <- c(1, which(h %in% c("0", "10", "20", "40", "120", "1000")))
stopifnot(length(keep) == 7)
out <- c("\\begin{tabular}{lrrrrrr}")
for (l in D[-1]) {
  if (is_block(l)) { out <- c(out, block(l, 7)); next }
  if (is_row(l)) { v <- cells(l)[keep]; if (!grepl("^Standardized", v[1])) v[-1] <- sprintf("$%s$", v[-1])   # minus signs
    out <- c(out, row(v)); next }
  out <- c(out, l)
}
wr(out, "tabA_prior_u.tex")
