# ==============================================================================
# 28_make_geuvadis_tables.R
# Tables of the GEUVADIS application (27_geuvadis_application.R).
#   note/tab_geu_paths.tex     number of variants selected per method, by gene and population
#   note/tab_geu_variants.tex  per selected variant: nu_j, nominal and exact log10 p-values
# Exact p-values below 1e-12 are reported as "< -12" (Monte Carlo resolution, B = 1000 draws).
# The exact value shown is the maximum over the 95% profile-REML interval of delta.
# Usage (from the repository root): Rscript calibration/28_make_geuvadis_tables.R
# ==============================================================================
N <- "calibration/note"
P <- read.csv("calibration/results/geuvadis_app_paths.csv")
V <- read.csv("calibration/results/geuvadis_app_variants.csv")
ord_gene <- c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266")
meth <- c(wald_h1 = "Wald", score_h0 = "Score", lrt = "RV", mbf_hc = "BF marg.", cbf_deployed = "BF eBIC",
          ols_ebic = "Régression", susie = "SuSiE")
P$key <- paste(P$pop, P$gene)
keys <- unique(P[order(match(P$pop, c("EUR", "YRI")), match(P$gene, ord_gene)), c("key", "pop", "gene", "n", "m")])
rows <- vapply(seq_len(nrow(keys)), function(k) { d <- P[P$key == keys$key[k], ]
  ks <- vapply(names(meth), function(mm) { x <- d$K[d$method == mm]; if (!length(x) || is.na(x)) "--" else sprintf("%d", x) }, "")
  sprintf("%s & %s & %d & %d & %s \\\\", keys$gene[k], keys$pop[k], keys$n[k], keys$m[k], paste(ks, collapse = " & ")) }, "")
pops <- keys$pop
rows <- unlist(lapply(seq_along(rows), function(i) if (i > 1 && pops[i] != pops[i - 1]) c("\\midrule", rows[i]) else rows[i]))
writeLines(c(sprintf("\\begin{tabular}{llrr%s}", strrep("c", length(meth))), "\\toprule",
  sprintf("Gène & pop. & $n$ & $m$ & %s \\\\", paste(meth, collapse = " & ")), "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_geu_paths.tex"))

f1 <- function(x) if (!is.finite(x)) "--" else if (x <= -12) "$<\\!-12$" else sprintf("$%.1f$", x)
tr_lab <- function(s) { p <- strsplit(s, "\\+")[[1]]
  paste(c(lrt = "RV", score_h0 = "S", wald_h1 = "W")[p], collapse = ", ") }
V <- V[order(match(V$pop, c("EUR", "YRI")), match(V$gene, ord_gene), V$step, V$pos), ]
rows <- vapply(seq_len(nrow(V)), function(i) { v <- V[i, ]
  sprintf("%s & %s & %d & %s & %s & $%.1f$ & %s & %s & %s & %s & %s & %s \\\\", v$gene, v$pop, v$step,
          gsub(",", "\\\\,", format(v$pos, big.mark = ",")), tr_lab(v$treatments), v$nu_j,
          f1(v$log10p_nom_wald), f1(v$log10p_exmax_wald), f1(v$log10p_nom_score), f1(v$log10p_exmax_score),
          f1(v$log10p_nom_lrt), f1(v$log10p_exmax_lrt)) }, "")
vp <- V$pop
rows <- unlist(lapply(seq_along(rows), function(i) if (i > 1 && vp[i] != vp[i - 1]) c("\\midrule", rows[i]) else rows[i]))
writeLines(c("\\begin{tabular}{llrrlrcccccc}", "\\toprule",
  "& & & & & & \\multicolumn{2}{c}{Wald, $\\hat\\delta$} & \\multicolumn{2}{c}{Score, $\\hat\\delta_0$} & \\multicolumn{2}{c}{RV} \\\\",
  "\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}",
  "Gène & pop. & pas & position & retenu par & $\\nu_j$ & nominal & exact & nominal & exact & nominal & exact \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_geu_variants.tex"))
cat(readLines(file.path(N, "tab_geu_paths.tex")), sep = "\n"); cat("\n"); cat(readLines(file.path(N, "tab_geu_variants.tex")), sep = "\n")
