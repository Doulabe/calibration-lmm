# ==============================================================================
# 34d_purenull_summary.R
# Pure-noise control of the selection study (29b_export_purenull.R): share of the 500
# datasets (y = noise, four GEUVADIS windows) in which each procedure selects at least one
# variant under its native rule, next to the same share in the 1,000 background-only
# datasets of the main study (selection_summary.csv of 34/34c).
# Sources: 14c_selection_purenull.R (paths of the article, conditional regression, SuSiE),
#          30_comparators_selection.R with COMP_ROOT=comparators_purenull (GEMMA paths,
#          QTLtools, FarmCPU), 31_susieinf_selection.py with the same root (SuSiE-inf).
# Output: results/selection_purenull_summary.csv, <article>/tables/tabS_purenull.tex
# Usage (from the repository root): ART=<article folder> Rscript calibration/34d_purenull_summary.R
# ==============================================================================
suppressPackageStartupMessages(library(jsonlite))
ART <- Sys.getenv("ART", "calibration/output")
RES <- "calibration/results"; ROOT <- file.path(RES, "comparators_purenull")
GENES <- c("ERAP2", "RPS26", "SLFN5", "PEX6-region")
D <- setNames(lapply(GENES, function(s) readRDS(file.path(ROOT, "data", paste0(s, ".rds")))), GENES)
lev_ebic_t <- function(n, m, nu) 2 * pt(-sqrt((nu - 1) * expm1((log(n) + 2 * log(m)) / n)), nu - 1)
prefix_at <- function(pp, lev) { if (is.null(pp) || !nrow(pp)) return(integer(0))
  lev <- rep_len(lev, nrow(pp)); ok <- pp$p <= lev
  if (!ok[1]) integer(0) else pp$idx[seq_len(if (all(ok)) nrow(pp) else which(!ok)[1] - 1)] }
rows <- list(); add <- function(id, method, k) rows[[length(rows) + 1]] <<- data.frame(id = id, method = method, K = k)
## paths of the article, conditional regression, SuSiE
for (f in list.files(file.path(RES, "selection_purenull"), "rds$", full.names = TRUE)) { r <- readRDS(f)
  for (mn in names(r$selections)) { s <- r$selections[[mn]]; if (!(length(s) == 1 && is.na(s))) add(r$id, mn, length(s)) } }
## GEMMA paths, QTLtools, FarmCPU
for (g in GENES) { d <- D[[g]]; n <- nrow(d$Z); m <- ncol(d$Z)
  for (id in d$meta$id) { f <- file.path(ROOT, "r_arms", paste0(id, ".rds")); if (!file.exists(f)) next; r <- readRDS(f)
    if (is.list(r$gemma)) for (tn in c("wald", "lrt", "score")) { pp <- r$gemma[[tn]]
      lev <- if (tn == "lrt") pchisq(log(n) + 2 * log(m), 1, lower.tail = FALSE) else lev_ebic_t(n, m, n - seq_len(max(nrow(pp), 1)))
      add(id, paste0("gemma_", tn), length(prefix_at(pp, lev))) }
    if (is.list(r$qtltools)) add(id, "qtltools", length(r$qtltools[["a0.05"]]))
    if (is.data.frame(r$farmcpu)) add(id, "farmcpu", as.integer(any(r$farmcpu$p <= 0.05 / m))) } }
## SuSiE-inf
for (g in GENES) { f <- file.path(ROOT, "susieinf", paste0(g, "_susie.jsonl")); if (!file.exists(f)) next
  for (l in readLines(f)) { r <- fromJSON(l); add(r$id, "susieinf", length(unlist(r$top))) } }
P <- do.call(rbind, rows)
S <- aggregate(cbind(any = K > 0) ~ method, P, function(v) c(n = length(v), k = sum(v), rate = mean(v)))
S <- data.frame(method = S$method, n = S$any[, "n"], k = S$any[, "k"], rate = S$any[, "rate"])
write.csv(S, file.path(RES, "selection_purenull_summary.csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)
B <- read.csv(file.path(RES, "comparators", "selection_summary.csv"))
bg <- function(mt) { d <- B[B$method == mt & B$rule == "native", ]; if (nrow(d) != 1) return(c(NA, NA)); c(d$null_any, d$n_null) }
pc <- function(x) if (!is.finite(x)) "--" else if (x == 0) "$0$" else if (x < 0.1) sprintf("$%.1f$", 100 * x) else sprintf("$%.0f$", 100 * x)
cp_up <- function(k, n) qbeta(0.975, k + 1, n - k)                   # exact (Clopper-Pearson) upper 95% bound
DEF <- list(c("rv", "lrt", "Likelihood-ratio path", "this study"), c("wald_h1", "wald_h1", "Wald path", "this study"),
            c("score_h0", "score_h0", "Null-estimate path, $t(\\widehat\\delta_0)$", "this study"),
            c("gemma_wald", "gemma_wald", "Wald path", "GEMMA"), c("gemma_lrt", "gemma_lrt", "Likelihood-ratio path", "GEMMA"),
            c("gemma_score", "gemma_score", "Score path", "GEMMA"),
            c("susie", "susie", "Top variant per credible set", "SuSiE"), c("susieinf", "susieinf", "Top variant per credible set", "SuSiE-inf"),
            c("qtltools", "qtltools", "Permutation, 5\\% per gene", "QTLtools"), c("farmcpu", "farmcpu", "Pseudo-QTN model, Bonferroni", "FarmCPU"),
            c("ols_ebic", "ols_ebic", "Conditional regression, eBIC", "this study"))
body <- vapply(DEF, function(d) { b <- bg(d[1]); z <- S[S$method == d[2], ]
  pn <- if (nrow(z)) sprintf("%s {\\scriptsize (%d)} & %s", pc(z$rate), z$n, pc(cp_up(z$k, z$n))) else "-- & --"
  sprintf("%s & %s & %s {\\scriptsize (%d)} & %s \\\\", d[3], d[4], pc(b[1]), as.integer(b[2]), pn) }, "")
writeLines(c("\\begin{tabular}{llccc}", "\\toprule",
             "& & Background only & \\multicolumn{2}{c}{Pure noise} \\\\", "\\cmidrule(lr){3-3}\\cmidrule(lr){4-5}",
             "Procedure & Software & any (\\%) & any (\\%) & upper bound (\\%) \\\\", "\\midrule", body,
             "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabS_purenull.tex"))
cat("wrote tabS_purenull.tex\n")
