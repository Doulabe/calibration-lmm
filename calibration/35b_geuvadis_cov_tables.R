# ==============================================================================
# 35b_geuvadis_cov_tables.R
# Tables of the GEUVADIS application with sample covariates (27c, 33e, 33f + 33b):
#   tab_geu_selection_main.tex  selections of eight procedures (main text)
#   tab_geu_selection.tex       selections of every procedure (supplement)
#   tab_geu_examples.tex        four illustrative variant x test cases (main text)
#   tab_geu_gemma_pvalues.tex   nominal, plug-in and protected p-values of the variants on GEMMA's paths
#   tabA_geu_variants.tex       the same for the paths of the article (candidate-specific kernel)
# p-values are log10; "plug-in" is the p-value from the exact law at delta_hat, "protected" the
# maximum of its Monte Carlo estimates over a grid of delta0 on [0, 1e4], an approximation of the valid
# supremum, replaced (with a dagger) by the controlled value of 33g/27d for the checked cases;
# values are truncated at -12. A circle marks an approximate protected value above -12, whose bound rests
# on a few draws (relative standard error above 0.5 at the point of the bound).
# Usage (from the repository root): ART=<article folder> Rscript calibration/35b_geuvadis_cov_tables.R
# ==============================================================================
suppressPackageStartupMessages(library(jsonlite))
ART <- Sys.getenv("ART", "calibration/output")
RES <- "calibration/results"
wr <- function(lines, f) { writeLines(lines, file.path(ART, "tables", f)); cat("wrote", f, "\n") }
ORD <- c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266")
TAGS <- c(ERAP2 = "ENSG00000164308", RPS26 = "ENSG00000197728", SLFN5 = "ENSG00000166750", SNHG5 = "ENSG00000203875",
          "FLVCR1-AS1" = "ENSG00000198468", "PEX6-region" = "ENSG00000124587", "TRA2A-AS" = "ENSG00000230658", ZNF266 = "ENSG00000174652")
SI <- lapply(readLines(file.path(RES, "geuvadis_comp_cov", "susieinf.jsonl")), fromJSON); names(SI) <- vapply(SI, `[[`, "", "window")
cap <- function(x) if (!is.finite(x)) "--" else if (x <= -12) "$<\\!-12$" else sprintf("$%.1f$", x)   # exact-law values (Monte Carlo)
nomf <- function(x) if (!is.finite(x)) "--" else if (x <= -99) "$<\\!-99$" else sprintf("$%.1f$", x)  # nominal values (closed form)
posf <- function(p) gsub(",", "\\\\,", format(p, big.mark = ","))
nuf <- function(v) if (v >= 1000) formatC(round(v), format = "d", big.mark = "{,}") else if (v >= 100) sprintf("%.0f", v) else sprintf("%.1f", v)

## ---- selections ----------------------------------------------------------------------------------
main_rows <- c(); supp_rows <- c(); long <- list()
for (pop in c("EUR", "YRI")) for (g in ORD) {
  a <- readRDS(file.path(RES, "geuvadis_app_cov", sprintf("%s_%s.rds", pop, TAGS[[g]])))
  b <- readRDS(file.path(RES, "geuvadis_comp_cov", sprintf("%s_%s.rds", pop, TAGS[[g]])))
  Z <- as.matrix(read.table(file.path(RES, "geuvadis_comp", "data", sprintf("%s_%s_Z.txt", pop, g))))
  stopifnot(ncol(Z) == b$m, identical(a$pos, b$pos)); lead <- b$lead_j
  si <- SI[[sprintf("%s_%s", pop, g)]]
  S <- list(lrt = a$sel$lrt, wald = a$sel$wald_h1, null = a$sel$score_h0,
            gw = b$gemma_sel$wald, gl = b$gemma_sel$lrt, gs = b$gemma_sel$score,
            mlmm = if (is.list(b$mlmm)) b$mlmm$extBIC else NA, susie = a$sel$susie, susieinf = unlist(si$susieinf_top),
            bslmm = if (is.list(b$bslmm)) b$bslmm$sel else NA, qtl = if (is.list(b$qtltools)) b$qtltools$sel else NA,
            farmcpu = if (is.list(b$farmcpu)) b$farmcpu$sel else NA, ols = a$sel$ols_ebic)
  cellk <- function(s) { if (length(s) == 1 && is.na(s)) return("err"); if (!length(s)) return("$0$")
    tag <- if (is.na(lead)) FALSE else max(cor(Z[, s, drop = FALSE], Z[, lead])^2) >= 0.8
    sprintf("$%d%s$", length(s), if (tag) "" else "^\\dagger") }
  cells <- vapply(S, cellk, "")
  main_rows <- c(main_rows, sprintf("%s & %s & %d & %s \\\\", g, pop, b$n,
                                    paste(cells[c("lrt", "gw", "gl", "gs", "mlmm", "susie", "susieinf", "qtl")], collapse = " & ")))
  supp_rows <- c(supp_rows, sprintf("%s & %s & %d & %s \\\\", g, pop, b$n, paste(cells, collapse = " & ")))
  for (nm in names(S)) long[[length(long) + 1]] <- data.frame(pop = pop, gene = g, method = nm,
    K = if (length(S[[nm]]) == 1 && is.na(S[[nm]])) NA else length(S[[nm]]),
    positions = if (length(S[[nm]]) && !all(is.na(S[[nm]]))) paste(b$pos[S[[nm]]], collapse = ";") else "")
}
mid <- function(r) unlist(lapply(seq_along(r), function(i) if (i == 9) c("\\midrule", r[i]) else r[i]))
wr(c("\\begin{tabular}{llr*{8}{c}}", "\\toprule",
     "& & & \\shortstack{Candidate-specific\\\\kernel} & \\multicolumn{4}{c}{Regional kernel held fixed} & \\multicolumn{2}{c}{Fine-mapping} & \\\\",
     "\\cmidrule(lr){4-4}\\cmidrule(lr){5-8}\\cmidrule(lr){9-10}",
     "Gene & Pop. & $n$ & LR path & GEMMA W & GEMMA L & GEMMA S & MLMM & SuSiE & SuSiE-inf & QTLtools \\\\",
     "\\midrule", mid(main_rows), "\\bottomrule", "\\end{tabular}"), "tab_geu_selection_main.tex")
wr(c("\\begin{tabular}{llr*{13}{c}}", "\\toprule",
     "& & & \\multicolumn{3}{c}{Candidate-specific kernel} & \\multicolumn{4}{c}{Regional kernel held fixed} & \\multicolumn{3}{c}{Sparse Bayesian} & \\multicolumn{3}{c}{No random effect} \\\\",
     "\\cmidrule(lr){4-6}\\cmidrule(lr){7-10}\\cmidrule(lr){11-13}\\cmidrule(lr){14-16}",
     "Gene & Pop. & $n$ & LR & Wald & Null & GEMMA W & GEMMA L & GEMMA S & MLMM & SuSiE & SuSiE-inf & BSLMM & QTLtools & FarmCPU & Cond. \\\\",
     "\\midrule", mid(supp_rows), "\\bottomrule", "\\end{tabular}"), "tab_geu_selection.tex")
write.csv(do.call(rbind, long), file.path(RES, "geuvadis_cov_selections.csv"), row.names = FALSE)

## ---- variants on GEMMA's paths ----------------------------------------------------------------
V <- do.call(rbind, lapply(list.files(file.path(RES, "geuvadis_comp_cov"), "^(EUR|YRI)_ENSG[0-9]+[.]rds$", full.names = TRUE), function(f) {
  r <- readRDS(f); if (is.null(r$variants)) return(NULL)
  v <- r$variants; v$selected <- vapply(v$j, function(j) any(vapply(r$gemma_sel, function(s) j %in% s, TRUE)), TRUE)
  cbind(pop = r$pop, gene = r$gene, n = r$n, m = r$m, ncov = r$n_cov, v) }))
V$gene <- factor(V$gene, ORD); V$pop <- factor(V$pop, c("EUR", "YRI")); V <- V[order(V$pop, V$gene, V$step, V$pos), ]
## protected values checked on the refined grid (33g_protected_check.R): the controlled value replaces the
## coarse-grid value, and the cell is marked with a dagger
fs_pck <- list.files(file.path(RES, "protected_check"), "[.]rds$", full.names = TRUE)
fs_pck <- c(fs_pck[!grepl("_D1e3", fs_pck)], fs_pck[grepl("_D1e3", fs_pck)])  # decision-domain versions last, so they prevail
## each value with the relative standard error at the point of its bound (above 0.5: the bound rests on a few draws)
rdpc <- function(f) { r <- readRDS(f); cbind(r$info[, c("pop", "gene", "pos", "stat", "step", "controlled")],
                                              relse = r$profile$rel_se[which.max(r$profile$log10up)]) }
PCK <- do.call(rbind, lapply(fs_pck, rdpc))
fs_pca <- list.files(file.path(RES, "protected_check_all"), "[.]rds$", full.names = TRUE)
PCA <- if (length(fs_pca)) do.call(rbind, lapply(fs_pca, rdpc)) else NULL
fs_pcb <- list.files(file.path(RES, "protected_check_bigB"), "[.]rds$", full.names = TRUE)
PCB <- if (length(fs_pcb)) do.call(rbind, lapply(fs_pcb, rdpc)) else NULL
PG <- rbind(PCB, PCA, PCK[PCK$stat %in% c("wald", "kr", "score"), ])   # 100,000 draws (33k, 33l), then 33i (all 38 x 4), then 33g
PG <- PG[!duplicated(PG[, c("pop", "gene", "pos", "step", "stat")]), ]
V$ctl_wald <- V$ctl_kr <- V$ctl_score <- V$ctl_lrt <- FALSE
V$rs_wald <- V$rs_kr <- V$rs_score <- V$rs_lrt <- NA_real_
for (i in seq_len(nrow(PG))) { q <- PG[i, ]
  h <- which(V$pop == q$pop & V$gene == q$gene & V$pos == q$pos & V$step == q$step); if (length(h) != 1) next
  V[h, paste0("log10p_sup_", q$stat)] <- q$controlled; V[h, paste0("ctl_", q$stat)] <- TRUE; V[h, paste0("rs_", q$stat)] <- q$relse }
ALLCTL <- all(V$ctl_wald & V$ctl_kr & V$ctl_score & V$ctl_lrt)        # then no dagger is needed in the table
cat("controlled protected values:", sum(V$ctl_wald) + sum(V$ctl_kr) + sum(V$ctl_score) + sum(V$ctl_lrt), "of", 4 * nrow(V), "\n")
capd <- function(x, flag, rs = NA) paste0(cap(x), if (isTRUE(flag) && !ALLCTL) "$^\\dagger$" else "",
                                          if (isTRUE(rs > 0.5) && x > -12) "$^\\circ$" else "")
lev <- function(n, m, step, ncov, test) { pen <- log(n) + 2 * log(m); nu <- n - ncov - step
  if (test == "lrt") pchisq(pen, 1, lower.tail = FALSE, log.p = TRUE) / log(10) else
    (log(2) + pt(-sqrt((nu - 1) * expm1(pen / n)), nu - 1, log.p = TRUE)) / log(10) }
vrows <- unique(vapply(seq_len(nrow(V)), function(q) { w <- V[q, ]
  sprintf("%s & %s & %d & %s%s & %s & $%s$ & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    w$gene, w$pop, w$step, posf(w$pos), if (w$selected) "" else "$^\\ast$",
    gsub("wald", "W", gsub("lrt", "L", gsub("score", "S", gsub("\\+", ", ", w$tests)))), nuf(w$nu_j),
    nomf(w$log10p_gemma_wald), cap(w$log10p_ex_wald), capd(w$log10p_sup_wald, w$ctl_wald, w$rs_wald),
    nomf(w$log10p_kr_full), cap(w$log10p_ex_kr), capd(w$log10p_sup_kr, w$ctl_kr, w$rs_kr),
    nomf(w$log10p_gemma_lrt), capd(w$log10p_sup_lrt, w$ctl_lrt, w$rs_lrt), nomf(w$log10p_gemma_score),
    capd(w$log10p_sup_score, w$ctl_score, w$rs_score)) }, ""))
wr(c("\\begin{tabular}{llrlcr*{10}{c}}", "\\toprule",
     "& & & & & & \\multicolumn{3}{c}{Wald} & \\multicolumn{3}{c}{Kenward--Roger} & \\multicolumn{2}{c}{Likelihood ratio} & \\multicolumn{2}{c}{Score} \\\\",
     "\\cmidrule(lr){7-9}\\cmidrule(lr){10-12}\\cmidrule(lr){13-14}\\cmidrule(lr){15-16}",
     "Gene & Pop. & Step & Position & Path & $\\nu_j$ & nom. & plug-in & prot. & nom. & plug-in & prot. & nom. & prot. & nom. & prot. \\\\",
     "\\midrule", vrows, "\\bottomrule", "\\end{tabular}"), "tab_geu_gemma_pvalues.tex")
## summaries for the text
fin <- function(x) pmax(x, -12)                              # beyond the resolution: compared at -12
cat("\nGEMMA reproduction: Wald median rel", signif(median(abs(V$log10p_gemma_wald - V$log10p_ours_wald) / abs(V$log10p_gemma_wald)), 2),
    " max abs (p > 1e-12)", signif(max(abs(V$log10p_gemma_wald - V$log10p_ours_wald)[V$log10p_gemma_wald > -12]), 2),
    " LR max", signif(max(abs(V$log10p_gemma_lrt - V$log10p_ours_lrt)), 2), " score max", signif(max(abs(V$log10p_gemma_score - V$log10p_ours_score)), 2), "\n")
d <- with(V, data.frame(gene, pop, step, nu_j, wald = fin(log10p_sup_wald) - fin(log10p_gemma_wald), kr = fin(log10p_sup_kr) - fin(log10p_kr_full),
                        lrt = fin(log10p_sup_lrt) - fin(log10p_gemma_lrt), score = fin(log10p_sup_score) - fin(log10p_gemma_score),
                        kr_plug = fin(log10p_ex_kr) - fin(log10p_kr_full), score_plug = fin(log10p_ex_score) - fin(log10p_gemma_score),
                        wald_plug = fin(log10p_ex_wald) - fin(log10p_gemma_wald)))
cat("protected - nominal (log10, truncated at -12): max wald", round(max(d$wald), 1), " range kr", round(range(d$kr), 2),
    " range lrt", round(range(d$lrt), 2), " range score", round(range(d$score), 2), "\n")
cat("plug-in - nominal: min kr", round(min(d$kr_plug), 1), " min score", round(min(d$score_plug), 1), " max wald", round(max(d$wald_plug), 1), "\n")
RS <- with(V, cbind(wald = rs_wald <= 0.5, kr = rs_kr <= 0.5, lrt = rs_lrt <= 0.5, score = rs_score <= 0.5))
BD <- with(V, cbind(wald = log10p_sup_wald, kr = log10p_sup_kr, lrt = log10p_sup_lrt, score = log10p_sup_score))
LV <- matrix(with(V, mapply(lev, n, m, step, ncov, "wald")), nrow(V), 4)   # Wald-equivalent level of the paths (35d)
cat("resolved bounds only (rel se <= 0.5): max wald", round(max(d$wald[RS[, "wald"]]), 1), " range kr", round(range(d$kr[RS[, "kr"]]), 2),
    " range lrt", round(range(d$lrt[RS[, "lrt"]]), 2), " range score", round(range(d$score[RS[, "score"]]), 2), "\n")
cat("bounds resting on a few draws:", sum(!RS), "of", length(RS), " largest", round(max(BD[!RS]), 2),
    " smallest distance below the level", round(min(LV[!RS] - BD[!RS]), 2), "\n")
cat("untruncated: max (prot - nom) wald", round(max(V$log10p_sup_wald - V$log10p_gemma_wald), 1),
    " max (nom - plug) kr", round(max(V$log10p_kr_full - V$log10p_ex_kr), 1), " max (nom - plug) score", round(max(V$log10p_gemma_score - V$log10p_ex_score), 1), "\n")
print(d[order(-d$wald), ][1:6, ], digits = 3, row.names = FALSE)
write.csv(V, file.path(RES, "geuvadis_cov_gemma_variants.csv"), row.names = FALSE)

## ---- four illustrative cases (main text) -----------------------------------------------------------
ex <- list(
  list(gene = "PEX6-region", pop = "EUR", pos = 42944850L, test = "wald", lab = "Wald (GEMMA)", what = "overstated"),
  list(gene = "ERAP2", pop = "EUR", pos = 96252432L, test = "kr", lab = "Kenward--Roger", what = "close to protected, understated"),
  list(gene = "TRA2A-AS", pop = "YRI", pos = 23134790L, test = "score", lab = "Score (GEMMA)", what = "signal lost"),
  list(gene = "ZNF266", pop = "YRI", pos = 9477868L, test = "wald", lab = "Wald (GEMMA)", what = "overstated, supported"))
erows <- vapply(ex, function(e) { w <- V[V$gene == e$gene & V$pop == e$pop & V$pos == e$pos, ][1, ]
  nom <- switch(e$test, wald = w$log10p_gemma_wald, kr = w$log10p_kr_full, score = w$log10p_gemma_score)
  plg <- switch(e$test, wald = w$log10p_ex_wald, kr = w$log10p_ex_kr, score = w$log10p_ex_score)
  pro <- switch(e$test, wald = w$log10p_sup_wald, kr = w$log10p_sup_kr, score = w$log10p_sup_score)
  rs <- switch(e$test, wald = w$rs_wald, kr = w$rs_kr, score = w$rs_score)
  sprintf("%s & %s & %s & $%s$ & %s & $%.1f$ & %s & %s & %s & %s & %s \\\\", e$gene, e$pop, posf(w$pos), nuf(w$nu_j),
          e$lab, lev(w$n, w$m, w$step, w$ncov, "wald"), nomf(nom), cap(plg), capd(pro, FALSE, rs), nomf(w$log10p_gemma_lrt), e$what) }, "")
wr(c("\\begin{tabular}{lllrlcccccl}", "\\toprule",
     "& & & & & & \\multicolumn{3}{c}{$\\log_{10}p$ of the test} & LR & \\\\", "\\cmidrule(lr){7-9}",
     "& & & & & & Reported & \\multicolumn{2}{c}{Recalibrated} & & \\\\", "\\cmidrule(lr){8-9}",
     "Gene & Pop. & Position & $\\nu_j$ & Test & Level & nominal & plug-in & protected & nominal & Consequence \\\\",
     "\\midrule", erows, "\\bottomrule", "\\end{tabular}"), "tab_geu_examples.tex")

## ---- variants of the paths of the article --------------------------------------------------------
A <- read.csv(file.path(RES, "geuvadis_app_cov_variants.csv"))
A$gene <- factor(A$gene, ORD); A$pop <- factor(A$pop, c("EUR", "YRI")); A <- A[order(A$pop, A$gene, A$step, A$pos), ]
A$ctl_lrt <- FALSE
capa <- function(x, flag) paste0(cap(x), if (isTRUE(flag)) "$^\\dagger$" else "")
PCK_path <- PCK[grepl("path", PCK$stat), ]                          # the decision-domain versions (_D1e3) come last
for (i in seq_len(nrow(PCK_path))) { q <- PCK_path[i, ]
  h <- which(A$pop == q$pop & A$gene == q$gene & A$pos == q$pos & A$step == q$step); stopifnot(length(h) == 1)
  A$log10p_sup_lrt[h] <- q$controlled; A$ctl_lrt[h] <- TRUE }
arows <- vapply(seq_len(nrow(A)), function(q) { w <- A[q, ]
  sprintf("%s & %s & %d & %s & %s & $%s$ & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", w$gene, w$pop, w$step, posf(w$pos),
          gsub("wald_h1", "W", gsub("lrt", "L", gsub("score_h0", "N", gsub("\\+", ", ", w$treatments)))), nuf(w$nu_j),
          nomf(w$log10p_nom_wald), cap(w$log10p_ex_wald), cap(w$log10p_sup_wald),
          nomf(w$log10p_nom_null), cap(w$log10p_ex_null), cap(w$log10p_sup_null),
          nomf(w$log10p_nom_lrt), cap(w$log10p_ex_lrt), capa(w$log10p_sup_lrt, w$ctl_lrt)) }, "")
wr(c("\\begin{tabular}{llrlcr*{9}{c}}", "\\toprule",
     "& & & & & & \\multicolumn{3}{c}{Wald, $t(\\widehat\\delta)$} & \\multicolumn{3}{c}{Null estimate, $t(\\widehat\\delta_0)$} & \\multicolumn{3}{c}{Likelihood ratio} \\\\",
     "\\cmidrule(lr){7-9}\\cmidrule(lr){10-12}\\cmidrule(lr){13-15}",
     "Gene & Pop. & Step & Position & Path & $\\nu_j$ & nom. & plug-in & prot. & nom. & plug-in & prot. & nom. & plug-in & prot. \\\\",
     "\\midrule", arows, "\\bottomrule", "\\end{tabular}"), "tabA_geu_variants.tex")
