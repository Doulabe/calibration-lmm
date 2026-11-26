# ==============================================================================
# 35_make_article_tables.R
# English tables of the article (exact laws, calibration, selection with
# published software, GEUVADIS EUR/YRI), written to the article's tables/ folder.
# Inputs: results of 06, 08, 14-28, 32-34 in calibration/results.
# Usage (from the repository root): Rscript calibration/35_make_article_tables.R [part]
#   part = calib | select | geuvadis | appendix | all (default)
# ==============================================================================
args <- commandArgs(TRUE); PART <- if (length(args)) args[1] else "all"
RES <- "calibration/results"
TAB <- Sys.getenv("ART_TABLES", "calibration/output/tables")
dir.create(TAB, showWarnings = FALSE)
wr <- function(lines, f) { writeLines(lines, file.path(TAB, f)); cat("wrote", f, "\n") }
## number formatting: two significant digits, powers of ten outside [0.01, 1e4)
fm <- function(v) {
  if (!is.finite(v)) return("--")
  if (v == 0) return("$0$")
  if (v < 0.01 || v >= 1e4) { e <- floor(log10(v)); m <- v / 10^e
    if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
    return(if (abs(round(m, 1) - 1) < 1e-9) sprintf("$10^{%d}$", e) else sprintf("$%.1f\\times10^{%d}$", m, e)) }
  if (v >= 10) return(sprintf("$%.0f$", v))
  sprintf("$%.2f$", v)
}
cellm <- function(v) sprintf("%s {\\scriptsize $[$%s, %s$]$}", fm(median(v)), sub("^\\$|\\$$", "", fm(min(v))) |> (\(x) paste0("$", x, "$"))(),
                             fm(max(v)))
cellr <- function(v) { a <- fm(median(v)); lo <- fm(min(v)); hi <- fm(max(v))
  sprintf("%s {\\scriptsize [%s, %s]}", a, lo, hi) }
## size ratios and powers below 1e-6 rest on a few draws and are censored (as in 37)
fmc <- function(v) if (is.finite(v) && v < 1e-6) "$<10^{-6}$" else fm(v)
cellrc <- function(v) sprintf("%s {\\scriptsize [%s, %s]}", fmc(median(v)), fmc(min(v)), fmc(max(v)))
AL <- c(1e-3, 1e-5, 5e-8)
ALAB <- c("$\\alpha=10^{-3}$", "$\\alpha=10^{-5}$", "$\\alpha=5\\cdot10^{-8}$")

## ---------------------------------------------------------------------------------------------------
if (PART %in% c("calib", "all")) {
  T <- read.csv(file.path(RES, "treatments_pilot_final.csv"))
  M <- read.csv(file.path(RES, "mbf_prior_fix.csv"))
  K <- read.csv(file.path(RES, "kr_exact.csv"))
  E <- read.csv(file.path(RES, "edgeworth_check.csv"))
  T <- T[!is.na(T$alpha), c("config", "design", "nu_delta", "method", "alpha", "size_ratio", "power")]
  M <- M[M$method == "mbf_hc", c("config", "design", "nu_delta", "method", "alpha", "size_ratio", "power")]
  K <- K[K$method %in% c("kr", "satterthwaite"), c("config", "design", "nu_delta", "method", "alpha", "size_ratio", "power")]
  E <- E[E$method %in% c("wald_edgeworth", "lrt_edgeworth"), c("config", "design", "nu_delta", "method", "alpha", "size_ratio", "power")]
  A <- rbind(T, M, K, E)
  rowsdef <- list(
    c("oracle_t", "Known $\\delta$, Student $t_{\\nu-1}$ (reference)"),
    c("wald_h1", "Wald, $\\widehat\\delta$ with the candidate (GEMMA, MM4LMM, GridLMM REML)"),
    c("score_h0", "Score $Q(\\widehat\\delta_0)$, Student reference (as GEMMA's score test)"),
    c("lrt_ml", "Likelihood ratio, $\\delta$ under each model (GEMMA LRT)"),
    c("kr", "Kenward--Roger (pbkrtest, MM4LMM, SAS)"),
    c("satterthwaite", "Satterthwaite degrees of freedom (lmerTest)"),
    c("wald_edgeworth", "Wald, Edgeworth threshold (GEMMA \\texttt{-crt})"),
    c("lrt_edgeworth", "Likelihood ratio, Edgeworth threshold"),
    c("signflip", "Sign-flipping score test (exact)"),
    c("bb_wald", "Berger--Boos $p$-value, exact Wald distribution"),
    c("bf_h1", "Bayes factor at $\\widehat\\delta$"),
    c("bf_h0", "Bayes factor at $\\widehat\\delta_0$"),
    c("mbf", "Marginal Bayes factor, $\\tau^2$ fixed"),
    c("mbf_hc", "Marginal Bayes factor, $\\tau$ half-Cauchy"))
  mk <- function(design, col) {
    d <- A[A$design == design, ]
    body <- vapply(rowsdef, function(rd) { z <- d[d$method == rd[1], ]
      cells <- vapply(AL, function(a) { v <- z[[col]][abs(z$alpha - a) < a * 1e-6]; if (!length(v)) "--" else cellrc(v) }, "")
      sprintf("%s & %s \\\\", rd[2], paste(cells, collapse = " & ")) }, "")
    c(body[1:4], "\\addlinespace", body[5:10], "\\addlinespace", body[11:14])
  }
  hdr <- c("\\begin{tabular}{p{6.6cm}ccc}", "\\toprule", sprintf("Treatment of $\\delta$ & %s \\\\", paste(ALAB, collapse = " & ")), "\\midrule")
  wr(c(hdr, "\\multicolumn{4}{l}{\\emph{GEUVADIS cis windows, published lead, 16 configurations}} \\\\",
       mk("GEUVADIS_cis", "size_ratio"), "\\midrule",
       "\\multicolumn{4}{l}{\\emph{AR(1) blocks with $\\rho=0.95$, 8 configurations}} \\\\",
       mk("AR1_rho0.95", "size_ratio"), "\\bottomrule", "\\end{tabular}"), "tab_calibration.tex")
  wr(c(hdr, "\\multicolumn{4}{l}{\\emph{GEUVADIS cis windows, published lead, 16 configurations}} \\\\",
       mk("GEUVADIS_cis", "power"), "\\midrule",
       "\\multicolumn{4}{l}{\\emph{AR(1) blocks with $\\rho=0.95$, 8 configurations}} \\\\",
       mk("AR1_rho0.95", "power"), "\\bottomrule", "\\end{tabular}"), "tab_power.tex")
  ## corrections by class of nu_delta
  A$cls <- cut(A$nu_delta, c(0, 20, 100, Inf), right = FALSE)
  cls_lab <- c("$\\nu_\\delta<20$", "$20\\le\\nu_\\delta<100$", "$\\nu_\\delta\\ge100$")
  cr <- list(c("wald_h1", "Wald, Student $t_{\\nu-1}$"), c("kr", "Kenward--Roger"), c("satterthwaite", "Satterthwaite"),
             c("wald_edgeworth", "Wald, Edgeworth"), c("lrt_ml", "Likelihood ratio, $\\chi^2_1$"),
             c("lrt_edgeworth", "Likelihood ratio, Edgeworth"), c("mbf_hc", "Marginal Bayes factor, $\\tau$ half-Cauchy"))
  body <- unlist(lapply(seq_along(levels(A$cls)), function(q) {
    d <- A[A$cls == levels(A$cls)[q], ]; nconf <- length(unique(d$config))
    c(if (q > 1) "\\midrule", sprintf("\\multicolumn{7}{l}{%s, %d configurations} \\\\", cls_lab[q], nconf),
      vapply(cr, function(rd) { z <- d[d$method == rd[1], ]
        sz <- vapply(AL, function(a) fm(median(z$size_ratio[abs(z$alpha - a) < a * 1e-6])), "")
        pw <- vapply(AL, function(a) fm(median(z$power[abs(z$alpha - a) < a * 1e-6])), "")
        sprintf("\\quad %s & %s & %s \\\\", rd[2], paste(sz, collapse = " & "), paste(pw, collapse = " & ")) }, "")) }))
  wr(c("\\begin{tabular}{lcccccc}", "\\toprule",
       "& \\multicolumn{3}{c}{Size relative to $\\alpha$} & \\multicolumn{3}{c}{Power} \\\\",
       "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
       "Test & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ & $10^{-3}$ & $10^{-5}$ & $5\\cdot10^{-8}$ \\\\", "\\midrule",
       body, "\\bottomrule", "\\end{tabular}"), "tab_corrections.tex")
  ## nu_j over all SNPs of the windows and versus n
  P <- read.csv(file.path(RES, "nu_delta_per_snp_geuvadis.csv"))
  cat("persnp columns:", names(P), "\n")
  NV <- read.csv(file.path(RES, "nu_vs_n.csv")); cat("nu_vs_n columns:", names(NV), "\n")
}

## ---------------------------------------------------------------------------------------------------
## appendix tables of the working note, translated (numbers unchanged)
if (PART %in% c("appendix", "all")) {
  NOTE <- "calibration/note"
  dict <- c(
    ## full labels first, before their parts are translated (leftovers in the article tables)
    "Retournement de signes sur le score" = "Sign-flipping score test",
    "BF marginal, $\\tau$ demi-Cauchy" = "Marginal BF, $\\tau$ half-Cauchy",
    "Tirages (\\%)" = "Draws (\\%)", "Écart de taille (\\%)" = "Size difference (\\%)",
    "dont saut de $\\hat\\delta_0$" = "at a jump of $\\hat\\delta_0$", "non monotones" = "non-monotone",
    "variance, 1\\textsuperscript{er} ordre" = "variance, first order", "variance exacte" = "exact variance",
    "Taille rapportée au nominal" = "Size relative to nominal", "A priori sur l'effet" = "Prior on the effect",
    "Effet standardisé $u$" = "Standardized effect $u$", "Non monotone" = "Non-monotone",
    "Oracle $t(\\delta_0)$, $\\delta$ connu" = "Oracle $t(\\delta_0)$, $\\delta$ known",
    "Oracle BF($\\delta_0$), $\\delta$ connu" = "Oracle BF($\\delta_0$), $\\delta$ known",
    "Rapport de vraisemblance conditionnel (exact)" = "Conditional likelihood ratio (exact)",
    "Rapport de vraisemblance ML" = "Likelihood ratio (ML)", "Rapport de vraisemblance" = "Likelihood ratio",
    "Retournement de signes" = "Sign flipping",
    "BF en $\\hat\\delta_0$ (bras REML)" = "BF at $\\hat\\delta_0$",
    "BF marginal (principal)" = "Marginal BF, $\\tau^2$ fixed",
    "BF en $\\hat\\delta$ sous $H_1$" = "BF at $\\hat\\delta$ under $H_1$",
    "$\\tau^2=0.04$ fixé" = "$\\tau^2=0.04$ fixed",
    "Échange, groupes de rapport" = "Exchange, groups of ratio",
    "médiane [déciles 1 et 9]" = "median [deciles 1 and 9]", "taille moyenne" = "mean size",
    "C1, $\\tilde\\delta$ fixé, ensemble libéral" = "C1, $\\tilde\\delta$ fixed, liberal set",
    "C2, $\\tilde\\delta$ fixé, ensemble SuSiE" = "C2, $\\tilde\\delta$ fixed, SuSiE set",
    "C3, RV conditionnel, ensemble libéral" = "C3, conditional LR, liberal set",
    "C4, RV conditionnel, ensemble SuSiE" = "C4, conditional LR, SuSiE set",
    "C5, point fixe" = "C5, fixed point", "RV, $\\delta$ réestimé" = "LR, $\\delta$ re-estimated",
    "ML, $\\delta$ ré-estimé" = "ML, $\\delta$ re-estimated", "modèle nul" = "null model",
    "ensemble libéral" = "liberal set", "ensemble SuSiE" = "SuSiE set", "causaux fixés" = "causals fixed",
    "h^2_{\\rm fond}" = "h^2_{\\rm bg}", "$\\tau^2$ fixé" = "$\\tau^2$ fixed", "demi-Cauchy" = "half-Cauchy",
    "premier ordre" = "first order", "Non li\\'e" = "Unlinked", "G\\`ene" = "Gene", "m\\'ediane" = "median",
    "médian" = "median", "na\\\"if" = "naive", "Statistique" = "Statistic", "Traitement" = "Treatment",
    "Procédure" = "Procedure", "Puissance" = "Power", "Classe" = "Class", "Taille" = "Size",
    "Rappel" = "Recall", "rappel" = "recall", "prédit" = "predicted", "départ" = "start", "écart" = "gap",
    "nul (\\%)" = "null (\\%)", "Sous $H_0$" = "Under $H_0$", "sous $H_0$" = "under $H_0$",
    "sous $H_1$" = "under $H_1$", "Plan" = "Design", "libéral" = "liberal", "queue" = "tail",
    "blocs" = "blocks", "& configs &" = "& Configurations &", "réestimé" = "re-estimated", "connu" = "known", "fixé" = "fixed",
    "retenu par" = "selected by", "Régression" = "OLS", "BF marg." = "Marg. BF", "& pas &" = "& step &",
    "Gène" = "Gene")
  tr <- function(x) { for (k in names(dict)) x <- gsub(k, dict[[k]], x, fixed = TRUE)
    x <- gsub("\\bRV\\b", "LR", x, perl = TRUE); x <- gsub("& fond &", "& $h^2_{\\\\rm bg}$ &", x, fixed = FALSE)
    x <- sub("^\\$K\\$ &", "$k_c$ &", x)        # number of causal variants k_c, as in the article (K is the kernel)
    x <- gsub("^nul &", "null &", x); x }
  for (f in c("tab_mono", "tab_laplace", "tab_firstorder", "tab_variance", "tab_synthetic", "tab_geuvadis", "tab_plugin",
              "tab_prior_compare", "tab_prior_diag", "tab_nonmono", "tab_ebic", "tab_clr_size", "tab_clr_power", "tab_clr_nu",
              "tab_ge", "tab_sel_delta", "tab_cond_frontier", "tab_cond_delta", "tab_nu_n", "tab_persnp")) {
    ## tab_geu_paths and tab_geu_variants removed: tabA_geu_variants.tex now comes from 35b
    x <- tr(readLines(file.path(NOTE, paste0(f, ".tex"))))
    left <- grep("[éèàùçÉ]|\\\\'e|\\\\`e", x, value = TRUE)
    if (length(left)) cat("  untranslated in", f, ":", substr(left, 1, 120), sep = "\n    ")
    wr(x, sub("^tab_", "tabA_", paste0(f, ".tex")))
  }
}

## ---------------------------------------------------------------------------------------------------
## GEUVADIS EUR and YRI: selections of every procedure, nominal and exact p-values of GEMMA's tests
if (PART %in% c("geuvadis", "all")) {
  suppressPackageStartupMessages(library(jsonlite))
  ORD <- c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266")
  TAGS <- c(ERAP2 = "ENSG00000164308", RPS26 = "ENSG00000197728", SLFN5 = "ENSG00000166750", SNHG5 = "ENSG00000203875",
            "FLVCR1-AS1" = "ENSG00000198468", "PEX6-region" = "ENSG00000124587", "TRA2A-AS" = "ENSG00000230658", ZNF266 = "ENSG00000174652")
  SI <- lapply(readLines(file.path(RES, "geuvadis_comp", "susieinf.jsonl")), fromJSON)
  names(SI) <- vapply(SI, `[[`, "", "window")
  rows <- c(); vrows <- c(); long <- list()
  for (pop in c("EUR", "YRI")) for (g in ORD) {
    a <- readRDS(file.path(RES, "geuvadis_app", sprintf("%s_%s.rds", pop, TAGS[[g]])))
    b <- readRDS(file.path(RES, "geuvadis_comp", sprintf("%s_%s.rds", pop, TAGS[[g]])))
    Z <- as.matrix(read.table(file.path(RES, "geuvadis_comp", "data", sprintf("%s_%s_Z.txt", pop, g))))
    stopifnot(ncol(Z) == b$m, identical(a$pos, b$pos))
    lead <- b$lead_j
    si <- SI[[sprintf("%s_%s", pop, g)]]
    S <- list(cbf = a$sel$cbf_deployed, lrt = a$sel$lrt, wald = a$sel$wald_h1, score = a$sel$score_h0,
              gw = b$gemma_sel$wald, gl = b$gemma_sel$lrt, gs = b$gemma_sel$score,
              mlmm = if (is.list(b$mlmm)) b$mlmm$extBIC else NA,
              susie = a$sel$susie, susieinf = unlist(si$susieinf_top), bslmm = if (is.list(b$bslmm)) b$bslmm$sel else NA,
              qtl = if (is.list(b$qtltools)) b$qtltools$sel else NA,
              farmcpu = if (is.list(b$farmcpu)) b$farmcpu$sel else NA, ols = a$sel$ols_ebic)   # order of the header
    cellk <- function(s) { if (length(s) == 1 && is.na(s)) return("err")
      if (!length(s)) return("$0$")
      tag <- if (is.na(lead)) FALSE else max(cor(Z[, s, drop = FALSE], Z[, lead])^2) >= 0.8
      sprintf("$%d%s$", length(s), if (tag) "" else "^\\dagger") }
    rows <- c(rows, sprintf("%s & %s & %d & %s \\\\", g, pop, b$n, paste(vapply(S, cellk, ""), collapse = " & ")))
    for (nm in names(S)) long[[length(long) + 1]] <- data.frame(pop = pop, gene = g, method = nm,
      K = if (length(S[[nm]]) == 1 && is.na(S[[nm]])) NA else length(S[[nm]]),
      positions = if (length(S[[nm]]) && !all(is.na(S[[nm]]))) paste(b$pos[S[[nm]]], collapse = ";") else "")
    ## exact p-values of GEMMA's statistics for the variants on GEMMA's paths
    if (!is.null(b$variants)) { v <- b$variants; v <- v[order(v$step), ]
      cap <- function(x) if (!is.finite(x)) "--" else if (x <= -12) "$<\\!-12$" else sprintf("$%.1f$", x)
      for (q in seq_len(nrow(v))) { w <- v[q, ]
        sel_any <- any(vapply(b$gemma_sel, function(s) w$j %in% s, TRUE))
        vrows <- c(vrows, sprintf("%s & %s & %d & %s%s & %s & $%s$ & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          g, pop, w$step, gsub(",", "\\\\,", format(w$pos, big.mark = ",")), if (sel_any) "" else "$^\\ast$",
          gsub("wald", "W", gsub("lrt", "L", gsub("score", "S", gsub("\\+", ", ", w$tests)))),
          if (w$nu_j >= 100) sprintf("%.0f", w$nu_j) else sprintf("%.1f", w$nu_j),
          cap(w$log10p_gemma_wald), cap(w$log10p_exmax_wald), cap(w$log10p_kr_full), cap(w$log10p_exmax_kr),
          cap(w$log10p_gemma_lrt), cap(w$log10p_exmax_lrt), cap(w$log10p_gemma_score), cap(w$log10p_exmax_score))) } }
  }
  rows <- unlist(lapply(seq_along(rows), function(i) if (i == 9) c("\\midrule", rows[i]) else rows[i]))
  wr(c("\\begin{tabular}{llr*{14}{c}}", "\\toprule",
       "& & & \\multicolumn{4}{c}{Candidate-specific kernel} & \\multicolumn{4}{c}{Regional kernel held fixed} & \\multicolumn{3}{c}{Sparse Bayesian} & \\multicolumn{3}{c}{No random effect} \\\\",
       "\\cmidrule(lr){4-7}\\cmidrule(lr){8-11}\\cmidrule(lr){12-14}\\cmidrule(lr){15-17}",
       "Gene & Pop. & $n$ & BF eBIC & LR & Wald & Score & GEMMA W & GEMMA L & GEMMA S & MLMM & SuSiE & SuSiE-inf & BSLMM & QTLtools & FarmCPU & OLS \\\\",
       "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_geu_selection.tex")
  write.csv(do.call(rbind, long), file.path(RES, "geuvadis_comp_selections.csv"), row.names = FALSE)
  vrows <- unique(vrows)
  wr(c("\\begin{tabular}{llrlcr*{8}{c}}", "\\toprule",
       "& & & & & & \\multicolumn{2}{c}{Wald} & \\multicolumn{2}{c}{Kenward--Roger} & \\multicolumn{2}{c}{Likelihood ratio} & \\multicolumn{2}{c}{Score} \\\\",
       "\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}\\cmidrule(lr){11-12}\\cmidrule(lr){13-14}",
       "Gene & Pop. & Step & Position & Path & $\\nu_j$ & nominal & exact & nominal & exact & nominal & exact & nominal & exact \\\\",
       "\\midrule", vrows, "\\bottomrule", "\\end{tabular}"), "tab_geu_gemma_pvalues.tex")
}

## ---------------------------------------------------------------------------------------------------
## selection-level study with published software (34, 34b, 34c)
if (PART %in% c("select", "all")) {
  S <- read.csv(file.path(RES, "comparators", "selection_summary.csv"), stringsAsFactors = FALSE)
  F <- read.csv(file.path(RES, "comparators", "selection_frontier.csv"), stringsAsFactors = FALSE)
  pick <- function(mt, rule = "native", thr = NULL) { d <- S[S$method == mt & S$rule == rule, ]
    if (!is.null(thr)) d <- d[d$thr == thr, ]; stopifnot(nrow(d) == 1); d }
  pc <- function(x) if (!is.finite(x)) "--" else if (x == 0) "$0$" else if (x < 0.1) sprintf("$%.1f$", 100 * x) else sprintf("$%.0f$", 100 * x)
  f2 <- function(x) if (!is.finite(x)) "--" else sprintf("$%.2f$", x)
  fr <- function(mt, col) { v <- F[[col]][F$method == mt]; if (!length(v) || !is.finite(v)) "--" else sprintf("$%.2f$", v) }
  row <- function(lab, sw, mt, rule = "native", thr = NULL, front = TRUE) { d <- pick(mt, rule, thr)
    sprintf("%s & %s & %s {\\scriptsize (%d)} & %s & %s & %s & %s & %s \\\\", lab, sw, pc(d$null_any), d$n_null,
            f2(d$recall), f2(d$precision), pc(d$any_off), if (front) fr(mt, "recall_at_1") else "--",
            if (front) fr(mt, "recall_at_5") else "--") }
  grp <- function(t) sprintf("\\multicolumn{8}{l}{\\emph{%s}} \\\\", t)
  body <- c(grp("Regional mixed model, candidate-specific kernel (this study)"),
    row("Likelihood ratio, $\\delta$ under each model", "this study", "rv"),
    row("Wald, $\\widehat\\delta$ with the candidate", "this study", "wald_h1", front = FALSE),
    row("Null estimate, $t(\\widehat\\delta_0)$", "this study", "score_h0", front = FALSE),
    row("Oracle, $\\delta$ known", "reference", "oracle"),
    "\\addlinespace", grp("Regional kernel held fixed along the path"),
    row("Wald path, eBIC level", "GEMMA", "gemma_wald"),
    row("Wald path, Bonferroni $0.05/m$", "GEMMA", "gemma_wald", "bonferroni", front = FALSE),
    row("Likelihood-ratio path, eBIC level", "GEMMA", "gemma_lrt"),
    row("Likelihood-ratio path, Bonferroni", "GEMMA", "gemma_lrt", "bonferroni", front = FALSE),
    row("Score path, eBIC level", "GEMMA", "gemma_score"),
    row("Forward--backward, extended BIC", "MLMM", "mlmm"),
    row("Forward--backward, multiple Bonferroni", "MLMM", "mlmm_mbonf", front = FALSE),
    "\\addlinespace", grp("Bayesian fine-mapping"),
    row("SuSiE, top variant of each 95\\% set", "susieR", "susie", front = FALSE),
    row("SuSiE-inf, top variant of each 95\\% set", "fine-mapping-inf", "susieinf"),
    row("BSLMM, PIP $>0.5$", "GEMMA", "bslmm"),
    "\\addlinespace", grp("No regional random effect"),
    row("Conditional regression, eBIC", "this study", "ols_ebic", front = FALSE),
    row("Forward--backward, permutation, 5\\% per gene", "QTLtools", "qtltools"),
    row("Pseudo-QTN model, Bonferroni", "FarmCPU (rMVP)", "farmcpu"),
    "\\addlinespace", grp("Variance ratio conditional on a provisional signal set"),
    row("Conditional LR, SuSiE set (C4)", "this study", "c4"),
    row("$\\tilde\\delta$ fixed, SuSiE set (C2)", "this study", "c2"))
  wr(c("\\begin{tabular}{p{5.4cm}lcccccc}", "\\toprule",
       "& & Background only & \\multicolumn{3}{c}{Datasets with causal variants} & \\multicolumn{2}{c}{Recall at a background-only rate of} \\\\",
       "\\cmidrule(lr){3-3}\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}",
       "Procedure & Software & any (\\%) & recall & precision & off-target (\\%) & $1\\%$ & $5\\%$ \\\\",
       "\\midrule", body, "\\bottomrule", "\\end{tabular}"), "tab_selection.tex")
}


## ---------------------------------------------------------------------------------------------------
## reduced selection table for the main text of the restructured article
if (PART %in% c("select_main", "all")) {
  S <- read.csv(file.path(RES, "comparators", "selection_summary.csv"), stringsAsFactors = FALSE)
  F <- read.csv(file.path(RES, "comparators", "selection_frontier.csv"), stringsAsFactors = FALSE)
  pick <- function(mt, rule = "native") { d <- S[S$method == mt & S$rule == rule, ]; stopifnot(nrow(d) == 1); d }
  pc <- function(x) if (!is.finite(x)) "--" else if (x == 0) "$0$" else if (x < 0.1) sprintf("$%.1f$", 100 * x) else sprintf("$%.0f$", 100 * x)
  f2 <- function(x) if (!is.finite(x)) "--" else sprintf("$%.2f$", x)
  fr <- function(mt) { v <- F$recall_at_1[F$method == mt]; if (!length(v) || !is.finite(v)) "--" else sprintf("$%.2f$", v) }
  row <- function(lab, sw, dt, mt, front = TRUE) { d <- pick(mt)
    sprintf("%s & %s & %s & %s {\\scriptsize (%d)} & %s & %s & %s \\\\", lab, sw, dt, pc(d$null_any), d$n_null,
            f2(d$recall), f2(d$precision), if (front) fr(mt) else "--") }
  body <- c(
    row("Likelihood-ratio path, candidate-specific kernel", "Section~\\ref{sec:paths}", "each model", "rv"),
    row("Wald path", "GEMMA", "with candidate", "gemma_wald"),
    row("Likelihood-ratio path", "GEMMA", "each model", "gemma_lrt"),
    row("Score path", "GEMMA", "null model", "gemma_score"),
    row("Forward--backward, extended BIC", "MLMM", "null, each step", "mlmm"),
    row("Top variant per credible set", "SuSiE", "no background", "susie", front = FALSE),
    row("Top variant per credible set", "SuSiE-inf", "moments, plug-in", "susieinf"),
    row("Permutation, 5\\% per gene", "QTLtools", "no background", "qtltools", front = FALSE),
    "\\addlinespace",
    row("Oracle", "reference", "known", "oracle"))
  wr(c("\\begin{tabular}{lllcccc}", "\\toprule",
       "& & & & \\multicolumn{2}{c}{Causal variants} & \\\\",
       "\\cmidrule(lr){5-6}",
       "Procedure & Implementation & $\\delta$ & Background only (\\%) & recall & precision & Recall at $1\\%$ \\\\",
       "\\midrule", body, "\\bottomrule", "\\end{tabular}"), "tab_selection_main.tex")
}

## ---------------------------------------------------------------------------------------------------
## GEUVADIS tables for the main text of the restructured article
if (PART %in% c("geu_main", "all")) {
  suppressPackageStartupMessages(library(jsonlite))
  ORD <- c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266")
  TAGS <- c(ERAP2 = "ENSG00000164308", RPS26 = "ENSG00000197728", SLFN5 = "ENSG00000166750", SNHG5 = "ENSG00000203875",
            "FLVCR1-AS1" = "ENSG00000198468", "PEX6-region" = "ENSG00000124587", "TRA2A-AS" = "ENSG00000230658", ZNF266 = "ENSG00000174652")
  SI <- lapply(readLines(file.path(RES, "geuvadis_comp", "susieinf.jsonl")), fromJSON); names(SI) <- vapply(SI, `[[`, "", "window")
  rows <- c()
  for (pop in c("EUR", "YRI")) for (g in ORD) {
    a <- readRDS(file.path(RES, "geuvadis_app", sprintf("%s_%s.rds", pop, TAGS[[g]])))
    b <- readRDS(file.path(RES, "geuvadis_comp", sprintf("%s_%s.rds", pop, TAGS[[g]])))
    Z <- as.matrix(read.table(file.path(RES, "geuvadis_comp", "data", sprintf("%s_%s_Z.txt", pop, g)))); lead <- b$lead_j
    S <- list(cbf = a$sel$cbf_deployed, lrt = a$sel$lrt, gw = b$gemma_sel$wald, gl = b$gemma_sel$lrt, gs = b$gemma_sel$score,
              mlmm = if (is.list(b$mlmm)) b$mlmm$extBIC else NA, susie = a$sel$susie,
              susieinf = unlist(SI[[sprintf("%s_%s", pop, g)]]$susieinf_top), qtl = if (is.list(b$qtltools)) b$qtltools$sel else NA)
    cellk <- function(s) { if (length(s) == 1 && is.na(s)) return("err"); if (!length(s)) return("$0$")
      tag <- if (is.na(lead)) FALSE else max(cor(Z[, s, drop = FALSE], Z[, lead])^2) >= 0.8
      sprintf("$%d%s$", length(s), if (tag) "" else "^\\dagger") }
    rows <- c(rows, sprintf("%s & %s & %d & %s \\\\", g, pop, b$n, paste(vapply(S, cellk, ""), collapse = " & ")))
  }
  rows <- unlist(lapply(seq_along(rows), function(i) if (i == 9) c("\\midrule", rows[i]) else rows[i]))
  wr(c("\\begin{tabular}{llr*{9}{c}}", "\\toprule",
       "& & & \\multicolumn{2}{c}{This paper} & \\multicolumn{4}{c}{Regional kernel held fixed} & \\multicolumn{2}{c}{Fine-mapping} & \\\\",
       "\\cmidrule(lr){4-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-11}",
       "Gene & Pop. & $n$ & BF eBIC & LR & GEMMA W & GEMMA L & GEMMA S & MLMM & SuSiE & SuSiE-inf & QTLtools \\\\",
       "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), "tab_geu_selection_main.tex")
  ## four illustrative variants, GEMMA's kernel and covariates
  V <- do.call(rbind, lapply(list.files(file.path(RES, "geuvadis_comp"), "^(EUR|YRI)_ENSG[0-9]+[.]rds$", full.names = TRUE), function(f) {
    r <- readRDS(f); if (is.null(r$variants)) return(NULL); v <- r$variants; v$ml_fixed <- NULL
    cbind(pop = r$pop, gene = r$gene, n = r$n, m = r$m, v) }))
  lev <- function(n, m, step, test) { pen <- log(n) + 2 * log(m); nu <- n - step
    if (test == "lrt") pchisq(pen, 1, lower.tail = FALSE, log.p = TRUE) / log(10) else
      (log(2) + pt(-sqrt((nu - 1) * expm1(pen / n)), nu - 1, log.p = TRUE)) / log(10) }
  ex <- list(
    list(gene = "PEX6-region", pop = "EUR", pos = 42944850L, test = "wald", lab = "Wald (GEMMA)", what = "overstates the evidence"),
    list(gene = "ERAP2", pop = "EUR", pos = 96252432L, test = "kr", lab = "Kenward--Roger", what = "understates the evidence"),
    list(gene = "TRA2A-AS", pop = "YRI", pos = 23134790L, test = "score", lab = "Score (GEMMA)", what = "signal lost"),
    list(gene = "SLFN5", pop = "YRI", pos = 33520965L, test = "wald", lab = "Wald (GEMMA)", what = "variant added"))
  f1 <- function(x) if (!is.finite(x)) "--" else if (x <= -12) "$<\\!-12$" else sprintf("$%.1f$", x)
  erows <- vapply(ex, function(e) { w <- V[V$gene == e$gene & V$pop == e$pop & V$pos == e$pos, ][1, ]
    nom <- switch(e$test, wald = w$log10p_gemma_wald, kr = w$log10p_kr_full, score = w$log10p_gemma_score)
    exa <- switch(e$test, wald = w$log10p_exmax_wald, kr = w$log10p_exmax_kr, score = w$log10p_exmax_score)
    sprintf("%s & %s & %d & %s & $%s$ & %s & $%.1f$ & $%.1f$ & %s & %s \\\\", e$gene, e$pop, w$step,
            gsub(",", "\\\\,", format(w$pos, big.mark = ",")), if (w$nu_j >= 100) sprintf("%.0f", w$nu_j) else sprintf("%.1f", w$nu_j),
            e$lab, lev(w$n, w$m, w$step, e$test), nom, f1(exa), e$what) }, "")
  wr(c("\\begin{tabular}{llrlrlcccl}", "\\toprule",
       "Gene & Pop. & Step & Position & $\\nu_j$ & Test & Level & Nominal & Exact & Consequence \\\\",
       "\\midrule", erows, "\\bottomrule", "\\end{tabular}"), "tab_geu_examples.tex")
}
