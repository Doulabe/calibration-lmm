# ==============================================================================
# 35c_protected_check_table.R
# Supplementary table of the numerical control of the protected p-values
# (33g_protected_check.R for GEMMA's statistics, 27d_protected_check_paths.R for the
# likelihood-ratio path of the article): level, nominal p-value, protected value of the
# coarse grid of 33e/27c and its argmax, maximum over the refined grid extended to
# delta0 = 1e10 and its argmax, and controlled value (largest p + 2 se, the "upper estimate").
# Usage (from the repository root): Rscript calibration/35c_protected_check_table.R
# Output: calibration/output/tables/tabS_protected_check.tex
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/output")
R <- lapply(list.files("calibration/results/protected_check", "[.]rds$", full.names = TRUE), readRDS)
S <- do.call(rbind, lapply(R, function(r) r$info[, c("id", "pop", "gene", "pos", "stat", "step", "level", "nominal",
                                                     "stored_protected", "stored_argmax", "new_max", "new_argmax", "controlled")]))
ord <- c("PEX6_EUR_wald", "ERAP2_EUR_kr", "TRA2A_YRI_score", "ZNF266_YRI_score", "ZNF266_YRI_wald", "PEX6_YRI_wald",
         "PEX6_YRI_lrpath_step2", "FLVCR1_YRI_lrpath_step2")
S$draws <- "4,000"
for (h in grep("_D1e3$", S$id)) {                          # decision-domain versions replace the [0, 1e4] ones
  b <- which(S$id == sub("_D1e3$", "", S$id[h]))
  if (length(b) == 1) { S[b, c("nominal", "new_max", "new_argmax", "controlled")] <- S[h, c("nominal", "new_max", "new_argmax", "controlled")]
    S$id[b] <- paste0(S$id[b], "_keep") } }
S <- S[!grepl("_D1e3$", S$id), ]; S$id <- sub("_keep$", "", S$id)
for (f in list.files("calibration/results/protected_check_bigB", "[.]rds$", full.names = TRUE)) { x <- readRDS(f)$info
  h <- which(S$id == x$id); if (length(h) != 1) next
  S[h, c("new_max", "new_argmax", "controlled")] <- x[, c("new_max", "new_argmax", "controlled")]; S$draws[h] <- "100,000" }
S <- S[order(match(S$id, ord)), ]
## independent replicate with fresh seeds (41_replicate_decisive.R), matched on pop, gene, position, statistic
REP <- do.call(rbind, lapply(list.files("calibration/results/replicate_check", "_rep[.]rds$", full.names = TRUE),
  function(f) readRDS(f)$info[, c("pop", "gene", "pos", "stat", "controlled")]))
S$replicate <- if (is.null(REP)) NA else REP$controlled[match(paste(S$pop, S$gene, S$pos, S$stat), paste(REP$pop, REP$gene, REP$pos, REP$stat))]
lab <- c(wald = "Wald (GEMMA)", kr = "Kenward--Roger", score = "Score (GEMMA)", lrt = "LR (GEMMA)")
f1 <- function(x) sprintf("$%.1f$", x)
fa <- function(x) if (x == 0) "$0$" else sprintf("$10^{%.1f}$", log10(x))
posf <- function(p) gsub(",", "\\\\,", format(p, big.mark = ","))
rows <- vapply(seq_len(nrow(S)), function(i) { s <- S[i, ]
  test <- if (grepl("path", s$stat)) "LR (path of this study)" else lab[[s$stat]]
  sprintf("%s & %s & %s & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", s$gene, s$pop, posf(s$pos), s$step, test,
          f1(s$level), f1(s$nominal), f1(s$stored_protected), fa(s$stored_argmax), s$draws, f1(s$new_max), fa(s$new_argmax),
          f1(s$controlled), if (is.na(s$replicate)) "--" else f1(s$replicate)) }, "")
writeLines(c("\\begin{tabular}{lllrlccccccccc}", "\\toprule",
             "& & & & & & & \\multicolumn{2}{c}{Coarse grid} & \\multicolumn{3}{c}{Refined grid} & \\multicolumn{2}{c}{Upper estimate} \\\\",
             "\\cmidrule(lr){8-9}\\cmidrule(lr){10-12}\\cmidrule(lr){13-14}",
             "Gene & Pop. & Position & Step & Test & Level & nominal & max & at $\\delta_0$ & draws & max & at $\\delta_0$ & first & replicate \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(ART, "tables", "tabS_protected_check.tex"))
cat("wrote tabS_protected_check.tex\n")
print(S, digits = 3, row.names = FALSE)
write.csv(S, "calibration/results/protected_check_summary_all.csv", row.names = FALSE)
