# ==============================================================================
# 35d_gemma_decisions_controlled.R
# Decisions of GEMMA's forward paths in the GEUVADIS application (33e) compared with the
# controlled protected p-values of the same statistics (33i_protected_check_all.R). For every
# entry of a path (selected variants and the first variant that is not selected), the nominal
# p-value of the path's test and its controlled protected value are compared with the
# per-test level of the extended BIC used by 33e (Wald-equivalent level for the three tests):
#   supported      nominal below the level and refined bound (largest p_hat + 2 se) below it
#   not supported  nominal below, and the largest p_hat - 2 se above the level
#   undetermined   nominal below, the level between the two (numerically undetermined)
#   reference / confirmed / undetermined stop: the same for the first variant that a path does not select
#   (reference = stop due to the reference distribution)
# Usage (from the repository root): Rscript calibration/35d_gemma_decisions_controlled.R
# Output: results/gemma_decisions_controlled.csv, calibration/output/tables/tabS_decisions.tex
# ==============================================================================
ART <- Sys.getenv("ART", "calibration/output")
RES <- "calibration/results"
## refined profiles: 100,000 draws (33k, 33l) where available, else 4,000 (33i). For each case, the bound is the largest
## p_hat + 2 se and the lower value the largest p_hat - 2 se over the grid of delta0 (relative se < 0.5)
summ <- function(r) { p <- r$profile; lo <- ifelse(p$rel_se < 0.5, p$log10p + log10(pmax(1 - 2 * p$rel_se, 1e-12)), -Inf)
  data.frame(r$info[, c("pop", "gene", "pos", "step", "stat")], controlled = max(p$log10up), lower = max(lo)) }
PCA <- do.call(rbind, lapply(list.files(file.path(RES, "protected_check_all"), "[.]rds$", full.names = TRUE), function(f) summ(readRDS(f))))
for (f in list.files(file.path(RES, "protected_check_bigB"), "[.]rds$", full.names = TRUE)) { x <- summ(readRDS(f))
  h <- which(PCA$pop == x$pop & PCA$gene == x$gene & PCA$pos == x$pos & PCA$step == x$step & PCA$stat == x$stat)
  if (length(h) == 1) PCA[h, c("controlled", "lower")] <- x[, c("controlled", "lower")] }
out <- list()
for (f in list.files(file.path(RES, "geuvadis_comp_cov"), "^(EUR|YRI)_ENSG[0-9]+[.]rds$", full.names = TRUE)) {
  r <- readRDS(f); if (is.null(r$variants)) next
  n <- r$n; m <- r$m; ncf <- 1L + r$n_cov; pen <- log(n) + 2 * log(m)
  lev_wald <- function(nu) log10(2 * pt(-sqrt((nu - 1) * expm1(pen / n)), nu - 1))
  for (tn in c("wald", "lrt", "score")) {
    s <- r$gemma_paths[[tn]]$idx; nsel <- length(r$gemma_sel[[tn]])
    s <- s[seq_len(min(length(s), nsel + 1L))]
    for (k in seq_along(s)) {
      pos <- r$pos[s[k]]; v <- r$variants[r$variants$pos == pos & r$variants$step == k, ]
      if (nrow(v) != 1) next
      nom <- v[[paste0("log10p_gemma_", tn)]]
      q <- PCA[PCA$pop == r$pop & PCA$gene == r$gene & PCA$pos == pos & PCA$step == k & PCA$stat == tn, ]
      if (nrow(q) != 1) next
      lv <- lev_wald(n - ncf - k + 1L); sel <- k <= nsel
      cls <- if (nom <= lv) { if (q$controlled <= lv) "supported" else if (q$lower > lv) "not supported" else "undetermined" } else
             { if (q$controlled <= lv) "reference" else if (q$lower > lv) "confirmed" else "undetermined stop" }
      out[[length(out) + 1]] <- data.frame(pop = r$pop, gene = r$gene, path = tn, step = k, pos = pos, selected = sel,
                                           level = lv, nominal = nom, protected = q$controlled, lower = q$lower, class = cls)
    }
  }
}
D <- do.call(rbind, out)
write.csv(D, file.path(RES, "gemma_decisions_controlled.csv"), row.names = FALSE)
print(table(D$path, D$class, D$pop))
print(D[!D$class %in% c("supported", "confirmed"), ], digits = 3, row.names = FALSE)
lab <- c(wald = "Wald", lrt = "Likelihood ratio", score = "Score")
rows <- unlist(lapply(c("EUR", "YRI"), function(pp) vapply(c("wald", "lrt", "score"), function(tn) {
  z <- D[D$pop == pp & D$path == tn, ]
  sprintf("%s & %s & %d & %d & %d & %d & %d & %d & %d \\\\", pp, lab[[tn]], sum(z$selected), sum(z$class == "supported"),
          sum(z$class == "not supported"), sum(z$class == "undetermined"), sum(z$class == "confirmed"),
          sum(z$class == "reference"), sum(z$class == "undetermined stop")) }, "")))
writeLines(c("\\begin{tabular}{llccccccc}", "\\toprule",
             "& & & \\multicolumn{3}{c}{Selections} & \\multicolumn{3}{c}{Stops} \\\\",
             "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
             "Pop. & Path & Selected & supported & not supported & undetermined & confirmed & reference & undetermined \\\\",
             "\\midrule", rows[1:3], "\\midrule", rows[4:6], "\\bottomrule", "\\end{tabular}"),
           file.path(ART, "tables", "tabS_decisions.tex"))
cat("wrote tabS_decisions.tex\n")
