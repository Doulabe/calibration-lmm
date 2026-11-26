# ==============================================================================
# 15_make_selection_tables.R
# Final tables of the selection-level demonstration (14_selection_demo.R, plus the
# delta-known arm of 16_selection_oracle.R when its CSV exists).
#   tab_sel_null.tex   global null: P(>= 1 selection) and mean number of selections
#   tab_sel_cells.tex  recall / precision per scenario and procedure
#   results/selection_demo_paired.txt  paired differences between procedures
# Usage (from the repository root): Rscript calibration/15_make_selection_tables.R
# ==============================================================================
R <- read.csv("calibration/results/selection_demo_replicates.csv")
f_or <- "calibration/results/selection_oracle_replicates.csv"
if (file.exists(f_or)) R <- rbind(R, read.csv(f_or)[, names(R)])
N <- "calibration/note"
cat("replicates:", length(unique(paste(R$cell, R$b))), " methods:", paste(unique(R$method), collapse = ", "), "\n")
lab <- c(wald_h1 = "Wald, $\\hat\\delta$ sous $H_1$", lrt = "Rapport de vraisemblance",
         mbf_fixed = "BF marginal, $\\tau^2$ fixé", mbf_hc = "BF marginal, $\\tau$ demi-Cauchy",
         score_h0 = "Score, $\\hat\\delta_0$ sous $H_0$", cbf_deployed = "BF marginal, eBIC",
         oracle_delta = "$\\delta$ connu (oracle)",
         ols_ebic = "Régression sans fond", susie = "SuSiE")
blocks <- list(c("wald_h1", "lrt", "mbf_fixed", "mbf_hc", "score_h0", "cbf_deployed"),
               "oracle_delta", c("ols_ebic", "susie"))
blocks <- Filter(length, lapply(blocks, intersect, unique(R$method)))
f2 <- function(v) sprintf("%.2f", v)
with_rules <- function(rows_by_block) head(unlist(lapply(rows_by_block, c, "\\midrule")), -1)

## ---- global null ----------------------------------------------------------------------
z <- R[R$Kc == 0, ]
null_row <- function(m) {
  cells <- vapply(c(0.10, 0.25), function(h) { w <- z[z$method == m & abs(z$h2bg - h) < 1e-9, ]
    p <- mean(w$any_off); se <- sqrt(p * (1 - p) / nrow(w))
    sprintf("$%s$ {\\scriptsize ($%s$)} & $%s$", f2(p), f2(se), f2(mean(w$fp025))) }, "")
  sprintf("%s & %s \\\\", lab[m], paste(cells, collapse = " & ")) }
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule",
  "& \\multicolumn{2}{c}{$h^2_{\\rm fond}=0.10$} & \\multicolumn{2}{c}{$h^2_{\\rm fond}=0.25$} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Procédure & $P(\\ge1)$ & nombre & $P(\\ge1)$ & nombre \\\\", "\\midrule",
  with_rules(lapply(blocks, function(b) vapply(b, null_row, ""))), "\\bottomrule", "\\end{tabular}"),
  file.path(N, "tab_sel_null.tex"))

## ---- scenarios with causal variants: recall / precision -----------------------------------
sc <- unique(R[R$Kc > 0, c("h2bg", "Kc", "h2c")]); sc <- sc[order(sc$h2c, sc$Kc, sc$h2bg), ]
## one row per scenario, one column per procedure (MBF columns equal the LRT ones, see paired file)
short <- c(wald_h1 = "Wald $H_1$", lrt = "RV", score_h0 = "Score $H_0$", cbf_deployed = "BF eBIC",
           oracle_delta = "oracle", ols_ebic = "régression", susie = "SuSiE")
cols <- intersect(names(short), unique(R$method)); n_mix <- sum(cols %in% names(short)[1:5])
rows <- apply(sc, 1, function(s) {
  cells <- vapply(cols, function(m) { w <- R[R$method == m & R$Kc == s["Kc"] &
      abs(R$h2c - s["h2c"]) < 1e-9 & abs(R$h2bg - s["h2bg"]) < 1e-9, ]
    sprintf("%s/%s", f2(mean(w$recall, na.rm = TRUE)), f2(mean(w$precision, na.rm = TRUE))) }, "")
  sprintf("$%d$ & $%.2f$ & $%.2f$ & %s \\\\", as.integer(s["Kc"]), s["h2c"], s["h2bg"], paste(cells, collapse = " & ")) })
writeLines(c(sprintf("\\begin{tabular}{rcc%s}", strrep("c", length(cols))), "\\toprule",
  sprintf("& & & \\multicolumn{%d}{c}{Modèle mixte} & \\multicolumn{%d}{c}{Sans modèle du fond} \\\\", n_mix, length(cols) - n_mix),
  sprintf("\\cmidrule(lr){4-%d}\\cmidrule(lr){%d-%d}", 3 + n_mix, 4 + n_mix, 3 + length(cols)),
  sprintf("$K$ & $h^2_c$ & fond & %s \\\\", paste(short[cols], collapse = " & ")), "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_sel_cells.tex"))

## ---- paired differences --------------------------------------------------------------------
W <- reshape(R[, c("cell", "b", "Kc", "h2bg", "h2c", "method", "recall", "any_off", "K_hat")],
             idvar = c("cell", "b", "Kc", "h2bg", "h2c"), timevar = "method", direction = "wide")
nn <- W[W$Kc > 0, ]
pd <- function(a, b, col = "recall", d = nn) { x <- d[[paste0(col, ".", a)]] - d[[paste0(col, ".", b)]]
  sprintf("%+.3f (SE %.3f, n = %d)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))), sum(!is.na(x))) }
out <- c(
  paste("recall LRT - score_h0     :", pd("lrt", "score_h0")),
  paste("recall wald_h1 - LRT      :", pd("wald_h1", "lrt")),
  paste("recall LRT - cbf_deployed :", pd("lrt", "cbf_deployed")),
  paste("recall mbf_hc - mbf_fixed :", pd("mbf_hc", "mbf_fixed")),
  paste("recall mbf_hc - LRT       :", pd("mbf_hc", "lrt")),
  paste("recall susie - LRT        :", pd("susie", "lrt")),
  paste("any_off wald_h1 - LRT     :", pd("wald_h1", "lrt", "any_off")),
  paste("identical K_hat mbf_hc / lrt (all reps):",
        sprintf("%.3f", mean(W$K_hat.mbf_hc == W$K_hat.lrt))))
if ("recall.oracle_delta" %in% names(W)) {
  out <- c(out, paste("recall oracle - LRT       :", pd("oracle_delta", "lrt")),
           paste("recall susie - oracle     :", pd("susie", "oracle_delta")))
  for (h in c(0.08, 0.15)) out <- c(out, sprintf("h2c %.2f: recall oracle %.3f, lrt %.3f, susie %.3f; K_hat=0 oracle %.2f, lrt %.2f",
    h, mean(nn$recall.oracle_delta[nn$h2c == h]), mean(nn$recall.lrt[nn$h2c == h]), mean(nn$recall.susie[nn$h2c == h]),
    mean(nn$K_hat.oracle_delta[nn$h2c == h] == 0), mean(nn$K_hat.lrt[nn$h2c == h] == 0)))
}
for (h in c(0.10, 0.25)) {
  z0 <- W[W$Kc == 0 & abs(W$h2bg - h) < 1e-9, ]
  ms <- intersect(c("wald_h1", "lrt", "score_h0", "oracle_delta", "ols_ebic", "susie"), unique(R$method))
  out <- c(out, sprintf("null bg %.2f: any_off %s (n=%d)", h,
    paste(sprintf("%s %.3f", ms, vapply(ms, function(m) mean(z0[[paste0("any_off.", m)]]), 0)), collapse = ", "), nrow(z0)))
}
writeLines(out, "calibration/results/selection_demo_paired.txt"); cat(out, sep = "\n")

## ---- absorption of the unselected causals into delta_hat (17_delta_inflation_check.R) ----------
f_d <- "calibration/results/selection_delta_inflation.csv"
if (file.exists(f_d) && "recall.oracle_delta" %in% names(W)) {
  D <- read.csv(f_d); D <- D[D$Kc > 0, ]
  rows <- apply(sc, 1, function(s) {
    d <- D[D$Kc == s["Kc"] & abs(D$h2c - s["h2c"]) < 1e-9 & abs(D$h2bg - s["h2bg"]) < 1e-9, ]
    w <- nn[nn$Kc == s["Kc"] & abs(nn$h2c - s["h2c"]) < 1e-9 & abs(nn$h2bg - s["h2bg"]) < 1e-9, ]
    sprintf("$%d$ & $%.2f$ & $%.2f$ & $%.1f$ & $%.1f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%+.2f$ \\\\",
            as.integer(s["Kc"]), s["h2c"], s["h2bg"], median(d$delta_absorb / d$delta_star),
            median(d$d0_reml / d$delta_star), median(d$dc_reml / d$delta_star),
            mean(w$recall.oracle_delta), mean(w$recall.lrt), mean(w$recall.oracle_delta - w$recall.lrt)) })
  writeLines(c("\\begin{tabular}{rcccccccc}", "\\toprule",
    "& & & \\multicolumn{3}{c}{$\\hat\\delta/\\delta^*$} & \\multicolumn{3}{c}{Rappel} \\\\",
    "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
    "$K$ & $h^2_c$ & fond & prédit & départ & causaux fixés & oracle & RV & écart \\\\", "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_sel_delta.tex"))
}
cat("\nselection tables written\n")
