# ==============================================================================
# 19_make_conditional_tables.R
# Tables of the conditional-delta lever (18_conditional_delta.R).
#   tab_cond_frontier.tex  global null P(>= 1 selection) and mean recall, per threshold
#   tab_cond_cells.tex     recall / precision at the eBIC threshold, per scenario
#   tab_cond_delta.tex     delta estimates relative to delta*, per scenario
#   results/conditional_delta_summary.txt  paired differences, per-test sizes, checks
# Usage (from the repository root): Rscript calibration/19_make_conditional_tables.R
# ==============================================================================
RES <- Sys.getenv("COND_RES", "calibration/results")      # overrides for partial checks
N <- Sys.getenv("COND_NOTE", "calibration/note")
S <- read.csv(file.path(RES, "conditional_delta_sel.csv"))
M <- read.csv(file.path(RES, "conditional_delta_meta.csv"))
P <- read.csv(file.path(RES, "conditional_delta_pertest.csv"))
D14 <- read.csv("calibration/results/selection_demo_replicates.csv")
 out <- character(0); say <- function(...) out <<- c(out, sprintf(...))
f2 <- function(v) sprintf("%.2f", v)
procs <- c("oracle", "rv", "c1", "c2", "c3", "c4", "c5")
lab <- c(oracle = "Oracle, $\\delta^*$", rv = "RV, $\\delta$ réestimé",
         c1 = "C1, $\\tilde\\delta$ fixé, ensemble libéral", c2 = "C2, $\\tilde\\delta$ fixé, ensemble SuSiE",
         c3 = "C3, RV conditionnel, ensemble libéral", c4 = "C4, RV conditionnel, ensemble SuSiE",
         c5 = "C5, point fixe")
say("datasets: %d (null %d, non-null %d); rv = 14 at eBIC: %d/%d", nrow(M), sum(M$Kc == 0), sum(M$Kc > 0),
    sum(M$check_rv, na.rm = TRUE), sum(!is.na(M$check_rv)))
S$key <- paste(S$id); nul <- S[S$Kc == 0, ]; nn <- S[S$Kc > 0, ]

## ---- frontier: null P(>= 1 selection) and non-null recall, per threshold -----------------------
TH <- c("a1e-04" = "$10^{-4}$", "a1e-05" = "$10^{-5}$", ebic = "eBIC", "a1e-06" = "$10^{-6}$")
TH <- TH[names(TH) %in% unique(S$thr_name)]
rows <- vapply(procs, function(p) {
  cells <- vapply(names(TH), function(t) {
    z <- nul[nul$proc == p & nul$thr_name == t, ]; w <- nn[nn$proc == p & nn$thr_name == t, ]
    sprintf("$%.1f$ & $%s$", 100 * mean(z$K_hat > 0), f2(mean(w$recall))) }, "")
  sprintf("%s & %s \\\\", lab[p], paste(cells, collapse = " & ")) }, "")
writeLines(c(sprintf("\\begin{tabular}{l%s}", strrep("cc", length(TH))), "\\toprule",
  sprintf("& %s \\\\", paste(sprintf("\\multicolumn{2}{c}{%s}", TH), collapse = " & ")),
  paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_along(TH), 2 * seq_along(TH) + 1), collapse = ""),
  sprintf("Procédure & %s \\\\", paste(rep("nul (\\%) & rappel", length(TH)), collapse = " & ")), "\\midrule",
  rows[1], "\\midrule", rows[-1], "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_cond_frontier.tex"))

## ---- recall / precision at eBIC, per scenario -------------------------------------------------------
sc <- unique(nn[, c("h2bg", "Kc", "h2c")]); sc <- sc[order(sc$h2c, sc$Kc, sc$h2bg), ]
eb <- nn[nn$thr_name == "ebic", ]
sus <- D14[D14$method == "susie" & D14$Kc > 0, ]
short <- c(rv = "RV", c1 = "C1", c2 = "C2", c3 = "C3", c4 = "C4", c5 = "C5", oracle = "oracle", susie = "SuSiE")
rows <- apply(sc, 1, function(s) {
  pick <- function(d) d[d$Kc == s["Kc"] & abs(d$h2c - s["h2c"]) < 1e-9 & abs(d$h2bg - s["h2bg"]) < 1e-9, ]
  cells <- vapply(names(short), function(p) { w <- if (p == "susie") pick(sus) else pick(eb[eb$proc == p, ])
    sprintf("%s/%s", f2(mean(w$recall, na.rm = TRUE)), f2(mean(w$precision, na.rm = TRUE))) }, "")
  sprintf("$%d$ & $%.2f$ & $%.2f$ & %s \\\\", as.integer(s["Kc"]), s["h2c"], s["h2bg"], paste(cells, collapse = " & ")) })
writeLines(c(sprintf("\\begin{tabular}{rcc%s}", strrep("c", length(short))), "\\toprule",
  sprintf("$K$ & $h^2_c$ & fond & %s \\\\", paste(short, collapse = " & ")), "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_cond_cells.tex"))

## ---- delta estimates relative to delta* -----------------------------------------------------------
sc2 <- unique(M[, c("h2bg", "Kc", "h2c")]); sc2 <- sc2[order(sc2$Kc > 0, sc2$h2c, sc2$Kc, sc2$h2bg), ]
rows <- apply(sc2, 1, function(s) {
  d <- M[M$Kc == s["Kc"] & abs(M$h2c - s["h2c"]) < 1e-9 & abs(M$h2bg - s["h2bg"]) < 1e-9, ]
  q <- function(v) sprintf("$%.2f$ {\\scriptsize [%.2f, %.2f]}", median(v), quantile(v, 0.1), quantile(v, 0.9))
  sprintf("%s & $%.2f$ & $%.2f$ & %s & %s & %s & $%.1f$ & $%.1f$ \\\\",
          if (s["Kc"] == 0) "nul" else sprintf("$%d$", as.integer(s["Kc"])), s["h2c"], s["h2bg"],
          q(d$delta_null / d$delta_star), q(d$delta_lib / d$delta_star), q(d$delta_sus / d$delta_star),
          mean(d$n_lib), mean(d$n_sus)) })
writeLines(c("\\begin{tabular}{rcccccrr}", "\\toprule",
  "& & & \\multicolumn{3}{c}{$\\hat\\delta/\\delta^*$, médiane [déciles 1 et 9]} & \\multicolumn{2}{c}{taille moyenne} \\\\",
  "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}",
  "$K$ & $h^2_c$ & fond & modèle nul & ensemble libéral & ensemble SuSiE & libéral & SuSiE \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}"), file.path(N, "tab_cond_delta.tex"))

## ---- paired summaries -----------------------------------------------------------------------------
W <- reshape(eb[, c("id", "proc", "recall", "any_off", "precision")], idvar = "id", timevar = "proc", direction = "wide")
pd <- function(a, b, col = "recall") { x <- W[[paste0(col, ".", a)]] - W[[paste0(col, ".", b)]]
  sprintf("%+.3f (SE %.3f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))) }
gap <- mean(W$recall.oracle) - mean(W$recall.rv)
for (p in procs) say("eBIC %-6s recall %.3f  prec %.3f  any_off(non-null) %.4f  null P(>=1) %.3f  d_recall vs rv %s  share of oracle gap %.2f",
  p, mean(W[[paste0("recall.", p)]]), mean(W[[paste0("precision.", p)]], na.rm = TRUE), mean(W[[paste0("any_off.", p)]]),
  mean(nul$K_hat[nul$proc == p & nul$thr_name == "ebic"] > 0), pd(p, "rv"),
  (mean(W[[paste0("recall.", p)]]) - mean(W$recall.rv)) / gap)
for (t in unique(S$thr_name)) say("null P(>=1) at %-7s: %s", t, paste(sprintf("%s %.3f", procs,
  vapply(procs, function(p) mean(nul$K_hat[nul$proc == p & nul$thr_name == t] > 0), 0)), collapse = ", "))
for (bg in c(0.10, 0.25)) say("null P(>=1) at eBIC, bg %.2f: %s", bg, paste(sprintf("%s %.3f", procs,
  vapply(procs, function(p) mean(nul$K_hat[nul$proc == p & nul$thr_name == "ebic" & abs(nul$h2bg - bg) < 1e-9] > 0), 0)), collapse = ", "))
pt <- aggregate(cbind(n_cand, rej_1e2, rej_1e3, rej_1e4) ~ proc, P, sum)
for (i in seq_len(nrow(pt))) say("per-test step 1 (null) %-6s size/alpha at 1e-2, 1e-3, 1e-4: %.2f, %.2f, %.2f",
  pt$proc[i], pt$rej_1e2[i] / pt$n_cand[i] / 1e-2, pt$rej_1e3[i] / pt$n_cand[i] / 1e-3, pt$rej_1e4[i] / pt$n_cand[i] / 1e-4)
fi <- M[, grep("^fixiter_", names(M))]; say("fixed point iterations: max %d, share hitting the cap %.3f", max(as.matrix(fi)), mean(as.matrix(fi) >= 8))
writeLines(out, file.path(RES, "conditional_delta_summary.txt")); cat(out, sep = "\n")
cat("\nconditional tables written\n")
