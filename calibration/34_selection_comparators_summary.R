# ==============================================================================
# 34_selection_comparators_summary.R
# Common metrics for every procedure of the selection-level study (1,800 datasets on
# four GEUVADIS cis windows: 1,000 global nulls and 800 datasets with causal variants).
# Sources:
#   18_conditional_delta.R      rv (LRT path, delta re-estimated under each model), oracle,
#                               levers c1-c5, at per-test thresholds and eBIC; SuSiE lead sets
#   14_selection_demo.R         wald_h1, score_h0, mbf_fixed, mbf_hc, cbf_deployed, ols_ebic
#                               (1,000 datasets of the demonstration)
#   14b_selection_extra_nulls.R wald_h1, score_h0, ols_ebic on the 800 extra background-only datasets
#   30_comparators_selection.R  GEMMA paths (Wald, LRT, score; fixed regional kernel),
#                               MLMM (native extBIC, mBonf, forward path), QTLtools (3 levels)
#   30_*  (ARMS=bslmm)          BSLMM inclusion probabilities (1,000 main datasets)
#   30b_farmcpu_selection.R     FarmCPU p-values (LD clumping r2 > 0.5)
#   31_susieinf_selection.py    SuSiE-inf credible sets
# Metrics (as in 14_selection_demo.R): a causal is recovered when a selection tags it at
# r2 >= 0.8; a selection is off-target when its r2 with every causal is below 0.25;
# under the global null every selection is off-target.
# Output: results/comparators/selection_long.csv, selection_summary.csv, selection_frontier.csv
# Usage (from the repository root): Rscript calibration/34_selection_comparators_summary.R
# ==============================================================================
suppressPackageStartupMessages(library(jsonlite))
ROOT <- "calibration/results/comparators"; RES <- "calibration/results"
GENES <- c("ERAP2", "RPS26", "SLFN5", "PEX6-region")
D <- setNames(lapply(GENES, function(s) readRDS(file.path(ROOT, "data", paste0(s, ".rds")))), GENES)
META <- do.call(rbind, lapply(GENES, function(s) cbind(D[[s]]$meta, gene = s, k = seq_len(nrow(D[[s]]$meta)))))
META$main <- META$Kc > 0 | META$b <= 25
metrics <- function(sel, causal, X) {                      # identical to 14_selection_demo.R
  kc <- length(causal); ks <- length(sel)
  if (ks == 0) return(c(K_hat = 0, tp = 0, fp08 = 0, fp025 = 0, recall = if (kc) 0 else NA,
                        precision = NA, f1 = if (kc) 0 else NA, any_off = 0, exact = 0))
  if (kc == 0) return(c(K_hat = ks, tp = 0, fp08 = ks, fp025 = ks, recall = NA, precision = 0,
                        f1 = NA, any_off = 1, exact = 0))
  r2 <- cor(X[, sel, drop = FALSE], X[, causal, drop = FALSE])^2
  best <- apply(r2, 1, max); tagged <- apply(r2 >= 0.8, 2, any)
  tp <- sum(best >= 0.8); prec <- tp / ks; rec <- mean(tagged)
  c(K_hat = ks, tp = tp, fp08 = sum(best < 0.8), fp025 = sum(best < 0.25), recall = rec,
    precision = prec, f1 = if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0,
    any_off = as.numeric(any(best < 0.25)), exact = sum(sel %in% causal))
}
clump <- function(cand, score, X, r2max = 0.5) {           # greedy LD clumping by decreasing score
  sel <- integer(0)
  for (j in cand[order(-score)]) if (!length(sel) || max(cor(X[, j], X[, sel, drop = FALSE])^2) <= r2max) sel <- c(sel, j)
  sel
}
rows <- list(); add <- function(id, method, rule, thr, sel) {
  mm <- META[META$id == id, ]; X <- D[[mm$gene]]$Z
  rows[[length(rows) + 1]] <<- data.frame(id = id, gene = mm$gene, h2bg = mm$h2bg, Kc = mm$Kc, h2c = mm$h2c, b = mm$b,
    main = mm$main, method = method, rule = rule, thr = thr, t(metrics(as.integer(sel), D[[mm$gene]]$causal[[mm$k]], X)))
}

## ---- 1. our engine (18 and 14) ---------------------------------------------------------------------
S <- read.csv(file.path(RES, "conditional_delta_sel.csv"))
S <- S[S$proc %in% c("rv", "oracle", "c1", "c2", "c3", "c4", "c5"), ]
S$main <- S$Kc > 0 | S$b <= 25
eng <- data.frame(id = S$id, gene = S$gene, h2bg = S$h2bg, Kc = S$Kc, h2c = S$h2c, b = S$b, main = S$main,
                  method = S$proc, rule = ifelse(S$thr_name == "ebic", "native", "threshold"),
                  thr = S$thr_name, S[, c("K_hat", "tp", "fp08", "fp025", "recall", "precision", "f1", "any_off", "exact")])
D14 <- read.csv(file.path(RES, "selection_demo_replicates.csv"))
D14$id <- sprintf("%s_b%03d", D14$cell, D14$b)
D14 <- D14[D14$method %in% c("wald_h1", "score_h0", "mbf_fixed", "mbf_hc", "cbf_deployed", "ols_ebic"), ]
d14 <- data.frame(id = D14$id, gene = D14$gene, h2bg = D14$h2bg, Kc = D14$Kc, h2c = D14$h2c, b = D14$b, main = TRUE,
                  method = D14$method, rule = "native", thr = "ebic",
                  D14[, c("K_hat", "tp", "fp08", "fp025", "recall", "precision", "f1", "any_off", "exact")])
## the 800 extra background-only datasets (b > 25) for the Wald and null-estimate paths of the
## article and conditional regression (14b_selection_extra_nulls.R), so that they cover 1,000 datasets
for (f in list.files(file.path(RES, "selection_extra_nulls"), "rds$", full.names = TRUE)) {
  r <- readRDS(f); for (mn in c("wald_h1", "score_h0", "ols_ebic")) { s <- r$selections[[mn]]
    if (!(length(s) == 1 && is.na(s))) add(r$id, mn, "native", "ebic", s) } }
## SuSiE lead sets for all 1,800 datasets (stored by 18_conditional_delta.R)
for (f in list.files(file.path(RES, "conditional_delta"), "rds$", full.names = TRUE)) {
  r <- readRDS(f); add(r$id, "susie", "native", "cs95", r$S_sus) }

## ---- 2. published software (30, 30b, 31) -----------------------------------------------------------
lev_ebic_t <- function(n, m, nu) 2 * pt(-sqrt((nu - 1) * expm1((log(n) + 2 * log(m)) / n)), nu - 1)
TH <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
prefix_at <- function(pp, lev) {                          # forward path: stop at the first p above the level
  if (is.null(pp) || !nrow(pp)) return(integer(0))
  lev <- rep_len(lev, nrow(pp)); ok <- pp$p <= lev
  if (!ok[1]) integer(0) else pp$idx[seq_len(if (all(ok)) nrow(pp) else which(!ok)[1] - 1)]
}
fc_all <- setNames(lapply(GENES, function(s) { f <- file.path(ROOT, "farmcpu", paste0(s, ".rds"))
  if (file.exists(f)) { L <- readRDS(f); setNames(L, vapply(L, `[[`, "", "id")) } else NULL }), GENES)
si_all <- setNames(lapply(GENES, function(s) { f <- file.path(ROOT, "susieinf", paste0(s, "_susie.jsonl"))
  if (!file.exists(f)) return(NULL); L <- lapply(readLines(f), fromJSON); setNames(L, vapply(L, `[[`, "", "id")) }), GENES)
for (i in seq_len(nrow(META))) {
  mm <- META[i, ]; id <- mm$id; X <- D[[mm$gene]]$Z; n <- nrow(X); m <- ncol(X)
  f <- file.path(ROOT, "r_arms", paste0(id, ".rds"))
  if (file.exists(f)) {
    r <- readRDS(f)
    if (is.list(r$gemma)) {
      for (tn in c("wald", "lrt", "score")) {
        pp <- r$gemma[[tn]]
        lev <- if (tn == "lrt") pchisq(log(n) + 2 * log(m), 1, lower.tail = FALSE) else lev_ebic_t(n, m, n - seq_len(max(nrow(pp), 1)))
        add(id, paste0("gemma_", tn), "native", "ebic", prefix_at(pp, lev))
        add(id, paste0("gemma_", tn), "bonferroni", "0.05/m", prefix_at(pp, 0.05 / m))
        for (a in TH) add(id, paste0("gemma_", tn), "threshold", format(a), prefix_at(pp, a))
      }
    }
    if (is.list(r$mlmm)) {
      add(id, "mlmm", "native", "extBIC", r$mlmm$extBIC); add(id, "mlmm_mbonf", "native", "mBonf", r$mlmm$mbonf)
      for (a in TH) add(id, "mlmm", "threshold", format(a), prefix_at(r$mlmm$fwd, a))
    }
    if (is.list(r$qtltools)) for (a in c("0.01", "0.05", "0.20"))
      add(id, "qtltools", if (a == "0.05") "native" else "threshold", a, r$qtltools[[paste0("a", a)]])
  }
  fb <- file.path(ROOT, "r_arms_bslmm", paste0(id, ".rds"))
  if (file.exists(fb)) { bs <- readRDS(fb)$bslmm
    if (is.list(bs)) for (th in c(0.1, 0.25, 0.5, 0.9)) {
      cand <- bs$idx[bs$pip >= th]; sel <- if (length(cand)) clump(cand, bs$pip[bs$pip >= th], X) else integer(0)
      add(id, "bslmm", if (th == 0.5) "native" else "threshold", format(th), sel) } }
  fc <- fc_all[[mm$gene]][[id]]
  if (!is.null(fc) && is.null(fc$error)) for (a in c(TH, 0.05 / m)) {
    cand <- fc$idx[fc$p <= a]; sel <- if (length(cand)) clump(cand, -log(fc$p[fc$p <= a]), X) else integer(0)
    add(id, "farmcpu", if (a == 0.05 / m) "native" else "threshold", if (a == 0.05 / m) "0.05/m" else format(a), sel) }
  si <- si_all[[mm$gene]][[id]]
  if (!is.null(si)) {
    top <- unlist(si$top); pipv <- if (length(top)) vapply(top, function(t) { w <- which(unlist(si$pip_idx) == t)
      if (length(w)) unlist(si$pip)[w] else 0 }, 0) else numeric(0)
    for (th in c(0, 0.5, 0.9, 0.99)) add(id, "susieinf", if (th == 0) "native" else "threshold", format(th), top[pipv >= th])
  }
}
L <- rbind(eng, d14, do.call(rbind, rows))
write.csv(L, file.path(ROOT, "selection_long.csv"), row.names = FALSE)

## ---- 3. summaries ---------------------------------------------------------------------------------------
key <- unique(L[, c("method", "rule", "thr")])
summ <- do.call(rbind, lapply(seq_len(nrow(key)), function(q) {
  d <- L[L$method == key$method[q] & L$rule == key$rule[q] & L$thr == key$thr[q], ]
  nul <- d[d$Kc == 0, ]; nn <- d[d$Kc > 0, ]
  data.frame(key[q, ], n_null = nrow(nul), null_any = mean(nul$K_hat > 0), null_K = mean(nul$K_hat),
             n_causal = nrow(nn), recall = mean(nn$recall), precision = mean(nn$precision, na.rm = TRUE),
             K_hat = mean(nn$K_hat), any_off = mean(nn$any_off), f1 = mean(nn$f1),
             null_any_main = mean(nul$K_hat[nul$main] > 0), n_null_main = sum(nul$main))
}))
write.csv(summ, file.path(ROOT, "selection_summary.csv"), row.names = FALSE)
## recall at a common rate of datasets with a selection under the global null (linear interpolation)
front <- do.call(rbind, lapply(unique(summ$method), function(mt) {
  d <- summ[summ$method == mt & summ$n_null >= 200, ]; d <- d[order(d$null_any, d$recall), ]
  at <- function(z) { if (nrow(d) < 2 || z < min(d$null_any) || z > max(d$null_any)) return(NA_real_)
    approx(d$null_any, d$recall, xout = z, ties = max)$y }
  data.frame(method = mt, recall_at_1 = at(0.01), recall_at_5 = at(0.05), null_min = min(d$null_any), null_max = max(d$null_any),
             n_points = nrow(d))
}))
write.csv(front, file.path(ROOT, "selection_frontier.csv"), row.names = FALSE)
options(width = 220)
print(summ[summ$rule %in% c("native", "bonferroni"), ], digits = 3, row.names = FALSE)
print(front, digits = 3, row.names = FALSE)
