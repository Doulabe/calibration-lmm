# ==============================================================================
# 34c_patch_mlmm_frontier.R
# Replaces the MLMM threshold rows of results/comparators/selection_long.csv by those
# built from the corrected entry p-values of 34b_mlmm_entry_p.R (checked against MLMM's
# own scans to 2e-13), then recomputes selection_summary.csv and selection_frontier.csv
# exactly as 34_selection_comparators_summary.R does.
# Usage (from the repository root): Rscript calibration/34c_patch_mlmm_frontier.R
# ==============================================================================
ROOT <- "calibration/results/comparators"
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
prefix_at <- function(pp, lev) {
  if (is.null(pp) || !nrow(pp)) return(integer(0))
  ok <- pp$p <= lev; if (!ok[1]) integer(0) else pp$idx[seq_len(if (all(ok)) nrow(pp) else which(!ok)[1] - 1)]
}
L <- read.csv(file.path(ROOT, "selection_long.csv"), stringsAsFactors = FALSE)
L <- L[!(L$method == "mlmm" & L$rule == "threshold"), ]
EN <- readRDS(file.path(ROOT, "mlmm_entry.rds"))
TH <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
add <- do.call(rbind, lapply(names(EN), function(id) { mm <- META[META$id == id, ]; X <- D[[mm$gene]]$Z
  do.call(rbind, lapply(TH, function(a) data.frame(id = id, gene = mm$gene, h2bg = mm$h2bg, Kc = mm$Kc, h2c = mm$h2c,
    b = mm$b, main = mm$main, method = "mlmm", rule = "threshold", thr = format(a),
    t(metrics(prefix_at(EN[[id]], a), D[[mm$gene]]$causal[[mm$k]], X))))) }))
L <- rbind(L, add[, names(L)])
write.csv(L, file.path(ROOT, "selection_long.csv"), row.names = FALSE)
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
front <- do.call(rbind, lapply(unique(summ$method), function(mt) {
  d <- summ[summ$method == mt & summ$n_null >= 190, ]; d <- d[order(d$null_any, d$recall), ]
  at <- function(z) { if (nrow(d) < 2 || z < min(d$null_any) || z > max(d$null_any)) return(NA_real_)
    approx(d$null_any, d$recall, xout = z, ties = max)$y }
  data.frame(method = mt, recall_at_1 = at(0.01), recall_at_5 = at(0.05), null_min = min(d$null_any),
             null_max = max(d$null_any), n_points = nrow(d))
}))
write.csv(front, file.path(ROOT, "selection_frontier.csv"), row.names = FALSE)
options(width = 200)
print(summ[summ$method %in% c("mlmm", "rv", "gemma_lrt"), c("method", "rule", "thr", "null_any", "recall", "precision", "any_off")],
      digits = 3, row.names = FALSE)
print(front, digits = 3, row.names = FALSE)
