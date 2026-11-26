# ==============================================================================
# 29b_export_purenull.R
# Datasets of the pure-noise control of the selection study: on the same four GEUVADIS
# windows (n = 358, real genotypes), y = e with e ~ N(0, I), no regional background and no
# causal variant, 125 datasets per window. Same format as 29_export_selection_datasets.R,
# for 30_comparators_selection.R and 31_susieinf_selection.py (COMP_ROOT=comparators_purenull)
# and 14c_selection_purenull.R.
# Output: results/comparators_purenull/data/<gene>.rds, <gene>_Z.txt, <gene>_Y.txt, <gene>_meta.csv
# Usage (from the repository root): Rscript calibration/29b_export_purenull.R
# ==============================================================================
src <- readLines("calibration/29_export_selection_datasets.R")
eval(parse(text = src[seq_len(grep("^grid <- rbind", src)[1] - 1L)]))   # GENES, load_window()
B_PURE <- 125L
OUT <- "calibration/results/comparators_purenull/data"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
for (gi in seq_len(nrow(GENES))) {
  W <- load_window(GENES$tag[gi]); Z <- W$Z; n <- nrow(Z); sym <- GENES$symbol[gi]
  Y <- vapply(seq_len(B_PURE), function(b) { set.seed(50000 + 1000 * gi + b); rnorm(n) }, numeric(n))
  meta <- data.frame(id = sprintf("%s_pure_b%03d", sym, seq_len(B_PURE)), cell = paste0(sym, "_pure"), h2bg = 0, Kc = 0L,
                     h2c = 0, b = seq_len(B_PURE), delta_star = 0, causal = "", stringsAsFactors = FALSE)
  saveRDS(list(gene = sym, chr = GENES$chr[gi], Z = Z, G = W$G, pos = W$pos, Y = Y, meta = meta,
               causal = rep(list(integer(0)), B_PURE)), file.path(OUT, paste0(sym, ".rds")))
  write.table(Z, file.path(OUT, paste0(sym, "_Z.txt")), row.names = FALSE, col.names = FALSE)
  write.table(Y, file.path(OUT, paste0(sym, "_Y.txt")), row.names = FALSE, col.names = FALSE)
  write.csv(meta, file.path(OUT, paste0(sym, "_meta.csv")), row.names = FALSE)
  cat(sprintf("%s: n = %d, m = %d, %d datasets\n", sym, n, ncol(Z), B_PURE))
}
