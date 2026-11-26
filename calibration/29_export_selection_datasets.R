# ==============================================================================
# 29_export_selection_datasets.R
# Regenerates the 1,800 datasets of the selection-level study (14_selection_demo.R,
# 18_conditional_delta.R: four GEUVADIS cis windows, n = 358, 1,000 datasets with
# K in {0,1,3} plus 800 extra global-null datasets) from the same seeds, checks the
# causal sets against the stored 14_selection_demo.R outputs, and exports them for
# the external software comparators (30_comparators_selection.R, 31_susieinf_selection.py).
# Output: results/comparators/data/<gene>.rds (Z, positions, Y, meta, causal sets)
#         results/comparators/data/<gene>_Z.txt, <gene>_Y.txt, <gene>_meta.csv (for Python)
# Usage (from the repository root): Rscript calibration/29_export_selection_datasets.R
# ==============================================================================
B <- 25L; B_NULL <- 125L
GENES <- data.frame(tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000124587"),
                    symbol = c("ERAP2", "RPS26", "SLFN5", "PEX6-region"), chr = c(5L, 12L, 17L, 6L),
                    stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_window <- function(tag) {                     # as load_Z() of 14_selection_demo.R, keeping positions
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(tag, ".gt"))), "\t", fixed = TRUE)
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  pos <- as.integer(vapply(L, `[`, "", 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  Z <- scale(G[, keep, drop = FALSE]); Z[!is.finite(Z)] <- 0
  list(Z = Z, G = G[, keep, drop = FALSE], pos = pos[keep])
}
grid <- rbind(
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = c(1L, 3L), h2c = c(0.08, 0.15),
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = (B + 1L):B_NULL, stringsAsFactors = FALSE))
grid$cell <- sprintf("%s_bg%.2f_K%d_c%.2f", GENES$symbol[grid$gene], grid$h2bg, grid$Kc, grid$h2c)
grid$id <- sprintf("%s_b%03d", grid$cell, grid$b)
OUT <- "calibration/results/comparators/data"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
n_checked <- 0L
for (gi in seq_len(nrow(GENES))) {
  W <- load_window(GENES$tag[gi]); Z <- W$Z; n <- nrow(Z); m <- ncol(Z)
  rows <- which(grid$gene == gi)
  Y <- matrix(NA_real_, n, length(rows)); causal <- vector("list", length(rows))
  for (k in seq_along(rows)) {
    g <- grid[rows[k], ]
    set.seed(30000 + 1000 * g$gene + 100 * g$Kc + round(100 * g$h2bg) + 7919 * g$b + round(1000 * g$h2c))
    cs <- integer(0)
    while (length(cs) < g$Kc) {
      cand <- sample.int(m, 1)
      if (!length(cs) || max(cor(Z[, cand], Z[, cs, drop = FALSE])^2) < 0.5) cs <- c(cs, cand)
    }
    bvec <- rnorm(m); u <- as.numeric(Z %*% bvec); u <- u / sd(u) * sqrt(g$h2bg)
    sig <- as.numeric(Z[, cs, drop = FALSE] %*% rep(sqrt(g$h2c), g$Kc))
    s2e <- max(1 - g$h2bg - g$Kc * g$h2c, 0.05)
    Y[, k] <- sig + u + rnorm(n) * sqrt(s2e); causal[[k]] <- cs
    f14 <- file.path("calibration/results/selection_demo", paste0(g$id, ".rds"))
    if (file.exists(f14)) {
      st <- readRDS(f14)
      if (!identical(as.integer(st$causal), as.integer(cs))) stop("causal set differs from 14: ", g$id)
      n_checked <- n_checked + 1L
    }
  }
  meta <- grid[rows, c("id", "cell", "h2bg", "Kc", "h2c", "b")]
  meta$delta_star <- meta$h2bg / pmax(1 - meta$h2bg - meta$Kc * meta$h2c, 0.05)
  meta$causal <- vapply(causal, function(v) paste(v, collapse = ";"), "")
  sym <- GENES$symbol[gi]
  saveRDS(list(gene = sym, chr = GENES$chr[gi], Z = Z, G = W$G, pos = W$pos, Y = Y, meta = meta, causal = causal),
          file.path(OUT, paste0(sym, ".rds")))
  write.table(Z, file.path(OUT, paste0(sym, "_Z.txt")), row.names = FALSE, col.names = FALSE)
  write.table(Y, file.path(OUT, paste0(sym, "_Y.txt")), row.names = FALSE, col.names = FALSE)
  write.csv(meta, file.path(OUT, paste0(sym, "_meta.csv")), row.names = FALSE)
  cat(sprintf("%s: n = %d, m = %d, %d datasets\n", sym, n, m, length(rows)))
}
cat("causal sets checked against 14_selection_demo.R:", n_checked, "\n")
