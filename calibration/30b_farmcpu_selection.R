# ==============================================================================
# 30b_farmcpu_selection.R
# FarmCPU (rMVP 1.4.6, MVP.FarmCPU) on the 1,800 datasets of the selection-level study,
# with the settings of 30_comparators_selection.R (static bins of 10, 50 and 100 kb,
# bin.selection 5..50, maxLoop 10, p.threshold 0.05/m). Unlike 30_*, it also records the
# markers whose final p-value is reported as NaN by rMVP. On real genotypes these are,
# apart from rare exceptions recorded here, exact duplicates (r2 = 1) of pseudo-QTNs,
# which carry their own p-value.
# Output: results/comparators/farmcpu/<gene>.rds (list over datasets: idx, p for p <= 1e-3,
#         NaN markers and their maximal r2 with a finite-p marker at p <= 1e-3)
# Usage (from the repository root): nice -n 10 Rscript calibration/30b_farmcpu_selection.R
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
suppressPackageStartupMessages({ library(bigmemory); library(rMVP) })
ROOT <- "calibration/results/comparators"; OUT <- file.path(ROOT, "farmcpu")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GENES <- data.frame(symbol = c("ERAP2", "RPS26", "SLFN5", "PEX6-region"), chr = c(5L, 12L, 17L, 6L))
for (gi in seq_len(nrow(GENES))) {
  s <- GENES$symbol[gi]; fo <- file.path(OUT, paste0(s, ".rds")); if (file.exists(fo)) next
  D <- readRDS(file.path(ROOT, "data", paste0(s, ".rds"))); m <- ncol(D$G); n <- nrow(D$G)
  geno <- bigmemory::as.big.matrix(t(D$G)); map <- data.frame(SNP = paste0("s", seq_len(m)), CHROM = GENES$chr[gi], POS = D$pos)
  res <- vector("list", ncol(D$Y)); t0 <- Sys.time()
  for (k in seq_len(ncol(D$Y))) {
    r <- tryCatch(MVP.FarmCPU(phe = data.frame(id = sprintf("i%03d", seq_len(n)), y = D$Y[, k]), geno = geno, map = map,
                              CV = NULL, maxLoop = 10, method.bin = "static", bin.size = c(1e4, 5e4, 1e5),
                              bin.selection = seq(5, 50, 5), p.threshold = 0.05 / m, ncpus = 1, verbose = FALSE),
                  error = function(e) NULL)
    if (is.null(r)) { res[[k]] <- list(id = D$meta$id[k], error = TRUE); next }
    p <- r[, ncol(r)]; keep <- which(is.finite(p) & p <= 1e-3); nan <- which(!is.finite(p))
    dup <- if (length(nan) && length(keep)) apply(cor(D$Z[, nan, drop = FALSE], D$Z[, keep, drop = FALSE])^2, 1, max) else rep(NA_real_, length(nan))
    res[[k]] <- list(id = D$meta$id[k], idx = keep, p = p[keep], nan = nan, nan_r2 = dup)
  }
  saveRDS(res, fo)
  cat(s, ncol(D$Y), "datasets", format(Sys.time() - t0), "\n")
}
