# ==============================================================================
# 34b_mlmm_entry_p.R
# Entry p-values of the forward path of MLMM 0.1.1 on the 1,800 datasets of the
# selection study. 30_comparators_selection.R stored MLMM's forward sequence and its
# per-step heritability, but took the minimum p-value of each scan as the entry p-value,
# which after the first step can be a cofactor whose projected column is numerically null.
# Here the p-value of the variant that enters at step k is recomputed as MLMM computes it:
# Gower-normalized kernel K_norm, V = h2 K_norm + (1 - h2) I with MLMM's h2 for the model
# with the first k - 1 cofactors (EMMA REML), whitening by the inverse Cholesky factor,
# F-test of the variant given the cofactors with (1, n - 1 - ncol(cof)) degrees of freedom.
# Part 1 checks the recomputation against fresh MLMM runs (entering variant read from
# MLMM's own scan) on 12 datasets.
# Output: results/comparators/mlmm_entry.rds (list id -> data.frame(idx, p))
# Usage (from the repository root): Rscript calibration/34b_mlmm_entry_p.R
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
.libPaths(c(Sys.getenv("COMPARATORS_RLIB", "comparators/Rlib"), .libPaths()))
suppressPackageStartupMessages(library(mlmm))
ROOT <- "calibration/results/comparators"
GENES <- c("ERAP2", "RPS26", "SLFN5", "PEX6-region")
D <- setNames(lapply(GENES, function(s) readRDS(file.path(ROOT, "data", paste0(s, ".rds")))), GENES)
KN <- lapply(D, function(d) { n <- nrow(d$Z); K <- tcrossprod(d$Z) / ncol(d$Z)
  (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * K) * K })
entry_p <- function(gene, y, seq, h2s) {
  Z <- D[[gene]]$Z; n <- nrow(Z); K <- KN[[gene]]
  vapply(seq_along(seq), function(k) {
    h2 <- h2s[k]; M <- solve(chol(h2 * K + (1 - h2) * diag(n)))
    C <- cbind(1, Z[, seq[seq_len(k - 1)], drop = FALSE]); Ct <- crossprod(M, C); Yt <- crossprod(M, y)
    Q <- qr.Q(qr(Ct)); R0 <- Yt - Q %*% crossprod(Q, Yt)
    xt <- crossprod(M %*% (diag(n) - tcrossprod(Q)), Z[, seq[k]])
    rss1 <- sum(lsfit(xt, R0, intercept = FALSE)$residuals^2); rss0 <- sum(R0^2)
    df2 <- n - 1 - ncol(C); pf((rss0 - rss1) / rss1 * df2, 1, df2, lower.tail = FALSE) }, 0)
}
## ---- part 1: check against fresh MLMM runs ------------------------------------------------------------
chk <- list(); set.seed(5)
for (gene in GENES[1:2]) for (k in c(which(D[[gene]]$meta$Kc == 3)[c(1, 30, 60)], which(D[[gene]]$meta$Kc == 1)[c(5, 50, 90)])) {
  d <- D[[gene]]; n <- nrow(d$Z); m <- ncol(d$Z); y <- d$Y[, k]; smp <- sprintf("i%03d", seq_len(n))
  X <- d$Z; colnames(X) <- paste0("s", seq_len(m)); rownames(X) <- smp
  K <- tcrossprod(d$Z) / m; dimnames(K) <- list(smp, smp)
  sink(tempfile()); fit <- mlmm(setNames(y, smp), X, K, nbchunks = 2, maxsteps = 11); sink()
  fw <- fit$step_table[grepl("^fwd", fit$step_table$step_) & fit$step_table$ncof > 0, ]
  seq <- as.integer(sub("^\\+s", "", fw$cof)); h2s <- c(fit$step_table$h2[1], fw$h2)[seq_along(seq)]
  native <- vapply(seq_along(seq), function(q) { ps <- fit$pval_step[[q]]$out
    ps$pval[match(paste0("s", seq[q]), as.character(ps$SNP))] }, 0)
  ours <- entry_p(gene, y, seq, h2s)
  chk[[length(chk) + 1]] <- data.frame(gene = gene, k = k, step = seq_along(seq), idx = seq, native = native, ours = ours)
}
C <- do.call(rbind, chk)
cat("check against MLMM's own scans: max |log10 ratio| =", max(abs(log10(C$native / C$ours)), na.rm = TRUE),
    " (", sum(is.na(C$native)), "unmatched )\n")
print(head(C, 12), digits = 4)
## ---- part 2: all datasets -------------------------------------------------------------------------------
out <- list()
for (gene in GENES) { d <- D[[gene]]
  for (k in seq_len(ncol(d$Y))) { id <- d$meta$id[k]; f <- file.path(ROOT, "r_arms", paste0(id, ".rds"))
    if (!file.exists(f)) next; r <- readRDS(f); if (!is.list(r$mlmm)) next
    fw <- r$mlmm$fwd; if (!nrow(fw)) { out[[id]] <- data.frame(idx = integer(0), p = numeric(0)); next }
    h2s <- c(r$mlmm$h2_null, fw$h2)[seq_len(nrow(fw))]
    out[[id]] <- data.frame(idx = fw$idx, p = entry_p(gene, d$Y[, k], fw$idx, h2s)) } }
saveRDS(out, file.path(ROOT, "mlmm_entry.rds")); cat("datasets:", length(out), "\n")
