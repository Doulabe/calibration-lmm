# ==============================================================================
# 44_selection_bootstrap.R
# Uncertainty of the recalls at equal error of the selection study.
# Stratified bootstrap of the datasets, resampled with replacement within each design cell
# (window and h2_bg, and for the datasets with causal variants their number and size; the
# background-only datasets of the main design and the 800 extra ones are separate strata). For
# each replicate and procedure, the rate of background-only datasets with a selection and the
# mean recall are recomputed at every threshold, and the recall at a rate of 1 (and 5) per cent
# is interpolated exactly as in 34_selection_comparators_summary.R (points with at least 200
# background-only datasets, linear interpolation, NA outside the range of the rates).
# Input: results/comparators/selection_long.csv (after 34c_patch_mlmm_frontier.R).
# Output: results/selection_bootstrap.csv (estimate and percentile interval of the recall at 1
# and 5 per cent per procedure) and results/selection_bootstrap_diff.csv (differences between
# the paths compared in the text).
# Usage (from the repository root): Rscript calibration/44_selection_bootstrap.R [R=2000]
# ==============================================================================
a0 <- commandArgs(TRUE)
NB <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 2000L
L <- read.csv("calibration/results/comparators/selection_long.csv", stringsAsFactors = FALSE)
L$h2c[is.na(L$h2c)] <- 0
## datasets and strata
ids <- unique(L[, c("id", "gene", "h2bg", "Kc", "h2c", "main")])
stopifnot(!anyDuplicated(ids$id))
ids$stratum <- with(ids, ifelse(Kc == 0, sprintf("%s_%g_null_%s", gene, h2bg, ifelse(main, "main", "extra")),
                                sprintf("%s_%g_K%d_%g", gene, h2bg, Kc, h2c)))
cat("datasets:", nrow(ids), " strata:", length(unique(ids$stratum)), "\n")
## one row per (method, rule, threshold): selection indicator on background-only datasets, recall otherwise
L$grp <- paste(L$method, L$rule, L$thr, sep = "|")
G <- unique(L$grp); ix <- match(L$id, ids$id); ig <- match(L$grp, G)
nul <- L$Kc == 0
Msel <- Mnul <- Mrec <- Mcau <- matrix(0, length(G), nrow(ids))
Msel[cbind(ig[nul], ix[nul])] <- as.numeric(L$K_hat[nul] > 0); Mnul[cbind(ig[nul], ix[nul])] <- 1
Mrec[cbind(ig[!nul], ix[!nul])] <- L$recall[!nul]; Mcau[cbind(ig[!nul], ix[!nul])] <- 1
gm <- do.call(rbind, strsplit(G, "|", fixed = TRUE)); n_null <- rowSums(Mnul)
## frontier of 34: recall at a common rate of background-only datasets with a selection
frontier <- function(null_any, recall, z) {
  vapply(unique(gm[, 1]), function(mt) {
    k <- which(gm[, 1] == mt & n_null >= 200)
    if (length(k) < 2) return(rep(NA_real_, length(z)))
    o <- order(null_any[k], recall[k]); x <- null_any[k][o]; y <- recall[k][o]
    vapply(z, function(zz) if (zz < min(x) || zz > max(x)) NA_real_ else approx(x, y, xout = zz, ties = max)$y, 0)
  }, numeric(length(z)))
}
est <- frontier(as.numeric(Msel %*% rep(1, nrow(ids)) / n_null),
                as.numeric(Mrec %*% rep(1, nrow(ids)) / rowSums(Mcau)), c(0.01, 0.05))
## bootstrap weights: multiplicity of each dataset, resampled within strata
set.seed(44000)
st <- split(seq_len(nrow(ids)), ids$stratum)
W <- matrix(0, nrow(ids), NB)
for (b in seq_len(NB)) for (s in st) { k <- s[sample.int(length(s), length(s), replace = TRUE)]
  W[, b] <- W[, b] + tabulate(k, nbins = nrow(ids)) }
NA_ <- (Msel %*% W) / (Mnul %*% W); RC <- (Mrec %*% W) / (Mcau %*% W)
BS <- array(NA_real_, c(2, ncol(est), NB), dimnames = list(c("at1", "at5"), colnames(est), NULL))
for (b in seq_len(NB)) BS[, , b] <- frontier(NA_[, b], RC[, b], c(0.01, 0.05))
q <- function(v) { v <- v[is.finite(v)]; if (length(v) < 0.9 * NB) c(NA, NA) else unname(quantile(v, c(0.025, 0.975))) }
OUT <- do.call(rbind, lapply(colnames(est), function(mt) data.frame(method = mt,
  recall_at_1 = est[1, mt], lo_1 = q(BS[1, mt, ])[1], hi_1 = q(BS[1, mt, ])[2], na_share_1 = mean(!is.finite(BS[1, mt, ])),
  recall_at_5 = est[2, mt], lo_5 = q(BS[2, mt, ])[1], hi_5 = q(BS[2, mt, ])[2])))
write.csv(OUT, "calibration/results/selection_bootstrap.csv", row.names = FALSE)
## differences between the paths compared in the text, on the same bootstrap replicates
PAIRS <- list(c("rv", "gemma_lrt"), c("rv", "mlmm"), c("rv", "gemma_score"), c("gemma_lrt", "gemma_score"),
              c("gemma_lrt", "mlmm"), c("rv", "gemma_wald"), c("gemma_lrt", "gemma_wald"), c("mlmm", "gemma_wald"),
              c("gemma_score", "gemma_wald"), c("oracle", "rv"), c("oracle", "gemma_lrt"), c("oracle", "mlmm"),
              c("oracle", "gemma_score"))
PAIRS <- Filter(function(p) all(p %in% colnames(est)), PAIRS)
DIF <- do.call(rbind, lapply(PAIRS, function(p) { d <- BS[1, p[1], ] - BS[1, p[2], ]
  data.frame(a = p[1], b = p[2], diff_at_1 = est[1, p[1]] - est[1, p[2]], lo = q(d)[1], hi = q(d)[2],
             share_positive = mean(d[is.finite(d)] > 0)) }))
write.csv(DIF, "calibration/results/selection_bootstrap_diff.csv", row.names = FALSE)
options(width = 200)
print(OUT[order(-OUT$recall_at_1), ], digits = 3, row.names = FALSE)
print(DIF, digits = 3, row.names = FALSE)
