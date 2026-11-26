# ==============================================================================
# 14b_selection_extra_nulls.R
# Completes the background-only arm of the selection study: the stepwise paths of the
# article on the candidate-specific kernel (Wald at delta_hat, test at the null estimate,
# likelihood ratio; 14_selection_demo.R) and conditional regression without random effect
# (stepwise_baseline, eBIC) on the 800 background-only datasets b = 26..125 exported by
# 29_export_selection_datasets.R, so that every procedure is evaluated on the same
# 1,000 background-only datasets. The marginal Bayes factor paths of 14 are not run.
# Usage (from the repository root): nice -n 10 Rscript calibration/14b_selection_extra_nulls.R [cores=2]
# Output: results/selection_extra_nulls/<id>.rds, results/selection_extra_nulls.csv
# ==============================================================================
src <- readLines("calibration/14_selection_demo.R")
stop_at <- grep("^## ---- data", src)[1] - 1L                # settings, engine algebra, paths, metrics
src <- sub("^NCORES <- .*$", "NCORES <- 2L", src)             # the core count is set below
eval(parse(text = src[seq_len(stop_at)]))
args <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L

## step statistics without the marginal Bayes factors (identical algebra to step_stats() of 14)
step_stats <- function(y, Xc, F_k, U, s_raw, m_k) {
  n <- length(y); k <- ncol(F_k); nu <- n - k; mc <- ncol(Xc); G <- length(GRIDD)
  s_rank <- pmax(s_raw, 0) / m_k
  Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F_k); UtX <- crossprod(U, Xc)
  Dm <- um <- R0 <- lZ <- lM <- matrix(0, mc, G)
  for (g in seq_len(G)) {
    pc <- pieces_x(Uty, UtF, UtX, s_rank, GRIDD[g], GRIDD[g] / m_k, n, k)
    Dm[, g] <- pc$D; um[, g] <- pc$u; R0[, g] <- pc$RSS0; lZ[, g] <- pc$ldZ; lM[, g] <- pc$ldM
  }
  R1 <- pmax(R0 - um^2 / Dm, 1e-300)
  l0 <- -0.5 * (lZ + lM) - 0.5 * nu * log(R0)
  l1 <- -0.5 * (lZ + lM + log(Dm)) - 0.5 * (nu - 1) * log(R1)
  g1 <- max.col(l1, ties.method = "first"); g0 <- max.col(l0, ties.method = "first")
  ix1 <- cbind(seq_len(mc), g1); ix0 <- cbind(seq_len(mc), g0)
  t2w <- (nu - 1) * (um[ix1]^2 / Dm[ix1]) / R1[ix1]
  Q2 <- nu * um[ix0]^2 / (Dm[ix0] * R0[ix0]); t2s <- (nu - 1) * Q2 / pmax(nu - Q2, 1e-12)
  lrt <- apply(-lZ - n * log(R1), 1, max) - apply(-lZ - n * log(R0), 1, max)
  list(t2w = t2w, t2s = t2s, lrt = lrt, nu = nu)
}
TREAT <- c("wald_h1", "score_h0", "lrt")

GEN <- c("ERAP2", "RPS26", "SLFN5", "PEX6-region")
D <- setNames(lapply(GEN, function(s) readRDS(file.path("calibration/results/comparators/data", paste0(s, ".rds")))), GEN)
JOBS <- do.call(rbind, lapply(GEN, function(s) { m <- D[[s]]$meta
  data.frame(gene = s, k = seq_len(nrow(m)), id = m$id, Kc = m$Kc, b = m$b, stringsAsFactors = FALSE) }))
JOBS <- JOBS[JOBS$Kc == 0 & JOBS$b > 25, ]
OUTD <- "calibration/results/selection_extra_nulls"; dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
one <- function(j) {
  g <- JOBS$gene[j]; id <- JOBS$id[j]; ck <- file.path(OUTD, paste0(id, ".rds"))
  if (file.exists(ck) && file.size(ck) > 0) return("cached")
  Z <- D[[g]]$Z; y <- D[[g]]$Y[, JOBS$k[j]]
  t0 <- Sys.time(); sels <- run_treatments(y, Z)
  ols <- tryCatch(stepwise_baseline(y, Z, kernel = "none", delta_mode = "once", stop = "ebic", K_max = K_MAX)$indices,
                  error = function(e) NA_integer_)
  sels$ols_ebic <- ols
  saveRDS(list(id = id, gene = g, selections = sels, secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), ck)
  "done"
}
suppressPackageStartupMessages(library(parallel))
ord <- order(JOBS$b, JOBS$gene)
st <- unlist(mclapply(ord, function(j) tryCatch(one(j), error = function(e) paste("FAIL", JOBS$id[j], conditionMessage(e))),
                      mc.cores = NCORES, mc.preschedule = FALSE))
cat("datasets:", sum(st == "done"), "new,", sum(st == "cached"), "cached,", sum(grepl("FAIL", st)), "failed\n"); print(head(st[grepl("FAIL", st)]))
R <- do.call(rbind, lapply(list.files(OUTD, "rds$", full.names = TRUE), function(f) { r <- readRDS(f)
  do.call(rbind, lapply(names(r$selections), function(mn) { s <- r$selections[[mn]]
    data.frame(id = r$id, gene = r$gene, method = mn, K_hat = if (length(s) == 1 && is.na(s)) NA else length(s)) })) }))
write.csv(R, "calibration/results/selection_extra_nulls.csv", row.names = FALSE)
print(aggregate(cbind(any = K_hat > 0) ~ method, R, mean))
