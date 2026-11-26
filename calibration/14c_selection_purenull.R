# ==============================================================================
# 14c_selection_purenull.R
# Pure-noise control of the selection study (y = e, 29b_export_purenull.R): the stepwise
# paths of the article on the candidate-specific kernel (Wald at delta_hat, test at the null
# estimate, likelihood ratio; eBIC level), conditional regression (eBIC) and SuSiE (susieR,
# L = 10, top PIP per 95% credible set), on the 500 datasets.
# Usage (from the repository root): nice -n 10 Rscript calibration/14c_selection_purenull.R [cores=2]
# Output: results/selection_purenull/<id>.rds, results/selection_purenull.csv
# ==============================================================================
src <- readLines("calibration/14b_selection_extra_nulls.R")
eval(parse(text = src[seq_len(grep("^GEN <- ", src)[1] - 1L)]))    # settings, algebra, paths (no Bayes factors)
GEN <- c("ERAP2", "RPS26", "SLFN5", "PEX6-region")
D <- setNames(lapply(GEN, function(s) readRDS(file.path("calibration/results/comparators_purenull/data", paste0(s, ".rds")))), GEN)
JOBS <- do.call(rbind, lapply(GEN, function(s) data.frame(gene = s, k = seq_len(nrow(D[[s]]$meta)), id = D[[s]]$meta$id,
                                                          b = D[[s]]$meta$b, stringsAsFactors = FALSE)))
OUTD <- "calibration/results/selection_purenull"; dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
one <- function(j) {
  g <- JOBS$gene[j]; id <- JOBS$id[j]; ck <- file.path(OUTD, paste0(id, ".rds"))
  if (file.exists(ck) && file.size(ck) > 0) return("cached")
  Z <- D[[g]]$Z; y <- D[[g]]$Y[, JOBS$k[j]]
  t0 <- Sys.time(); sels <- run_treatments(y, Z)
  sels$ols_ebic <- tryCatch(stepwise_baseline(y, Z, kernel = "none", delta_mode = "once", stop = "ebic", K_max = K_MAX)$indices,
                            error = function(e) NA_integer_)
  sels$susie <- tryCatch({ f <- susie(Z, y, L = 10, verbose = FALSE); cs <- f$sets$cs
    if (is.null(cs)) integer(0) else unique(vapply(cs, function(ix) ix[which.max(f$pip[ix])], 1L)) }, error = function(e) NA_integer_)
  saveRDS(list(id = id, gene = g, selections = sels, secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), ck)
  "done"
}
suppressPackageStartupMessages(library(parallel))
st <- unlist(mclapply(order(JOBS$b, JOBS$gene), function(j) tryCatch(one(j), error = function(e) paste("FAIL", JOBS$id[j], conditionMessage(e))),
                      mc.cores = NCORES, mc.preschedule = FALSE))
cat("datasets:", sum(st == "done"), "new,", sum(st == "cached"), "cached,", sum(grepl("FAIL", st)), "failed\n"); print(head(st[grepl("FAIL", st)]))
R <- do.call(rbind, lapply(list.files(OUTD, "rds$", full.names = TRUE), function(f) { r <- readRDS(f)
  do.call(rbind, lapply(names(r$selections), function(mn) { s <- r$selections[[mn]]
    data.frame(id = r$id, gene = r$gene, method = mn, K_hat = if (length(s) == 1 && is.na(s)) NA else length(s)) })) }))
write.csv(R, "calibration/results/selection_purenull.csv", row.names = FALSE)
print(aggregate(cbind(any = K_hat > 0) ~ method, R, mean))
