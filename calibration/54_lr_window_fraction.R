# ==============================================================================
# 54_lr_window_fraction.R
# Workload of step 3 of the procedure in the Discussion. Replays the LR paths of the GEUVADIS
# application (27c_geuvadis_application_cov.R, same data, kernels and eBIC-type level) and counts,
# at every step, the candidates whose nominal chi2_1 p-value of LR lies between level/10 and level,
# the range in which the protected critical value can change a decision.
# Usage (from the repository root): nice -n 10 Rscript calibration/54_lr_window_fraction.R [cores=4]
# Output: results/lr_window_fraction.csv (one row per population, gene and step)
# ==============================================================================
src <- readLines("calibration/27c_geuvadis_application_cov.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # data, step_stats, load_window, covariates
one <- function(pop, i) {
  g <- GENES[i, ]; W <- load_window(pop, g$tag); if (is.null(W)) return(NULL)
  X <- scale(W$G); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); C <- covariates(pop, W$samples)
  n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m); a_n <- cbind(1, C)
  lev <- pchisq(pen, 1, lower.tail = FALSE); lo <- qchisq(lev, 1, lower.tail = FALSE); hi <- qchisq(lev / 10, 1, lower.tail = FALSE)
  sel <- integer(0); out <- list()
  for (k in seq_len(K_MAX)) {
    rem <- setdiff(seq_len(m), sel); eg <- eigen(tcrossprod(X[, rem, drop = FALSE]), symmetric = TRUE)
    st <- step_stats(y, X[, rem, drop = FALSE], cbind(a_n, X[, sel, drop = FALSE]), eg$vectors, pmax(eg$values, 0), length(rem) - 1L)
    b <- which.max(st$lrt); ok <- st$lrt[b] > pen
    out[[k]] <- data.frame(pop = pop, gene = g$symbol, step = k, n = n, m = m, tested = length(rem), level = lev,
                           in_window = sum(st$lrt > lo & st$lrt <= hi), below_tenth = sum(st$lrt > hi),
                           max_lrt = st$lrt[b], accepted = ok)
    if (!ok) break
    sel <- c(sel, rem[b])
  }
  do.call(rbind, out)
}
jobs <- expand.grid(i = seq_len(nrow(GENES)), pop = names(POPS), stringsAsFactors = FALSE)
O <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(q) one(jobs$pop[q], jobs$i[q]), mc.cores = NCORES, mc.preschedule = FALSE))
write.csv(O, "calibration/results/lr_window_fraction.csv", row.names = FALSE)
options(width = 200); print(O, digits = 3)
for (p in c("EUR", "YRI", "ALL")) { z <- if (p == "ALL") O else O[O$pop == p, ]; f <- z[z$step == 1, ]
  cat(sprintf("%s: first steps %d tests, %d in window (%.2f%%), %d below level/10; all steps %d tests, %d in window (%.2f%%)\n", p,
              sum(f$tested), sum(f$in_window), 100 * sum(f$in_window) / sum(f$tested), sum(f$below_tenth),
              sum(z$tested), sum(z$in_window), 100 * sum(z$in_window) / sum(z$tested))) }
