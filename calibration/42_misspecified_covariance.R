# ==============================================================================
# 42_misspecified_covariance.R
# Exact sizes when the true covariance is not the working one.
# The working model is that of the calibration study
# (36_calibration_enlarged.R): intercept, regional kernel K of the window without the
# candidate, a single ratio delta. The truth is Gaussian with covariance Sigma, and the
# decomposition behind the exact distributions still holds: with b the generalized
# least-squares estimate of beta under the TRUE covariance and r = s - xi b, b and r are
# independent and the working statistics depend on s = r + xi b, the Wald statistics
# affinely in b given r. Only the draw of r changes, the conditional probabilities being
# computed as in 36 with S0 = xi' Sigma^-1 xi.
# Scenarios (background heritability 0.20 in total, kernels scaled to unit mean diagonal):
#   correct  Sigma = I + delta0 K, delta0 = 0.25 (the calibration study)
#   sparse   Sigma = I + 0.25 K_S, K_S from 20 markers of the window drawn at random
#            (seed per window, the candidate excluded): a sparse regional background
#   outside  Sigma = I + 0.125 K + 0.125 K_out, K_out from the markers of the seven other
#            GEUVADIS windows (other chromosomes): half of the background outside the window
# Configurations: the eight GEUVADIS windows (European sample, n = 358), three candidates
# per window at the quantiles 0.10, 0.50 and 0.90 of nu_j (h2 = 0.20), 1,000 draws.
# Treatments: Wald at delta_hat (REML), t at the null estimate (REML, as EMMAX), likelihood
# ratio (ML), Kenward-Roger, and the Wald test at the pseudo-true ratio delta* of the working
# null model (maximizer of the expected restricted log-likelihood), whose size measures the
# misspecification alone; each treatment is referred to its usual reference distribution.
# Usage (from the repository root): nice -n 10 Rscript calibration/42_misspecified_covariance.R [cores=3] [--smoke]
# Output: results/misspecified/<config>_<scenario>.rds, results/misspecified.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # GENES, load_geu, NUS, pick_grid, crt_factors, GRID2, AL, NBIS
a0 <- commandArgs(TRUE)
NC <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 3L
SMOKE42 <- "--smoke" %in% a0
B42 <- if (SMOKE42) 200L else 1000L
OUTM <- paste0("calibration/results/misspecified", if (SMOKE42) "_smoke" else "")
dir.create(OUTM, showWarnings = FALSE, recursive = TRUE)
QM <- c(0.10, 0.50, 0.90); H2 <- 0.20; NSP <- 20L

## ---- candidates: three per window at quantiles of nu_j (h2 = 0.20), as in 36 --------------------------
CF <- do.call(rbind, lapply(seq_len(nrow(GENES)), function(i) {
  d <- NUS[NUS$gene == GENES$symbol[i] & abs(NUS$h2 - H2) < 1e-9, ]; lv <- log(d$nu_delta); picks <- integer(0)
  for (q in quantile(lv, QM)) { o <- order(abs(lv - q)); o <- o[!d$pos[o] %in% d$pos[picks]]; picks <- c(picks, o[1]) }
  data.frame(i = i, gene = GENES$symbol[i], pos = d$pos[picks], q = QM) }))
SCEN <- c("correct", "sparse", "outside")
JOBS <- merge(CF, data.frame(scen = SCEN), by = NULL)
JOBS$id <- sprintf("%s_%d_%s", gsub("[^A-Za-z0-9]", "", JOBS$gene), JOBS$pos, JOBS$scen)
if (SMOKE42) JOBS <- JOBS[JOBS$i == 1 & JOBS$q == 0.5, ]

mis_config <- function(sp, Sig, dstar, B) {
  ## working statistics of 36 under a Gaussian truth with covariance Sig in the coordinates of sp (full rank)
  ell <- sp$ell; xi <- sp$xi; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  stopifnot(sp$nperp == 0, sp$nrest == 0, r == nu)
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2)
  R <- chol(Sig); g <- backsolve(R, forwardsolve(t(R), xi)); S0 <- sum(xi * g)   # S0 = xi' Sigma^-1 xi
  z <- crossprod(R, matrix(rnorm(r * B), r, B))              # draws of s under H0, covariance Sig
  bz <- colSums(g * z) / S0; Rk <- z - outer(xi, bz); rm(z)  # r = s - xi b, b the GLS estimate under Sig
  Rp <- numeric(B); ssr <- numeric(B); Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2); B_r <- crossprod(xi * Rk, W2); RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1  <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, GRID2)
  L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+"); ml1 <- apply(L1ml, 1, max)
  rm(A_r, B_r, RSS_r, l1, L1ml)
  at_r <- function(dl) { iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Rk2 * iD); Bx <- colSums(xi * Rk * iD); S2 <- colSums(xi^2 * iD)
    list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2) }
  Wd <- at_r(dh1); sig_h <- sqrt(Wd$RSS / (nu - 1))
  wa <- sqrt(Wd$S2 / S0) / sig_h; wg <- (Wd$Bx / Wd$S2) * sqrt(Wd$S2) / sig_h
  Os <- at_r(rep(dstar, B)); sig_s <- sqrt(Os$RSS / (nu - 1))                    # Wald at the pseudo-true ratio
  was <- sqrt(Os$S2 / S0) / sig_s; wgs <- (Os$Bx / Os$S2) * sqrt(Os$S2) / sig_s
  ev <- function(u, what) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); sk2 <- sk * sk
    A2 <- crossprod(sk2, W2); out <- list()
    if ("lrt" %in% what) out$lrt <- ml1 - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max)
    if ("Q" %in% what) {
      dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2); iD <- 1 / (1 + outer(ell, dh0))
      Ah <- colSums(sk2 * iD); Bh <- colSums(xi * sk * iD); Sh <- colSums(xi^2 * iD)
      out$Q <- Bh / sqrt(Sh * Ah / nu); out$dh0 <- dh0
    }
    out
  }
  lo <- rep(-30, B); hi <- rep(30, B)
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  dh0_null <- ev(rnorm(B), "Q")$dh0                           # null estimate on one null dataset per draw
  rootsU <- function(c0) {                                    # as in 36
    f <- function(u) ev(u, "lrt")$lrt
    always <- f(split) >= c0
    lo <- split; hi <- rep(120, B); never_r <- f(hi) < c0
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- f(mid) < c0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    ur <- (lo + hi) / 2; ur[never_r] <- Inf
    lo <- rep(-120, B); hi <- split; never_l <- f(lo) < c0
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; dn <- f(mid) < c0; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
    ul <- (lo + hi) / 2; ul[never_l] <- -Inf
    ur[always] <- split[always]; ul[always] <- split[always]
    list(up = ur, dn = ul)
  }
  rootsQ <- function(c0) {                                    # as in 36b
    one <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
      for (it in seq_len(NBIS + 6)) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      (lo + hi) / 2 }
    list(up = one(c0), dn = one(-c0))
  }
  lin <- function(a, g0, c0) list(up = (c0 - g0) / a, dn = (-c0 - g0) / a)
  perdraw <- function(rt) pnorm(-rt$up) + pnorm(rt$dn)
  fW <- crt_factors(sp, dh1, reml = TRUE)
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  res <- list()
  for (k in seq_along(AL)) {
    qc <- tq_[k] * sqrt(nu / (nu - 1 + tq_[k]^2))
    P <- list(oracle_star = perdraw(lin(was, wgs, tq_[k])), wald = perdraw(lin(wa, wg, tq_[k])),
              emmax = perdraw(rootsQ(qc)), lrt = perdraw(rootsU(cq[k])),
              kr = perdraw(lin(wa, wg, qt(1 - AL[k] / 2, df = fW$df) * sqrt(fW$ratio))))
    for (mt in names(P)) res[[length(res) + 1]] <- data.frame(method = mt, alpha = AL[k], size_ratio = mean(P[[mt]]) / AL[k],
      se_log10 = sd(P[[mt]]) / (sqrt(B) * mean(P[[mt]]) * log(10)))
  }
  list(sizes = do.call(rbind, res), med_dh0 = median(dh0_null), med_dh1 = median(dh1))
}

## pseudo-true ratio of the working null model: maximizer of -1/2 {sum log(1 + delta l) + nu log tr(W(delta) Sigma)}
pseudo_true <- function(sp, Sig) {
  grid <- 10^seq(-4, 4, by = 0.01); dS <- diag(Sig)
  f <- vapply(grid, function(dl) -0.5 * (sum(log1p(dl * sp$ell)) + sp$nu0 * log(sum(dS / (1 + dl * sp$ell)))), 0)
  grid[which.max(f)]
}
WIN <- list()
get_win <- function(i) { if (is.null(WIN[[as.character(i)]])) WIN[[as.character(i)]] <<- load_geu(i); WIN[[as.character(i)]] }
run_job <- function(k) {
  jb <- JOBS[k, ]; ck <- file.path(OUTM, paste0(jb$id, ".rds")); if (file.exists(ck)) return(paste("cached", jb$id))
  t0 <- Sys.time(); W <- get_win(jb$i); j <- match(jb$pos, W$pos); stopifnot(!is.na(j))
  sp <- spectral_setup(W$Z, j = j); U <- sp$U                 # eigenvectors of the working kernel, centred space
  Lw <- sp$ell
  Sig <- switch(jb$scen,
    correct = diag(1 + 0.25 * Lw),
    sparse  = { set.seed(4200 + jb$i); pool <- setdiff(seq_len(ncol(W$Z)), j); S <- sample(pool, NSP)
                ZS <- W$Z[, S, drop = FALSE]; KS <- crossprod(U, ZS) ; diag(length(Lw)) + 0.25 * tcrossprod(KS) / NSP },
    outside = { Zo <- do.call(cbind, lapply(setdiff(seq_len(nrow(GENES)), jb$i), function(i2) get_win(i2)$Z))
                Ko <- crossprod(U, Zo); diag(Lw * 0.125 + 1) + 0.125 * tcrossprod(Ko) / ncol(Zo) })
  Sig <- (Sig + t(Sig)) / 2
  ds <- pseudo_true(sp, Sig)
  set.seed(42000 + k)
  m <- mis_config(sp, Sig, ds, B42)
  ix <- indices(sp, ds)
  out <- cbind(data.frame(id = jb$id, gene = jb$gene, pos = jb$pos, q = jb$q, scen = jb$scen, nu = sp$nu0,
                          delta_star = ds, med_dh0 = m$med_dh0, med_dh1 = m$med_dh1, nu_j_star = ix[["nu_delta"]],
                          nu_j_true = if (jb$scen == "correct") indices(sp, 0.25)[["nu_delta"]] else NA), m$sizes)
  saveRDS(list(sizes = out, secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))), ck)
  sprintf("done %s (%.0f s)", jb$id, as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
cat("jobs:", nrow(JOBS), "\n")
suppressPackageStartupMessages(library(parallel))
print(unlist(mclapply(seq_len(nrow(JOBS)), function(k) tryCatch(run_job(k), error = function(e) paste("FAIL", JOBS$id[k], conditionMessage(e))),
                      mc.cores = NC, mc.preschedule = FALSE)))
M <- do.call(rbind, lapply(list.files(OUTM, "rds$", full.names = TRUE), function(f) readRDS(f)$sizes))
write.csv(M, paste0(OUTM, ".csv"), row.names = FALSE)
options(width = 200)
print(aggregate(size_ratio ~ scen + method + alpha, M, function(z) round(c(med = median(z), min = min(z), max = max(z)), 4)))
