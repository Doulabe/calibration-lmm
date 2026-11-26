# ==============================================================================
# 33c_fix_ml_term.R
# Corrects the maximum-likelihood statistics of the variant reports at conditional
# steps (k >= 2). With fixed effects F, log|V| = log|H'VH| - log|F'V^{-1}F| + log|F'F|:
#   - 33_geuvadis_comparators.R added log|F'V^{-1}F| instead of subtracting it, which
#     changed GEMMA's likelihood-ratio test and the ML null estimate of its score test at
#     steps >= 2 (step 1 is unaffected: intercept only and centred kernel, so F'V^{-1}F is
#     constant in delta);
#   - 27_geuvadis_application.R omitted the term, so its likelihood ratio at steps >= 2 was
#     the reduced-spectrum statistic rather than the maximum-likelihood one.
# The restricted likelihood is unaffected (the term cancels), hence so are the Wald and
# Kenward-Roger columns, which are copied unchanged. This script recomputes, for every
# report row with step >= 2, the observed likelihood ratio (and GEMMA's score) and their
# exact p-values with the correct term, from the same data, seeds and delta grids.
# Usage (from the repository root): nice -n 10 Rscript calibration/33c_fix_ml_term.R [cores=2]
# Output: overwrites the $variants of results/geuvadis_comp/*.rds and results/geuvadis_app/*.rds
#         (originals kept as *_preMLfix.rds), and checks GEMMA agreement at steps >= 2.
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(TRUE); NCORES <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
suppressPackageStartupMessages(library(parallel))
source("calibration/00_functions.R")
B_EX <- 1000L; GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
RES <- "calibration/results"
TAGS <- c(ERAP2 = "ENSG00000164308", RPS26 = "ENSG00000197728", SLFN5 = "ENSG00000166750", SNHG5 = "ENSG00000203875",
          "FLVCR1-AS1" = "ENSG00000198468", "PEX6-region" = "ENSG00000124587", "TRA2A-AS" = "ENSG00000230658", ZNF266 = "ENSG00000174652")
pick_grid <- function(O, grid) {
  g <- max.col(O, ties.method = "first"); dh <- grid[g]
  int <- g > 2 & g < length(grid)
  if (any(int)) {
    ii <- which(int); lg <- log(grid)
    x1 <- lg[g[ii] - 1]; x2 <- lg[g[ii]]; x3 <- lg[g[ii] + 1]
    f1 <- O[cbind(ii, g[ii] - 1)]; f2 <- O[cbind(ii, g[ii])]; f3 <- O[cbind(ii, g[ii] + 1)]
    den <- (x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1)
    num <- (x2 - x1)^2 * (f2 - f3) - (x2 - x3)^2 * (f2 - f1)
    xs <- ifelse(abs(den) > 1e-14, x2 - 0.5 * num / den, x2)
    dh[ii] <- exp(pmin(pmax(xs, x1), x3))
  }
  dh
}
## reduced model at a step, kernel either fixed (GEMMA) or candidate-specific,
## with ldF(delta) = log|F'V(delta)^{-1}F| on GRID2 (to be SUBTRACTED from log|H'VH| for ML)
spectral_any <- function(X, y, sel, j, kernel = c("fixed", "specific")) {
  kernel <- match.arg(kernel); n <- nrow(X); m <- ncol(X)
  Kn <- if (kernel == "fixed") tcrossprod(X) / m else { Rk <- setdiff(seq_len(m), c(sel, j)); tcrossprod(X[, Rk, drop = FALSE]) / length(Rk) }
  F <- cbind(1, X[, sel, drop = FALSE]); p <- ncol(F)
  H <- qr.Q(qr(F), complete = TRUE)[, -(1:p), drop = FALSE]
  M <- crossprod(H, Kn %*% H); e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values); U <- e$vectors[, keep, drop = FALSE]; ell <- e$values[keep]
  xt <- as.numeric(crossprod(H, X[, j])); yt <- as.numeric(crossprod(H, y))
  xi <- as.numeric(crossprod(U, xt)); xperp <- xt - as.numeric(U %*% xi); xp <- sqrt(sum(xperp^2))
  sk <- as.numeric(crossprod(U, yt)); nperp <- as.integer(xp > 1e-6 * sqrt(sum(xt^2)))
  sp_ <- if (nperp) sum(xperp * yt) / xp else 0
  ssr <- max(sum(yt^2) - sum(sk^2) - sp_^2, 0); nu <- ncol(H)
  ek <- eigen(Kn, symmetric = TRUE); UF <- crossprod(ek$vectors, F); lk <- pmax(ek$values, 0)
  ldF <- vapply(GRID2, function(dl) as.numeric(determinant(crossprod(UF, UF / (1 + dl * lk)))$modulus), 0)
  list(ell = ell, xi = xi, xp = if (nperp) xp else 0, nu0 = nu, r = length(ell), nperp = nperp,
       nrest = nu - length(ell) - nperp, n = n, ldF = ldF - ldF[1], obs = list(sk = sk, sp = sp_, ssr = ssr))
}
## ML statistics with the correct fixed-effect term: likelihood ratio, ML null delta and GEMMA's score
ml_stats <- function(sp, Sk, Sp, Ssr) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- sp$n
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
  RSS2 <- A2 - sweep(B2^2, 2, S2g, "/")
  ml0 <- -0.5 * sweep(n * log(A2), 2, ld2 - sp$ldF, "+"); ml1 <- -0.5 * sweep(n * log(RSS2), 2, ld2 - sp$ldF, "+")
  dm0 <- pick_grid(ml0, GRID2)
  iD <- 1 / (1 + outer(ell, dm0))
  A <- colSums(Sk^2 * iD) + Sp^2 + Ssr; Bx <- colSums(xi * Sk * iD) + xp * Sp; S <- colSums(xi^2 * iD) + xp^2
  l1 <- -0.5 * sweep((nu - 1) * log(RSS2), 2, ld2 + log(S2g), "+")
  list(lrt = 2 * (apply(ml1, 1, max) - apply(ml0, 1, max)), Fs = n * Bx^2 / (S * A), dh1 = pick_grid(l1, GRID2), l1 = l1)
}
lmean_exp <- function(lp) (max(lp) + log(mean(exp(lp - max(lp))))) / log(10)
exact_ml <- function(sp, delta0, obs_lrt, obs_Fs, B = B_EX, NBIS = 40L) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); ml_stats(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  ## split point: where the weighted covariance of the candidate with the data changes sign, at the
  ## ML null delta (the likelihood ratio and the score are U-shaped in u around it)
  tq <- function(u) { bv <- u / sqrt(S0); z <- ev(u); iD <- 1 / (1 + outer(ell, pmax(z$dh1, 0)))
    colSums(xi * (Rk + outer(xi, bv)) * iD) + xp * (Rp + xp * bv) }
  split <- { lo <- rep(-60, B); hi <- rep(60, B)
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- tq(mid) < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
  two <- function(f, level) {
    ur <- { lo <- split; hi <- rep(250, B)
      for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- f(mid) < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
    ul <- { lo <- rep(-250, B); hi <- split
      for (it in 1:NBIS) { mid <- (lo + hi) / 2; dn <- f(mid) < level; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }; (lo + hi) / 2 }
    a1 <- pnorm(-ur, log.p = TRUE); a2 <- pnorm(ul, log.p = TRUE); lmean_exp(pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))) }
  c(lrt = two(function(u) ev(u)$lrt, obs_lrt), score = if (is.null(obs_Fs)) NA else two(function(u) ev(u)$Fs, obs_Fs))
}
## ---- data ----------------------------------------------------------------------------------------
loadZ <- function(pop, g) { Z <- as.matrix(read.table(file.path(RES, "geuvadis_comp", "data", sprintf("%s_%s_Z.txt", pop, g))))
  y <- scan(file.path(RES, "geuvadis_comp", "data", sprintf("%s_%s_y.txt", pop, g)), quiet = TRUE); list(Z = Z, y = y) }
jobs <- expand.grid(pop = c("EUR", "YRI"), g = names(TAGS), what = c("comp", "app"), stringsAsFactors = FALSE)
fix_one <- function(q) {
  pop <- jobs$pop[q]; g <- jobs$g[q]; what <- jobs$what[q]
  f <- file.path(RES, if (what == "comp") "geuvadis_comp" else "geuvadis_app", sprintf("%s_%s.rds", pop, TAGS[[g]]))
  r <- readRDS(f); v <- r$variants; if (is.null(v)) return(paste("no variants", what, pop, g))
  rows <- which(v$step >= 2); if (!length(rows)) return(paste("nothing to fix", what, pop, g))
  bak <- sub("\\.rds$", "_preMLfix.rds", f); if (!file.exists(bak)) file.copy(f, bak)
  W <- loadZ(pop, g); X <- W$Z; y <- W$y; stopifnot(ncol(X) == length(r$pos))
  set.seed(3390 + match(g, names(TAGS)) + if (pop == "YRI") 50L else 0L + if (what == "app") 500L else 0L)
  for (i in rows) {
    key <- rownames(v)[i]; parts <- strsplit(sub("^\\|", "|", key), ",")[[1]]; bar <- which(parts == "|")
    before <- as.integer(parts[seq_len(bar - 1)]); j <- as.integer(parts[bar + 1]); stopifnot(j == v$j[i])
    sp <- spectral_any(X, y, before, j, if (what == "comp") "fixed" else "specific")
    o <- ml_stats(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
    dh <- o$dh1; prof <- 2 * (max(o$l1) - o$l1[1, ]); ci <- range(GRID2[prof <= qchisq(0.95, 1)])
    dgrid <- unique(c(ci[1], dh, ci[2], if (what == "app") exp(mean(log(pmax(ci, 1e-4)))) else NULL))
    lp <- t(vapply(dgrid, function(d0) exact_ml(sp, d0, o$lrt, if (what == "comp") o$Fs else NULL), numeric(2)))
    if (what == "comp") {
      v$log10p_ours_lrt[i] <- pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
      v$log10p_ours_score[i] <- pf(o$Fs, 1, sp$nu0 - 1, lower.tail = FALSE, log.p = TRUE) / log(10)
      v$log10p_ex_lrt[i] <- lp[dgrid == dh, "lrt"][1]; v$log10p_exmax_lrt[i] <- max(lp[, "lrt"])
      v$log10p_ex_score[i] <- lp[dgrid == dh, "score"][1]; v$log10p_exmax_score[i] <- max(lp[, "score"])
    } else {
      v$lrt[i] <- o$lrt; v$log10p_nom_lrt[i] <- pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
      v$log10p_ex_lrt[i] <- lp[dgrid == dh, "lrt"][1]; v$log10p_exmax_lrt[i] <- max(lp[, "lrt"])
    }
  }
  v$ml_fixed <- TRUE; r$variants <- v; saveRDS(r, f)
  sprintf("fixed %d rows %s %s %s", length(rows), what, pop, g)
}
st <- mclapply(seq_len(nrow(jobs)), function(q) tryCatch(fix_one(q), error = function(e) paste("FAIL", jobs$what[q], jobs$pop[q], jobs$g[q], conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st))
## agreement with GEMMA after the fix (all steps)
V <- do.call(rbind, lapply(list.files(file.path(RES, "geuvadis_comp"), "^(EUR|YRI)_ENSG[0-9]+\\.rds$", full.names = TRUE),
  function(f) { r <- readRDS(f); if (is.null(r$variants)) NULL else cbind(pop = r$pop, gene = r$gene, r$variants) }))
cat("GEMMA vs reduction, max |delta log10 p|: LRT", max(abs(V$log10p_gemma_lrt - V$log10p_ours_lrt)),
    " score", max(abs(V$log10p_gemma_score - V$log10p_ours_score)), " Wald", max(abs(V$log10p_gemma_wald - V$log10p_ours_wald)), "\n")
print(V[V$step >= 2, c("pop", "gene", "step", "pos", "log10p_gemma_lrt", "log10p_ours_lrt", "log10p_gemma_score", "log10p_ours_score")], digits = 4)
