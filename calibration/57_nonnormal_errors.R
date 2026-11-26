# ==============================================================================
# 57_nonnormal_errors.R
# Sensitivity of the sizes to non-Gaussian errors, by direct simulation on four configurations
# of the calibration study (36_calibration_enlarged.R). Each null dataset is y = v + e with a
# Gaussian background v ~ N(0, delta_0 K) and independent errors e of unit variance drawn from
#   normal   N(0, 1)
#   t5       Student with 5 degrees of freedom, scaled to unit variance (heavy tails)
#   gamma2   gamma with shape 2, centred and scaled to unit variance (skewness 1.41)
# The errors are drawn in R^n and projected on the spectral coordinates of the design, so that
# nothing in the computation uses their distribution. Every statistic is computed on the full
# dataset, without the decomposition (b, r):
#   known     t(delta_0), Student nu-1
#   wald      t(delta_hat), REML with the candidate, Student nu-1
#   null      t(delta_hat_0), REML of the null model, Student nu-1
#   kr        Kenward-Roger
#   lrt       likelihood ratio, ML under each model, chi2_1
#   rlr       restricted likelihood ratio on a common residual space (56), chi2_1
# Rejections are counted at alpha in {1e-3, 1e-4, 1e-5}.
# Usage (from the repository root): nice -n 10 Rscript calibration/57_nonnormal_errors.R [cores=1] [datasets=4e6]
# Output: results/nonnormal_errors/<config>_<error>.rds, results/nonnormal_errors.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # GENES, cfgs, load_geu, gen_blocks, GRID2, pick_grid, crt_factors
a0 <- commandArgs(TRUE); NC <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 1L
NDIR <- if (length(a0) >= 2) as.integer(as.numeric(a0[2])) else 4000000L
suppressPackageStartupMessages(library(parallel))
ALS <- c(1e-3, 1e-4, 1e-5); CH <- 4000L
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else
  sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)
ids <- vapply(cfgs, cfg_id, "")
WANT <- c("geu_SNHG5_86504582_h0.05", "geu_ERAP2_96252589_h0.05", "syn_n200_b5_h0.05", "geu_ERAP2_95943006_h0.05")
ERRS <- list(normal = function(n) rnorm(n),
             t5     = function(n) rt(n, 5) * sqrt(3 / 5),
             gamma2 = function(n) (rgamma(n, shape = 2) - 2) / sqrt(2))
TESTS <- c("known", "wald", "null", "kr", "lrt", "rlr")
OUT <- "calibration/results/nonnormal_errors"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

one <- function(id, err) {
  ck <- file.path(OUT, sprintf("%s_%s.rds", id, err)); if (file.exists(ck)) return(readRDS(ck))
  k <- match(id, ids); z <- cfgs[[k]]
  if (z$kind == "syn") { set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95); j <- 25L
  } else { w <- load_geu(z$i); Z <- w$Z; j <- match(z$pos, w$pos) }
  sp <- spectral_setup(Z, j = j); n <- nrow(Z)
  d0 <- z$h2 / (1 - z$h2); ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; U <- sp$U
  stopifnot(max(abs(colSums(U))) < 1e-6)                    # kernel directions orthogonal to the intercept
  x <- Z[, j] - mean(Z[, j]); q <- if (sp$nperp) (x - U %*% xi) / xp else NULL
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  tq <- qt(1 - ALS / 2, nu - 1); cq <- qchisq(1 - ALS, 1)
  stat_at <- function(Sk, Sk2, Sp, ssr, dl) {               # t statistic at the ratios dl (one per dataset)
    iD <- 1 / (1 + outer(ell, dl)); A <- colSums(Sk2 * iD) + Sp^2 + ssr; Bx <- colSums(xi * Sk * iD) + xp * Sp
    S2 <- colSums(xi^2 * iD) + xp^2; (Bx / sqrt(S2)) / sqrt((A - Bx^2 / S2) / (nu - 1))
  }
  chunk <- function(qq) {
    set.seed(57000000 + 10000 * k + 1000 * match(err, names(ERRS)) + qq)
    E <- matrix(ERRS[[err]](n * CH), n, CH); E <- sweep(E, 2, colMeans(E))
    Ue <- crossprod(U, E); Sp <- if (sp$nperp) drop(crossprod(q, E)) else numeric(CH)
    ssr <- pmax(colSums(E^2) - colSums(Ue^2) - Sp^2, 0); if (sp$nrest == 0) ssr[] <- 0
    Sk <- Ue + matrix(rnorm(r * CH), r, CH) * sqrt(d0 * ell); rm(E, Ue); Sk2 <- Sk * Sk
    A <- crossprod(Sk2, W2) + (Sp^2 + ssr); Bm <- crossprod(xi * Sk, W2) + xp * Sp
    lA <- log(A); lR <- log(A - sweep(Bm^2, 2, S2g, "/"))
    l1 <- -0.5 * sweep((nu - 1) * lR, 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, GRID2)
    dh0 <- pick_grid(-0.5 * sweep(nu * lA, 2, ld2, "+"), GRID2)
    lrt <- apply(-sweep((nu + 1) * lR, 2, ld2, "+"), 1, max) - apply(-sweep((nu + 1) * lA, 2, ld2, "+"), 1, max)
    rlr <- 2 * (apply(l1, 1, max) - apply(-0.5 * sweep((nu - 1) * lA, 2, ld2 + log(S2g), "+"), 1, max))
    tw <- abs(stat_at(Sk, Sk2, Sp, ssr, dh1)); fW <- crt_factors(sp, dh1, reml = TRUE)
    t0 <- abs(stat_at(Sk, Sk2, Sp, ssr, dh0)); tk <- abs(stat_at(Sk, Sk2, Sp, ssr, rep(d0, CH)))
    vapply(seq_along(ALS), function(a) c(sum(tk > tq[a]), sum(tw > tq[a]), sum(t0 > tq[a]),
      sum(tw > qt(1 - ALS[a] / 2, df = fW$df) * sqrt(fW$ratio)), sum(lrt > cq[a]), sum(rlr > cq[a])), numeric(length(TESTS)))
  }
  nch <- NDIR %/% CH
  R <- Reduce(`+`, mclapply(seq_len(nch), chunk, mc.cores = NC))
  N <- nch * CH; p <- as.vector(R) / N
  out <- data.frame(config = id, nu_delta = indices(sp, d0)[["nu_delta"]], n = n, error = err,
                    test = rep(TESTS, length(ALS)), alpha = rep(ALS, each = length(TESTS)), datasets = N,
                    rejections = as.vector(R), size_ratio = p / rep(ALS, each = length(TESTS)),
                    se = sqrt(p * (1 - p) / N) / rep(ALS, each = length(TESTS)))
  saveRDS(out, ck); out
}
for (id in WANT) for (err in names(ERRS)) { t0 <- Sys.time(); o <- one(id, err)
  cat(sprintf("%s %s done (%.1f min)\n", id, err, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  write.csv(do.call(rbind, lapply(list.files(OUT, "rds$", full.names = TRUE), readRDS)),
            "calibration/results/nonnormal_errors.csv", row.names = FALSE) }
