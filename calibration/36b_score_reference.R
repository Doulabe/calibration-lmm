# ==============================================================================
# 36b_score_reference.R
#
# Companion of 36_calibration_enlarged.R. In 36, the null-estimate treatment "score_t"
# compares the score Q(delta_hat_0), whose residual variance is that of the null model,
# with the Student quantile t_{nu-1}. At a fixed delta, Q = t * sqrt(nu / (nu - 1 + t^2))
# is smaller than the Wald statistic t of the same fit, so this reference is conservative
# even when delta is known (by about 0.06 at n = 200 and 0.6 at n = 1000 at 5e-8). That
# statistic is GEMMA's score test (n R^2 referred to F(1, n - c - 1)). EMMAX, P3D and each
# step of MLMM instead fit the candidate by generalized least squares at delta_hat_0 and
# refer t(delta_hat_0) to t_{nu-1}; since t and Q are monotone transforms of each other at
# the same delta, that test rejects when |Q(delta_hat_0)| exceeds
# q_c = c * sqrt(nu / (nu - 1 + c^2)), c = t_{nu-1, 1 - alpha/2}.
# This script recomputes, on exactly the same draws of the error contrasts as 36 (same
# seeds, same order of random draws), the exact size and power of
#   emmax_t   t(delta_hat_0) referred to t_{nu-1}   (EMMAX, P3D, MLMM step)
#   score_t   Q(delta_hat_0) referred to t_{nu-1}   (GEMMA score; reproduces 36 as a check)
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/36b_score_reference.R [cores=3]
# Output: results/calib_enlarged_score/<config>.rds, results/calib_enlarged_score.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
stop_at <- grep("^OUT <- ", src)[1] - 1L               # functions, genes and configurations of 36
eval(parse(text = src[seq_len(stop_at)]))
SMOKE <- FALSE

calib_score <- function(sp, delta0, tag, B = B_DRAW) {
  ## identical to the beginning of calib_config() in 36: same draws in the same order
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  evQ <- function(u) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    A2 <- crossprod(sk2, W2) + (s2p^2 + ssr)
    dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
    iD <- 1 / (1 + outer(ell, dh0))
    Ah <- colSums(sk2 * iD) + s2p^2 + ssr; Bh <- colSums(xi * sk * iD) + xp * s2p
    Sh <- colSums(xi^2 * iD) + xp^2
    Bh / sqrt(Sh * Ah / nu)
  }
  rootsQ <- function(c0) {                                   # as in 36
    one <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
      for (it in seq_len(NBIS + 6)) { mid <- (lo + hi) / 2; up <- evQ(mid) < level
        lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      (lo + hi) / 2 }
    list(up = one(c0), dn = one(-c0))
  }
  rej <- function(rt, lam) mean(pnorm(lam - rt$up) + pnorm(rt$dn - lam))
  tq_ <- qt(1 - AL / 2, df = nu - 1)
  res <- list()
  for (k in seq_along(AL)) {
    qc <- tq_[k] * sqrt(nu / (nu - 1 + tq_[k]^2))           # threshold of Q equivalent to |t| > c
    for (m in c("emmax_t", "score_t")) {
      rt <- rootsQ(if (m == "emmax_t") qc else tq_[k])
      res[[length(res) + 1]] <- data.frame(tag, method = m, alpha = AL[k], size_ratio = rej(rt, 0) / AL[k],
                                           power = rej(rt, tq_[k]), lambda = tq_[k])
    }
  }
  do.call(rbind, res)
}

OUT2 <- "calibration/results/calib_enlarged_score"
dir.create(OUT2, showWarnings = FALSE, recursive = TRUE)
run_score <- function(z) {                                   # seeds and set-up exactly as run_one() in 36
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "")
    ck <- file.path(OUT2, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L)
  } else {
    g <- GENES[z$i, ]; id <- sprintf("geu_%s_%d_h%.2f", g$symbol, z$pos, z$h2)
    ck <- file.path(OUT2, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973))
  }
  out <- calib_score(sp, z$h2 / (1 - z$h2), tag)
  saveRDS(out, ck); ck
}
## order: the synthetic configurations first (they contain the unlinked controls)
cfgs <- c(cfgs[vapply(cfgs, function(z) z$kind == "syn", TRUE)], cfgs[vapply(cfgs, function(z) z$kind == "geu", TRUE)])
cat("configurations:", length(cfgs), "\n")
st <- mclapply(cfgs, function(z) tryCatch(run_score(z), error = function(e) paste("ERROR", conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st)[grepl("ERROR", unlist(st))])
fs <- list.files(OUT2, "rds$", full.names = TRUE)
S <- do.call(rbind, lapply(fs, readRDS))
write.csv(S, "calibration/results/calib_enlarged_score.csv", row.names = FALSE)
## check: score_t must reproduce the values of 36 (same draws)
C <- read.csv("calibration/results/calib_enlarged.csv")
m <- merge(S[S$method == "score_t", c("config", "alpha", "size_ratio", "power")],
           C[C$method == "score_t", c("config", "alpha", "size_ratio", "power")], by = c("config", "alpha"))
cat("configs done:", length(fs), " max relative difference with 36 (size, power):",
    signif(max(abs(m$size_ratio.x / m$size_ratio.y - 1)), 3), signif(max(abs(m$power.x / m$power.y - 1)), 3), "\n")
