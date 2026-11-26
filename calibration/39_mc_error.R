# ==============================================================================
# 39_mc_error.R
#
# Monte Carlo error of the exact tail probabilities of the calibration study
# (36_calibration_enlarged.R, 36b_score_reference.R), on seven difficult configurations.
# Each exact probability is the mean over B draws of the error contrasts r of a
# conditional probability (Rao-Blackwellized estimator); its standard error is the
# standard deviation of these per-draw probabilities divided by sqrt(B). For each
# configuration, the computation is repeated with B = 10,000 fresh draws, and the
# script reports, at alpha = 1e-3, 1e-5, 5e-8:
#   ratio_1k   size ratio of the stored B = 1,000 run (36, 36b)
#   ratio_10k  size ratio with B = 10,000
#   se_log10   standard error of log10(size ratio) for B = 1,000, from the per-draw
#              variance of the B = 10,000 run (delta method: sd / (sqrt(1000) mean ln 10))
#   ess        effective number of draws, (sum w)^2 / sum w^2 over the 10,000 draws
# for the Wald test at delta_hat, the test at the null estimate t(delta_hat_0), the
# likelihood ratio and the Kenward-Roger test.
# Usage (from the repository root): nice -n 10 Rscript calibration/39_mc_error.R [cores=2]
# Output: results/mc_error/<config>.rds, results/mc_error.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
stop_at <- grep("^OUT <- ", src)[1] - 1L                    # functions, genes and configurations of 36
eval(parse(text = src[seq_len(stop_at)]))
args <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
SMOKE <- FALSE
B_BIG <- 10000L

draws_config <- function(sp, delta0, B = B_BIG) {
  ## same algebra as calib_config() in 36; returns the per-draw conditional rejection probabilities
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1  <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+")
  dh1 <- pick_grid(l1, GRID2)
  L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+"); ml1 <- apply(L1ml, 1, max)
  rm(A_r, B_r, RSS_r, l1, L1ml)
  at_r <- function(dl) {
    iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Rk2 * iD) + Rp^2 + ssr; Bx <- colSums(xi * Rk * iD) + xp * Rp
    S2 <- colSums(xi^2 * iD) + xp^2
    list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
  }
  Wd <- at_r(dh1); sig_h <- sqrt(Wd$RSS / (nu - 1))
  wa <- sqrt(Wd$S2 / S0) / sig_h; wg <- (Wd$Bx / Wd$S2) * sqrt(Wd$S2) / sig_h
  ev <- function(u, what) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    A2 <- crossprod(sk2, W2) + (s2p^2 + ssr); out <- list()
    if ("lrt" %in% what) out$lrt <- ml1 - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max)
    if ("Q" %in% what) {
      dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
      iD <- 1 / (1 + outer(ell, dh0))
      Ah <- colSums(sk2 * iD) + s2p^2 + ssr; Bh <- colSums(xi * sk * iD) + xp * s2p
      Sh <- colSums(xi^2 * iD) + xp^2
      out$Q <- Bh / sqrt(Sh * Ah / nu)
    }
    out
  }
  lo <- rep(-30, B); hi <- rep(30, B)
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  rootsU <- function(c0) {                                   # likelihood ratio, as in 36
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
  rootsQ <- function(c0) {                                   # as in 36b
    one <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
      for (it in seq_len(NBIS + 6)) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < level
        lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      (lo + hi) / 2 }
    list(up = one(c0), dn = one(-c0))
  }
  lin <- function(a, g, c0) list(up = (c0 - g) / a, dn = (-c0 - g) / a)
  perdraw <- function(rt) pnorm(-rt$up) + pnorm(rt$dn)      # conditional rejection probability under H0
  fW <- crt_factors(sp, dh1, reml = TRUE)
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  out <- list()
  for (k in seq_along(AL)) {
    qc <- tq_[k] * sqrt(nu / (nu - 1 + tq_[k]^2))
    out[[sprintf("wald_t_%g", AL[k])]]   <- perdraw(lin(wa, wg, tq_[k]))
    out[[sprintf("kr_%g", AL[k])]]       <- perdraw(lin(wa, wg, qt(1 - AL[k] / 2, df = fW$df) * sqrt(fW$ratio)))
    out[[sprintf("lrt_chi2_%g", AL[k])]] <- perdraw(rootsU(cq[k]))
    out[[sprintf("emmax_t_%g", AL[k])]]  <- perdraw(rootsQ(qc))
  }
  out
}

CHOSEN <- c("geu_SNHG5_85909462_h0.05", "geu_ERAP2_96252589_h0.05", "syn_n1000_b5_h0.05", "syn_n200_b20_h0.05_unlinked")
## complete with GEUVADIS configurations near nu_delta = 25, 60 and 150, read from the stored results
C1 <- read.csv("calibration/results/calib_enlarged.csv")
U <- unique(C1[grepl("^geu_", C1$config), c("config", "nu_delta")])
for (tg in c(25, 60, 150)) CHOSEN <- c(CHOSEN, U$config[which.min(abs(log(U$nu_delta / tg)))])
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else
  sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)
ids <- vapply(cfgs, cfg_id, "")
sel <- cfgs[match(CHOSEN, ids)]; stopifnot(!anyNA(match(CHOSEN, ids)))
OUTD <- "calibration/results/mc_error"; dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
run_mc <- function(z) {
  id <- cfg_id(z); ck <- file.path(OUTD, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
  if (z$kind == "syn") {
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    sp <- if (z$unlinked) { set.seed(77 + z$n); spectral_setup(Z, x = rnorm(z$n)) } else spectral_setup(Z, j = 25L)
  } else {
    W <- load_geu(z$i); sp <- spectral_setup(W$Z, j = match(z$pos, W$pos))
  }
  set.seed(39000 + nchar(id) * 7L + sum(utf8ToInt(id)))   # fresh draws, independent of 36
  t0 <- Sys.time(); P <- draws_config(sp, z$h2 / (1 - z$h2))
  res <- do.call(rbind, lapply(names(P), function(nm) { w <- P[[nm]]; parts <- regmatches(nm, regexec("^(.*)_([^_]+)$", nm))[[1]]
    a <- as.numeric(parts[3]); data.frame(config = id, method = parts[2], alpha = a, ratio_10k = mean(w) / a,
      se_log10_1k = sd(w) / (sqrt(1000) * mean(w) * log(10)), se_log10_10k = sd(w) / (sqrt(length(w)) * mean(w) * log(10)),
      ess = sum(w)^2 / sum(w^2)) }))
  attr(res, "secs") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(res, ck); ck
}
suppressPackageStartupMessages(library(parallel))
st <- mclapply(sel, function(z) tryCatch(run_mc(z), error = function(e) paste("ERROR", conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st)[grepl("ERROR", unlist(st))])
R <- do.call(rbind, lapply(list.files(OUTD, "rds$", full.names = TRUE), readRDS))
S1 <- read.csv("calibration/results/calib_enlarged_score.csv")
ref <- rbind(C1[, c("config", "nu_delta", "method", "alpha", "size_ratio")], S1[S1$method == "emmax_t", c("config", "nu_delta", "method", "alpha", "size_ratio")])
M <- merge(R, ref, by = c("config", "method", "alpha"))
names(M)[names(M) == "size_ratio"] <- "ratio_1k"
M$diff_log10 <- log10(M$ratio_10k) - log10(M$ratio_1k)
write.csv(M, "calibration/results/mc_error.csv", row.names = FALSE)
options(width = 200)
print(M[order(M$nu_delta, M$method, M$alpha), c("config", "nu_delta", "method", "alpha", "ratio_1k", "ratio_10k", "diff_log10", "se_log10_1k", "ess")],
      digits = 3, row.names = FALSE)
