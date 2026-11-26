# ==============================================================================
# 10_clr_test.R
#
# Conditional likelihood-ratio (CLR) test, analogue of Moreira (2003).
# Under H0 the spectral coordinates s_i are independent and symmetric whatever
# (delta, sigma2), so conditionally on (|s_i|) and on the sum of squares of the
# directions free of the candidate, the signs are iid Rademacher. The profile-ML
# likelihood-ratio statistic
#   Lambda(s) = max_d[-sum log(1+d l) - n log RSS(d)] - max_d[-sum log(1+d l) - n log A(d)]
# depends on the signs only through B(d) = sum_i w_i(d) xi_i s_i (A is sign-invariant),
# so its conditional null law is that of Lambda(eps) with eps Rademacher. The CLR test
# rejects when P_eps(Lambda(eps) >= Lambda_obs) <= alpha: EXACT and similar for any delta.
#
# Tail probabilities down to 5e-8 are estimated by importance sampling: a symmetric
# mixture of two exponential tilts of the signs, P(eps_i = +-1) prop. to exp(+-theta a_i),
# with a_i = w_i(d*) xi_i |s_i| at the delta d* maximising the H1 profile on the observed
# data, and theta the saddlepoint of sum a_i eps_i at its observed absolute value.
#
# Power is estimated by direct simulation of datasets under H1 (u ~ N(lambda, 1),
# lambda = t_{nu-1, alpha}, as in 06/08), with the other treatments computed on the SAME
# datasets (paired): oracle t(delta0), LRT with the delta-known threshold, marginal BF with
# the half-Cauchy prior on tau (08), sign-flip on the score numerator (06), Wald plug-in.
# Size under H0 is checked at alpha = 0.05 and 0.01 (exactness holds by construction).
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/10_clr_test.R [cores=2] [--smoke]
# Output: results/clr_test/<config>.rds, results/clr_test.csv
# ==============================================================================

args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
SMOKE  <- "--smoke" %in% args
src <- readLines("calibration/08_mbf_prior_fix.R")
cut <- grep("^## ---- same configurations", src)
eval(parse(text = src[1:(cut - 1)]))           # helpers of 06 and 08 (node_pairs, logp1_mix, ...)
GRIDL <- c(0, 10^seq(-4, 4, length.out = 41))   # delta grid of the LRT statistic
N0 <- if (SMOKE) 200L else 1000L                 # datasets under H0 (size check)
N1 <- if (SMOKE) 100L else 400L                  # datasets per power point
MIS <- if (SMOKE) 2000L else 3000L               # importance samples per p-value

clr_config <- function(sp, delta0, tag, h0_cache = NULL, h1_seed = 1L) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  has_p <- sp$nperp == 1L
  WL <- 1 / (1 + outer(ell, GRIDL)); ldL <- colSums(log(1 + outer(ell, GRIDL)))
  SL <- colSums(xi^2 * WL) + xp^2
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  WQ <- 1 / (1 + outer(ell, QG)); ldQ <- colSums(log(1 + outer(ell, QG)))
  SQg <- colSums(xi^2 * WQ) + xp^2
  tq <- log(QG); dq <- tq[2] - tq[1]
  lwq <- log(c(dq / 2, rep(dq, length(tq) - 2), dq / 2)) + 0.5 * tq - log(pi) - log1p(QG)
  NPhc <- node_pairs(SQg, "hc")
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  tq_ <- qt(1 - AL / 2, df = nu - 1); qa <- tq_ * sqrt(nu / (nu - 1 + tq_^2))
  cLR <- n * log1p(tq_^2 / (nu - 1))
  cHC <- vapply(qa, function(q) lbf_oracle(q^2 * S0 / nu, S0, nu, "hc"), 0)
  ## coefficient matrix of B(delta) in the signs, K x G (kernel directions, then x_perp)
  coefL <- if (has_p) rbind(xi * WL, rep(xp, length(GRIDL))) else xi * WL
  one <- function(s, spv, ssr) {
    as_ <- abs(c(s, if (has_p) spv))
    sg  <- sign(c(s, if (has_p) spv))
    aM  <- coefL * as_                                   # K x G, a_i(delta)
    Ag  <- colSums(s^2 * WL) + spv^2 + ssr
    Bo  <- colSums(sg * aM)
    h1o <- -ldL - n * log(Ag - Bo^2 / SL)
    m0  <- max(-ldL - n * log(Ag))
    lam_obs <- max(h1o) - m0
    ## ---- importance sampling of the sign-flip tail ----
    gs <- which.max(h1o); a <- aM[, gs]; tobs <- abs(Bo[gs]); K <- length(a)
    th <- if (tobs > 0) { t <- tobs / sum(a^2)
      for (it in 1:60) { thv <- tanh(a * t); K1 <- sum(a * thv); K2 <- sum(a^2 * (1 - thv^2))
        tn <- t - (K1 - tobs) / max(K2, 1e-300); if (!is.finite(tn) || tn <= 0) tn <- t * 0.5
        if (K1 < tobs && tn < t) tn <- 2 * t
        tn <- min(tn, 30 / max(abs(a)))
        if (abs(tn - t) < 1e-10 * max(t, 1e-12)) { t <- tn; break }; t <- tn }
      t } else 0
    pp <- plogis(2 * th * a)                             # P(eps_i = +1) under the + tilt
    lpp <- plogis(2 * th * a, log.p = TRUE)              # stable logs: never -Inf, so no 0 * -Inf
    lpm <- plogis(-2 * th * a, log.p = TRUE)             # in the indicator products below
    half <- MIS %/% 2
    U <- matrix(runif(MIS * K), MIS, K)
    E <- ifelse(U < rep(pp, each = MIS), 1, -1)
    E[(half + 1):MIS, ] <- -E[(half + 1):MIS, ]          # mirrored half = the - tilt
    lqp <- (E > 0) %*% lpp + (E < 0) %*% lpm             # log q_+ (eps)
    lqm <- (E < 0) %*% lpp + (E > 0) %*% lpm             # log q_- (eps)
    mx <- pmax(lqp, lqm); lq <- log(0.5) + mx + log(exp(lqp - mx) + exp(lqm - mx))
    lw <- -K * log(2) - lq
    Bm <- E %*% aM                                       # MIS x G
    lamE <- apply(sweep(-n * log(pmax(rep(Ag, each = MIS) - sweep(Bm^2, 2, SL, "/"), 1e-300)),
                        2, ldL, "-"), 1, max) - m0
    hit <- lamE >= lam_obs - 1e-9
    lp <- if (any(hit)) { v <- lw[hit]; mv <- max(v); mv + log(sum(exp(v - mv))) - log(MIS) } else -Inf
    ## ---- comparators on the same dataset ----
    A2 <- colSums(s^2 * W2) + spv^2 + ssr; B2 <- colSums(xi * s * W2) + xp * spv
    l1 <- -0.5 * ((nu - 1) * log(A2 - B2^2 / S2g) + ld2 + log(S2g))
    l0 <- -0.5 * (nu * log(A2) + ld2)
    at <- function(dl) { w <- 1 / (1 + dl * ell)
      A <- sum(s^2 * w) + spv^2 + ssr; Bx <- sum(xi * s * w) + xp * spv; S <- sum(xi^2 * w) + xp^2
      c(A = A, B = Bx, S = S) }
    dh1 <- pick_grid(matrix(l1, 1), GRID2); dh0 <- pick_grid(matrix(l0, 1), GRID2)
    q1 <- at(dh1); t_w <- (q1["B"] / q1["S"]) * sqrt(q1["S"] * (nu - 1) / (q1["A"] - q1["B"]^2 / q1["S"]))
    q0 <- at(delta0); t_or <- (q0["B"] / q0["S"]) * sqrt(q0["S"] * (nu - 1) / (q0["A"] - q0["B"]^2 / q0["S"]))
    qh <- at(dh0); w0 <- 1 / (1 + dh0 * ell)
    lsf <- log_signflip_p(matrix(c(abs(xi * w0 * s), abs(xp * spv)), ncol = 1), abs(qh["B"]))
    AQ <- colSums(s^2 * WQ) + spv^2 + ssr; BQ <- colSums(xi * s * WQ) + xp * spv
    l0w <- -0.5 * (nu * log(AQ) + ldQ) + lwq
    lmbf <- logp1_mix(matrix(BQ^2 / AQ, 1), matrix(l0w, 1), NPhc, nu) - rowLSE(matrix(l0w, 1))
    c(clr_lp = lp, lrt = lam_obs, t_wald = unname(t_w), t_or = unname(t_or),
      lsf = unname(lsf), lmbf = unname(lmbf), ess = sum(exp(lw[hit]))^2 / max(sum(exp(2 * lw[hit])), 1e-300))
  }
  sim <- function(lam, N) {
    beta <- lam / sqrt(S0)
    Z <- matrix(rnorm(r * N), r, N) * sqrt(d0) + xi * beta
    P <- if (has_p) rnorm(N) + xp * beta else numeric(N)
    R <- if (sp$nrest > 0) rchisq(N, sp$nrest) else numeric(N)
    t(vapply(seq_len(N), function(j) one(Z[, j], P[j], R[j]),
             setNames(numeric(7), c("clr_lp", "lrt", "t_wald", "t_or", "lsf", "lmbf", "ess"))))
  }
  res <- list()
  if (!is.null(h0_cache) && file.exists(h0_cache)) {
    res[[1]] <- readRDS(h0_cache)
  } else {
    H0 <- sim(0, N0); r0 <- list()
    for (a0 in c(0.05, 0.01)) {
      ct <- qt(1 - a0 / 2, nu - 1)
      r0[[length(r0) + 1]] <- data.frame(tag, lambda = 0, alpha = a0,
        method = c("clr", "lrt_ml", "oracle_t", "wald_h1", "signflip"),
        reject = c(mean(H0[, "clr_lp"] <= log(a0)), mean(H0[, "lrt"] > n * log1p(ct^2 / (nu - 1))),
                   mean(abs(H0[, "t_or"]) > ct), mean(abs(H0[, "t_wald"]) > ct), mean(H0[, "lsf"] <= log(a0))),
        n_sim = N0, na_share = mean(!is.finite(H0[, "clr_lp"]) & !is.infinite(H0[, "clr_lp"])))
    }
    res[[1]] <- do.call(rbind, r0)
    if (!is.null(h0_cache)) saveRDS(res[[1]], h0_cache)
  }
  set.seed(h1_seed)
  for (k in seq_along(AL)) {
    H1 <- sim(tq_[k], N1)
    res[[length(res) + 1]] <- data.frame(tag, lambda = tq_[k], alpha = AL[k],
      method = c("clr", "lrt_ml", "mbf_hc", "oracle_t", "wald_h1", "signflip"),
      reject = c(mean(H1[, "clr_lp"] <= log(AL[k])), mean(H1[, "lrt"] > cLR[k]),
                 mean(H1[, "lmbf"] > cHC[k]), mean(abs(H1[, "t_or"]) > tq_[k]),
                 mean(abs(H1[, "t_wald"]) > tq_[k]), mean(H1[, "lsf"] <= log(AL[k]))),
      n_sim = N1, na_share = mean(is.nan(H1[, "clr_lp"])))
  }
  out <- do.call(rbind, res)
  attr(out, "ess_median") <- NA
  out
}

## ---- configurations (same seeds for the genotype designs as 06) --------------------
cfgs <- list()
for (i in seq_len(8)) for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "geu", i = i, h2 = h2)
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2)
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t",
                fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, j = match(g$lead_pos, pos))
}
OUT <- "calibration/results/clr_test"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", n = z$n, blocks = z$nb,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(9100 + 1000 * z$n + z$nb + round(100 * z$h2))
  } else {
    id <- sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", n = nrow(L$Z), blocks = NA,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(9900 + z$i + round(100 * z$h2))
  }
  h0c <- file.path(OUT, "h0", paste0(id, ".rds")); dir.create(dirname(h0c), showWarnings = FALSE)
  out <- clr_config(sp, z$h2 / (1 - z$h2), tag, h0_cache = if (SMOKE) NULL else h0c,
                    h1_seed = 777L + (if (z$kind == "syn") 1000L * z$n + z$nb else z$i) + round(100 * z$h2))
  if (!SMOKE) saveRDS(out, ck)
  cat(sprintf("[%7.1f s] %s done (median IS ess under H0: %.0f)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), id, attr(out, "ess_median")))
  if (SMOKE) out else ck
}
if (SMOKE) {
  out <- run_one(cfgs[[1]])
  op <- options(width = 160); print(out[, c("lambda", "alpha", "method", "reject", "n_sim")], digits = 3)
  options(op)
} else {
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e)
    paste("ERROR", conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fx <- do.call(rbind, lapply(list.files(OUT, "rds$", full.names = TRUE), readRDS))
  write.csv(fx, "calibration/results/clr_test.csv", row.names = FALSE)
  cat("configs done:", length(unique(fx$config)), "of", length(cfgs), "\n")
}
cat("done in", format(Sys.time() - t0), "\n")
