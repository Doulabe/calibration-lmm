# ==============================================================================
# 56_restricted_lr_adjusted_power.R
# Two additions to the calibration study of 36_calibration_enlarged.R, on the same configurations
# and the same draws of the error contrasts (same seeds, same order of the random numbers).
#  1. Restricted likelihood ratio on a common residual space (Welham and Thompson 1997). Both
#     models use the error contrasts of the model with the candidate,
#       RLR = 2 { max_delta l1(delta) - max_delta l0R(delta) },
#       l0R(delta) = -1/2 { L_F(delta) + log I_j(delta) + (nu - 1) log RSS_0(delta) },
#     with l1 the restricted profile of the model with the candidate. Its profiles are bounded in
#     delta. Size relative to alpha and power with the chi2_1 reference.
#  2. Size-adjusted power. For the Wald statistic t(delta_hat), the null-estimate statistic
#     Q(delta_hat_0), the likelihood ratio and the restricted likelihood ratio, the threshold is
#     the exact null quantile at the true ratio, and the power is evaluated at the same
#     noncentrality lambda = t_{nu-1, 1-alpha/2} as in 36.
# Each statistic is evaluated on a grid of the standardized effect u for every draw of r, and the
# rejection boundaries at any threshold are interpolated on that grid (one crossing per side of
# the zero of the score, as in 36). The row lrt_chi2_grid repeats the nominal likelihood ratio of
# 36 with this interpolation, as a check against the bisection of 36.
# Also written: the size of the likelihood ratio and of the restricted likelihood ratio referred
# to chi2_1 when the ratio is known (Student law of t(delta_0), no simulation).
# Usage (from the repository root): nice -n 10 Rscript calibration/56_restricted_lr_adjusted_power.R [cores=2] [--smoke]
# Output: results/rlr_adjusted/<config>.rds, results/rlr_adjusted.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # functions, GENES, load_geu, cfgs, GRID2, AL
a0 <- commandArgs(TRUE)
NCORES <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 2L
SMOKE <- "--smoke" %in% a0
B_DRAW <- 1000L
UG <- c(seq(-40, -13, by = 1), seq(-12, 12, by = 0.1), seq(13, 40, by = 1))

calib_rlr <- function(sp, delta0, tag, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)             # same draws as calib_config of 36
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1  <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+")
  rl1 <- apply(l1, 1, max); dh1 <- pick_grid(l1, GRID2)
  ml1 <- apply(-sweep(n * log(RSS_r), 2, ld2, "+"), 1, max)
  iD <- 1 / (1 + outer(ell, dh1))
  A1 <- colSums(Rk2 * iD) + Rp^2 + ssr; B1 <- colSums(xi * Rk * iD) + xp * Rp; S21 <- colSums(xi^2 * iD) + xp^2
  sig_h <- sqrt((A1 - B1^2 / S21) / (nu - 1))
  wa <- sqrt(S21 / S0) / sig_h; wg <- (B1 / S21) * sqrt(S21) / sig_h
  ev <- function(u, what = c("lrt", "rlr", "Q")) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    lA <- log(crossprod(sk2, W2) + (s2p^2 + ssr)); out <- list()
    if ("lrt" %in% what) out$lrt <- ml1 - apply(-sweep(n * lA, 2, ld2, "+"), 1, max)
    if ("rlr" %in% what) out$rlr <- 2 * (rl1 - apply(-0.5 * sweep((nu - 1) * lA, 2, ld2 + log(S2g), "+"), 1, max))
    if ("Q" %in% what) {
      dh0 <- pick_grid(-0.5 * sweep(nu * lA, 2, ld2, "+"), GRID2)
      iD0 <- 1 / (1 + outer(ell, dh0))
      Ah <- colSums(sk2 * iD0) + s2p^2 + ssr; Bh <- colSums(xi * sk * iD0) + xp * s2p
      out$Q <- Bh / sqrt((colSums(xi^2 * iD0) + xp^2) * Ah / nu)
    }
    out
  }
  lo <- rep(-30, B); hi <- rep(30, B)                        # zero of the score, as in 36
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < 0
    lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  f0 <- ev(split); f0$Q2 <- f0$Q^2
  nU <- length(UG); Fm <- list(lrt = matrix(0, B, nU), rlr = matrix(0, B, nU), Q2 = matrix(0, B, nU))
  for (k in seq_len(nU)) { e <- ev(rep(UG[k], B)); Fm$lrt[, k] <- e$lrt; Fm$rlr[, k] <- e$rlr; Fm$Q2[, k] <- e$Q^2 }
  right <- outer(split, UG, "<"); kR0 <- rowSums(!right) + 1L   # first grid point to the right of the zero
  cross <- function(key, c0) {                               # boundaries (dn, up) of {statistic >= c0} in u
    Fk <- Fm[[key]]; fz <- f0[[key]]
    M <- (Fk >= c0) & right; any_r <- rowSums(M) > 0; k <- max.col(M, ties.method = "first")
    k[!any_r] <- nU; first <- k <= kR0
    ua <- ifelse(first, split, UG[pmax(k - 1L, 1L)]); fa <- ifelse(first, fz, Fk[cbind(seq_len(B), pmax(k - 1L, 1L))])
    ub <- UG[k]; fb <- Fk[cbind(seq_len(B), k)]
    up <- ua + (c0 - fa) / pmax(fb - fa, 1e-300) * (ub - ua); up <- pmin(pmax(up, ua), ub); up[!any_r] <- Inf
    M <- (Fk >= c0) & !right; any_l <- rowSums(M) > 0; k <- max.col(M, ties.method = "last")
    k[!any_l] <- 1L; last <- k >= kR0 - 1L
    ua <- ifelse(last, split, UG[pmin(k + 1L, nU)]); fa <- ifelse(last, fz, Fk[cbind(seq_len(B), pmin(k + 1L, nU))])
    ub <- UG[k]; fb <- Fk[cbind(seq_len(B), k)]
    dn <- ua + (c0 - fa) / pmax(fb - fa, 1e-300) * (ub - ua); dn <- pmax(pmin(dn, ua), ub); dn[!any_l] <- -Inf
    always <- fz >= c0; up[always] <- split[always]; dn[always] <- split[always]
    list(up = up, dn = dn)
  }
  lin <- function(c0) list(up = (c0 - wg) / wa, dn = (-c0 - wg) / wa)
  rej <- function(rt, lam) mean(pnorm(lam - rt$up) + pnorm(rt$dn - lam))
  solve_c <- function(size_fun, alpha, hi) {                 # threshold with exact size alpha (size decreasing in c)
    lo <- 0
    for (it in 1:60) { mid <- (lo + hi) / 2; if (size_fun(mid) > alpha) lo <- mid else hi <- mid }
    (lo + hi) / 2
  }
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  res <- list()
  add <- function(method, k, rt, thr) res[[length(res) + 1]] <<- data.frame(tag, method = method, alpha = AL[k],
      size_ratio = rej(rt, 0) / AL[k], power = rej(rt, tq_[k]), lambda = tq_[k], threshold = thr)
  for (k in seq_along(AL)) {
    add("lrt_chi2_grid", k, cross("lrt", cq[k]), cq[k])
    add("rlr_chi2", k, cross("rlr", cq[k]), cq[k])
    cw <- solve_c(function(c0) rej(lin(c0), 0), AL[k], 2000); add("wald_adj", k, lin(cw), cw)
    for (key in c("Q2", "lrt", "rlr")) {
      ck <- solve_c(function(c0) rej(cross(key, c0), 0), AL[k], max(Fm[[key]]))
      add(paste0(c(Q2 = "score", lrt = "lrt", rlr = "rlr")[[key]], "_adj"), k, cross(key, ck), ck)
    }
    ## known ratio: t(delta_0) ~ t_{nu-1}; T_ML = n log{1 + t^2/(nu-1)}, T_R = (nu-1) log{1 + t^2/(nu-1)}
    for (m in c("lrt_known", "rlr_known")) {
      dd <- if (m == "lrt_known") n else nu - 1
      tc <- sqrt((nu - 1) * (exp(cq[k] / dd) - 1))
      res[[length(res) + 1]] <- data.frame(tag, method = m, alpha = AL[k],
        size_ratio = 2 * pt(-tc, nu - 1) / AL[k],
        power = pt(-tc, nu - 1, ncp = tq_[k]) + pt(tc, nu - 1, ncp = tq_[k], lower.tail = FALSE), lambda = tq_[k], threshold = cq[k])
    }
  }
  do.call(rbind, res)
}

OUT <- "calibration/results/rlr_adjusted"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (SMOKE) cfgs <- cfgs[c(1, 9, 161)]
run_one <- function(z) {                                    # designs, tags and seeds of run_one of 36
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "")
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L)
  } else {
    g <- GENES[z$i, ]; id <- sprintf("geu_%s_%d_h%.2f", g$symbol, z$pos, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973))
  }
  out <- calib_rlr(sp, z$h2 / (1 - z$h2), tag)
  if (SMOKE) return(out)
  saveRDS(out, ck); ck
}
if (SMOKE) {
  old <- read.csv("calibration/results/calib_enlarged.csv")
  for (z in cfgs) { t0 <- Sys.time(); out <- run_one(z)
    cat(sprintf("\n%s  nu_delta = %.1f  (%.0f s)\n", out$config[1], out$nu_delta[1], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(out[, c("method", "alpha", "size_ratio", "power", "threshold")], row.names = FALSE, digits = 4)
    print(old[old$config == out$config[1] & old$method %in% c("lrt_chi2", "wald_t", "score_t", "oracle_t"),
              c("method", "alpha", "size_ratio", "power")], row.names = FALSE, digits = 4) }
} else {
  cat("configurations:", length(cfgs), "\n")
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) paste("ERROR", conditionMessage(e))),
                 mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT, "rds$", full.names = TRUE)
  write.csv(do.call(rbind, lapply(fs, readRDS)), "calibration/results/rlr_adjusted.csv", row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
