# ==============================================================================
# 49_protected_plugin_v2.R
#
# Second pass of 48_protected_plugin_calibration.R on the 174 configurations of
# 36_calibration_enlarged.R, with the plug-in evaluated as in the application and the
# protected threshold computed with more draws.
#   plug-in    reject when T >= c_alpha(delta_hat(r)), delta_hat the REML estimate of the model
#              with the candidate for both statistics (as 27c_geuvadis_application_cov.R; 48 used
#              the ML estimate for the LR). c_alpha without margin, interpolated on the coarse
#              grid of true ratios refined by 0.1 decade within 0.5 decade of each maximizer,
#              with the B1 = 1,000 common draws of 48 (stage 1).
#   protected  reject when T >= c* = max_delta c+_alpha(delta), c+ solving p_delta(c) + 2 se = alpha.
#              Stage 2 recomputes c+ with B2 = 10,000 common draws at delta0 and at the three
#              largest stage-1 values of each statistic and level. 48 used B1 = 1,000 draws, whose
#              threshold underestimates the quantile at 5e-8 when delta0 is near the maximizer.
#   evaluation size and power at delta0 with B3 = 10,000 independent draws, for the nominal,
#              plug-in and protected thresholds (and the stage-1 protected threshold of 48), with
#              Monte Carlo standard errors and the effective number of draws carrying each tail.
# Quantiles are solved by Brent's method on the log tail (48: bisection) and the LR is tabulated
# on u in [-10, 15] by 0.025 (48: [-15, 15] by 0.02).
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/49_protected_plugin_v2.R [cores=3] [--smoke] [--only=ID,...]
# Output: results/protected_plugin_v2/<config>.rds, results/protected_plugin_v2.csv,
#         results/protected_plugin_v2_thresholds.csv
# ==============================================================================
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
SMOKE  <- "--smoke" %in% args
ONLY   <- { a <- grep("^--only=", args, value = TRUE); if (length(a)) strsplit(sub("^--only=", "", a), ",")[[1]] else NULL }
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
B1 <- 1000L; B2 <- 10000L; B3 <- 10000L; CH2 <- 2500L
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))          # estimation grid of 36
AL    <- c(1e-3, 1e-5, 5e-8)
DG    <- c(0, 10^seq(-4, 4, by = 0.5), 1e6, 1e10)       # coarse grid of true ratios
UG    <- seq(-10, 15, by = 0.025)                       # u grid for the LR

pick_grid <- function(O, grid) {                          # as in 36
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
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L, 213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {                                 # as in 36
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, pos = pos)
}

## standard normal draws of the kernel and perpendicular parts and residual chi-square (as 36, 48)
base_draws <- function(sp, B, seed) {
  set.seed(seed)
  Zk <- matrix(rnorm(sp$r * B), sp$r, B); zp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  list(Zk = Zk, zp = zp, ssr = ssr)
}
## error contrasts at true ratio delta (as calib_config of 36)
contrasts <- function(sp, delta, Zk, zp, ssr) {
  d <- 1 + delta * sp$ell; S0 <- sum(sp$xi^2 / d) + sp$xp^2
  s <- Zk * sqrt(d)
  b <- (colSums((sp$xi / d) * s) + sp$xp * zp) / S0
  list(Rk = s - outer(sp$xi, b), Rp = zp - sp$xp * b, ssr = ssr, S0 = S0)
}
## Wald coefficients (t = wa u + wg given r), REML and ML estimates, LR tabulated on UG (as 48)
pieces <- function(sp, C) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- nu + 1L; B <- ncol(C$Rk)
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  Rk2 <- C$Rk * C$Rk
  A_r <- crossprod(Rk2, W2) + (C$Rp^2 + C$ssr); B_r <- crossprod(xi * C$Rk, W2) + xp * C$Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  dh1 <- pick_grid(-0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"), GRID2)
  L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+")
  ml1 <- L1ml[cbind(seq_len(B), max.col(L1ml, "first"))]; dh1ml <- pick_grid(L1ml, GRID2)
  iD <- 1 / (1 + outer(ell, dh1))
  A <- colSums(Rk2 * iD) + C$Rp^2 + C$ssr; Bx <- colSums(xi * C$Rk * iD) + xp * C$Rp
  S2 <- colSums(xi^2 * iD) + xp^2; sig_h <- sqrt((A - Bx^2 / S2) / (nu - 1))
  S2m <- matrix(S2g, B, length(GRID2), byrow = TRUE); ldm <- matrix(ld2, B, length(GRID2), byrow = TRUE)
  LRt <- matrix(0, B, length(UG))
  for (k in seq_along(UG)) {
    bv <- UG[k] / sqrt(C$S0)
    O <- -n * log(A_r + (2 * bv) * B_r + (bv * bv) * S2m) - ldm
    LRt[, k] <- ml1 - O[cbind(seq_len(B), max.col(O, "first"))]
  }
  list(dh1 = dh1, dh1ml = dh1ml, wa = sqrt(S2 / C$S0) / sig_h, wg = (Bx / S2) * sqrt(S2) / sig_h, LRt = LRt)
}
pieces_all <- function(sp, delta, D) {                    # by chunks of CH2 draws
  B <- ncol(D$Zk); idx <- split(seq_len(B), ceiling(seq_len(B) / CH2))
  lapply(idx, function(i) pieces(sp, contrasts(sp, delta, D$Zk[, i, drop = FALSE], D$zp[i], D$ssr[i])))
}
## per-draw rejection probabilities with u ~ N(lam, 1); c scalar or one threshold per draw
tw_draw <- function(c, lam, P) pnorm(lam - (c - P$wg) / P$wa) + pnorm((-c - P$wg) / P$wa - lam)
tl_draw <- function(c, lam, LRt) {
  K <- length(UG); z <- UG - lam
  PU <- pnorm(z, lower.tail = FALSE); PL <- pnorm(z)
  cell <- ifelse(z[-K] >= 0, PU[-K] - PU[-1], PL[-1] - PL[-K])
  G <- LRt >= c
  L1 <- G[, -K, drop = FALSE]; L2 <- G[, -1, drop = FALSE]
  tot <- drop((L1 & L2) %*% cell) + G[, 1] * PL[1] + G[, K] * PU[K]
  tr <- which(xor(L1, L2), arr.ind = TRUE)
  if (nrow(tr)) {
    b <- tr[, 1]; k <- tr[, 2]
    v1 <- LRt[cbind(b, k)]; v2 <- LRt[cbind(b, k + 1)]
    cb <- if (length(c) == 1) rep(c, length(b)) else c[b]
    zs <- UG[k] + (cb - v1) / (v2 - v1) * (UG[k + 1] - UG[k]) - lam
    m <- ifelse(L1[tr],
                ifelse(z[k] >= 0, PU[k] - pnorm(zs, lower.tail = FALSE), pnorm(zs) - PL[k]),
                ifelse(zs >= 0, pnorm(zs, lower.tail = FALSE) - PU[k + 1], PL[k + 1] - pnorm(zs)))
    rs <- rowsum(m, b); ix <- as.integer(rownames(rs)); tot[ix] <- tot[ix] + rs[, 1]
  }
  tot
}
tail_v <- function(stat, c, lam, PL) {                    # per-draw values over the chunks
  off <- 0L; v <- vector("list", length(PL))
  for (i in seq_along(PL)) {
    P <- PL[[i]]; nb <- length(P$wa)
    cc <- if (length(c) == 1L) c else c[off + seq_len(nb)]
    v[[i]] <- if (stat == "W") tw_draw(cc, lam, P) else tl_draw(cc, lam, P$LRt)
    off <- off + nb
  }
  unlist(v, use.names = FALSE)
}
tail_m  <- function(stat, c, lam, PL) mean(tail_v(stat, c, lam, PL))
tail_mp <- function(stat, c, lam, PL) { v <- tail_v(stat, c, lam, PL); mean(v) + 2 * sd(v) / sqrt(length(v)) }
## threshold where the (decreasing) tail estimate f crosses alpha: bracket from c0, then Brent
qsolve <- function(f, alpha, c0) {
  g <- function(c) log(max(f(c), 1e-300) / alpha)
  if (!is.finite(c0) || c0 <= 0) c0 <- 1
  g0 <- g(c0)
  if (g0 > 0) {
    lo <- c0; glo <- g0; hi <- 1.25 * c0; ghi <- g(hi)
    while (ghi > 0) { if (hi > 1e4) return(NA_real_); lo <- hi; glo <- ghi; hi <- 1.25 * hi; ghi <- g(hi) }
  } else {
    hi <- c0; ghi <- g0; lo <- c0 / 1.25; glo <- g(lo)
    while (glo <= 0) {
      hi <- lo; ghi <- glo; lo <- lo / 1.25
      if (lo < 1e-3) { lo <- 0; glo <- g(0); break }
      glo <- g(lo)
    }
    if (glo <= 0) return(NA_real_)
  }
  uniroot(g, c(lo, hi), f.lower = glo, f.upper = ghi, tol = 1e-7 * hi + 1e-9)$root
}
interp_c <- function(dg, Q, dh) {                         # c_alpha(delta_hat), log scale in delta
  sapply(seq_len(ncol(Q)), function(k) {
    out <- rep(Q[1, k], length(dh)); pos <- dh > 0
    out[pos] <- approx(log10(dg[-1]), Q[-1, k], log10(dh[pos]), rule = 2)$y
    out })
}

## stage 1: plain and 2-se quantiles on the coarse grid and its refinement, B1 common draws
stage1 <- function(sp, D, tq, cq) {
  q_at <- function(d) {
    PL <- pieces_all(sp, d, D)
    out <- matrix(NA_real_, 4, length(AL), dimnames = list(c("W", "L", "Wp", "Lp"), NULL))
    for (k in seq_along(AL)) {
      out["W", k]  <- qsolve(function(c) tail_m("W", c, 0, PL), AL[k], tq[k])
      out["L", k]  <- qsolve(function(c) tail_m("L", c, 0, PL), AL[k], cq[k])
      out["Wp", k] <- qsolve(function(c) tail_mp("W", c, 0, PL), AL[k], out["W", k])
      out["Lp", k] <- qsolve(function(c) tail_mp("L", c, 0, PL), AL[k], out["L", k])
    }
    out
  }
  Qc <- lapply(DG, q_at)
  arg <- unlist(lapply(c("W", "L", "Wp", "Lp"), function(row)
    sapply(seq_along(AL), function(k) DG[which.max(sapply(Qc, function(q) q[row, k]))])))
  lr <- unique(round(unlist(lapply(log10(arg[arg > 0 & arg <= 1e4]), function(a) a + seq(-0.5, 0.5, by = 0.1))), 6))
  lr <- lr[!lr %in% round(log10(DG[DG > 0]), 6)]
  Qr <- lapply(10^lr, q_at)
  allD <- c(DG, 10^lr); o <- order(allD)
  list(D = allD[o], Q = c(Qc, Qr)[o])
}
getq <- function(Q, row) t(sapply(Q, function(q) q[row, ]))

calib_v2 <- function(sp, delta0, seed_calib, tag) {
  nu <- sp$nu0; tq <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  t0 <- Sys.time()
  s1 <- stage1(sp, base_draws(sp, B1, seed_calib + 1L), tq, cq)      # the draws of 48
  QW1 <- getq(s1$Q, "W"); QL1 <- getq(s1$Q, "L"); QWp1 <- getq(s1$Q, "Wp"); QLp1 <- getq(s1$Q, "Lp")
  cW1 <- apply(QWp1, 2, max, na.rm = TRUE); cL1 <- apply(QLp1, 2, max, na.rm = TRUE)
  aW1 <- s1$D[apply(QWp1, 2, which.max)]; aL1 <- s1$D[apply(QLp1, 2, which.max)]
  t1 <- Sys.time()
  ## stage 2: delta0 and the three largest stage-1 values of each statistic and level, B2 draws
  top3 <- function(Q) unlist(lapply(seq_len(ncol(Q)), function(k) s1$D[order(-Q[, k])[1:3]]))
  S <- sort(unique(c(delta0, top3(QWp1), top3(QLp1))))
  near <- function(d) which.min(abs(log10(s1$D + 1e-12) - log10(d + 1e-12)))
  D2 <- base_draws(sp, B2, seed_calib + 2000000L)
  Q2 <- matrix(NA_real_, length(S), 2 * length(AL), dimnames = list(NULL, c(paste0("Wp", 1:3), paste0("Lp", 1:3))))
  P0W <- P0L <- rep(NA_real_, length(AL))
  for (i in seq_along(S)) {
    PL <- pieces_all(sp, S[i], D2); j <- near(S[i])
    for (k in seq_along(AL)) {
      Q2[i, k] <- qsolve(function(c) tail_mp("W", c, 0, PL), AL[k], QWp1[j, k])
      Q2[i, 3 + k] <- qsolve(function(c) tail_mp("L", c, 0, PL), AL[k], QLp1[j, k])
    }
    if (S[i] == delta0) for (k in seq_along(AL)) {
      P0W[k] <- qsolve(function(c) tail_m("W", c, 0, PL), AL[k], Q2[i, k])
      P0L[k] <- qsolve(function(c) tail_m("L", c, 0, PL), AL[k], Q2[i, 3 + k])
    }
    rm(PL)
  }
  rm(D2)
  cW2 <- apply(Q2[, 1:3, drop = FALSE], 2, max, na.rm = TRUE); cL2 <- apply(Q2[, 4:6, drop = FALSE], 2, max, na.rm = TRUE)
  aW2 <- S[apply(Q2[, 1:3, drop = FALSE], 2, which.max)]; aL2 <- S[apply(Q2[, 4:6, drop = FALSE], 2, which.max)]
  t2 <- Sys.time()
  ## evaluation at delta0 with B3 independent draws
  PL0 <- pieces_all(sp, delta0, base_draws(sp, B3, seed_calib + 3000000L))
  dh <- unlist(lapply(PL0, `[[`, "dh1"), use.names = FALSE)
  pW <- interp_c(s1$D, QW1, dh); pL <- interp_c(s1$D, QL1, dh)
  res <- list()
  for (k in seq_along(AL)) {
    add <- function(m, stat, thr) {
      v0 <- tail_v(stat, thr, 0, PL0); v1 <- tail_v(stat, thr, tq[k], PL0)
      res[[length(res) + 1]] <<- data.frame(tag, method = m, alpha = AL[k],
        size_ratio = mean(v0) / AL[k], size_se = sd(v0) / sqrt(length(v0)) / AL[k],
        size_ess = sum(v0)^2 / sum(v0^2), power = mean(v1), power_se = sd(v1) / sqrt(length(v1)),
        lambda = tq[k], threshold = if (length(thr) == 1L) thr else median(thr))
    }
    add("wald_t", "W", tq[k]); add("lrt_chi2", "L", cq[k])
    add("wald_plugin", "W", pW[, k]); add("lr_plugin", "L", pL[, k])
    add("wald_protected", "W", cW2[k]); add("lr_protected", "L", cL2[k])
    add("wald_protected_b1", "W", cW1[k]); add("lr_protected_b1", "L", cL1[k])
  }
  t3 <- Sys.time()
  out <- do.call(rbind, res)
  sec <- function(a, b) as.numeric(difftime(b, a, units = "secs"))
  attr(out, "thresholds") <- data.frame(tag, alpha = AL, t_nominal = tq, chi2_nominal = cq,
    wald_b1 = cW1, wald_argmax_b1 = aW1, lr_b1 = cL1, lr_argmax_b1 = aL1,
    wald_protected = cW2, wald_argmax = aW2, lr_protected = cL2, lr_argmax = aL2,
    wald_delta0 = P0W, lr_delta0 = P0L, n_delta1 = length(s1$D), n_delta2 = length(S),
    sec_stage1 = sec(t0, t1), sec_stage2 = sec(t1, t2), sec_eval = sec(t2, t3))
  attr(out, "profiles") <- list(D1 = s1$D, Q1 = s1$Q, S = S, Q2 = Q2, dh_quantiles = quantile(dh, c(0, .1, .25, .5, .75, .9, 1)),
                                dh_zero = mean(dh == 0))
  out
}

## ---- configurations of 36 (as 48) -------------------------------------------------------------
NUS <- read.csv("calibration/results/nu_delta_per_snp_geuvadis.csv")
QS  <- c(0.02, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.98)
cfgs <- list()
for (i in seq_len(nrow(GENES))) for (h2 in c(0.05, 0.20)) {
  d <- NUS[NUS$gene == GENES$symbol[i] & abs(NUS$h2 - h2) < 1e-9, ]
  lv <- log(d$nu_delta); target <- quantile(lv, QS)
  picks <- integer(0)
  for (q in target) { o <- order(abs(lv - q)); o <- o[!d$pos[o] %in% c(GENES$lead_pos[i], d$pos[picks])]; picks <- c(picks, o[1]) }
  cand <- c(GENES$lead_pos[i], d$pos[picks])
  for (k in seq_along(cand)) cfgs[[length(cfgs) + 1]] <- list(kind = "geu", i = i, h2 = h2, pos = cand[k],
    role = if (k == 1) "lead" else sprintf("q%02d", round(100 * QS[k - 1])),
    id = sprintf("geu_%s_%d_h%.2f", GENES$symbol[i], cand[k], h2))
}
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2, unlinked = FALSE,
                                   id = sprintf("syn_n%d_b%d_h%.2f", n, nb, h2))
for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = 1000L, nb = 100L, h2 = h2, unlinked = FALSE,
                                                           id = sprintf("syn_n1000_b100_h%.2f", h2))
for (n in c(200L, 1000L)) for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = 20L, h2 = h2,
  unlinked = TRUE, id = sprintf("syn_n%d_b20_h%.2f_unlinked", n, h2))
if (!is.null(ONLY)) cfgs <- cfgs[vapply(cfgs, `[[`, "", "id") %in% ONLY]

OUT2 <- "calibration/results/protected_plugin_v2"
dir.create(OUT2, showWarnings = FALSE, recursive = TRUE)
run_one <- function(z) {
  ck <- file.path(OUT2, paste0(z$id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
  if (z$kind == "syn") {
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = z$id, design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    seed <- 4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L
  } else {
    g <- GENES[z$i, ]
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = z$id, design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    seed <- 9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973)
  }
  out <- calib_v2(sp, z$h2 / (1 - z$h2), seed, tag)
  if (SMOKE) return(out)
  saveRDS(out, ck); ck
}
if (SMOKE) {
  for (z in cfgs) {
    out <- run_one(z)
    cat(sprintf("\n%s  nu_delta = %.1f\n", out$config[1], out$nu_delta[1]))
    print(out[, c("method", "alpha", "size_ratio", "size_se", "size_ess", "power", "threshold")], digits = 3, row.names = FALSE)
    print(attr(out, "thresholds")[, c("alpha", "wald_b1", "wald_protected", "wald_delta0", "wald_argmax", "lr_b1", "lr_protected",
                                      "lr_delta0", "lr_argmax", "n_delta1", "n_delta2", "sec_stage1", "sec_stage2", "sec_eval")],
          digits = 4, row.names = FALSE)
    pr <- attr(out, "profiles"); cat("stage-2 ratios:", signif(pr$S, 3), "\n"); print(round(pr$dh_quantiles, 4))
  }
} else {
  cat("configurations:", length(cfgs), "\n")
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) paste("ERROR", z$id, conditionMessage(e))),
                 mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT2, "rds$", full.names = TRUE); L <- lapply(fs, readRDS)
  write.csv(do.call(rbind, L), "calibration/results/protected_plugin_v2.csv", row.names = FALSE)
  write.csv(do.call(rbind, lapply(L, attr, "thresholds")), "calibration/results/protected_plugin_v2_thresholds.csv",
            row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
