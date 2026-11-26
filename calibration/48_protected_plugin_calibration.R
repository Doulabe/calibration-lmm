# ==============================================================================
# 48_protected_plugin_calibration.R
#
# Size and power of the protected and plug-in tests built on the Wald and likelihood-ratio
# statistics, on the 174 configurations of 36_calibration_enlarged.R (rows of Table 2).
#   protected  reject when T >= c* = max_delta c+_alpha(delta), where c+_alpha(delta) solves
#              p_delta(c) + 2 se = alpha (the margin of the numerical protected p-value of the
#              article). The maximum is taken on a coarse grid of true ratios, refined by
#              0.1 decade within 0.5 decade of each maximizer. c* depends on the design only.
#   plug-in    reject when T >= c_alpha(delta_hat(r)), c_alpha without margin on the coarse
#              grid, delta_hat the estimate with the candidate (REML for Wald, ML for LR).
# The LR is tabulated on a grid of u (the null residual sum of squares is quadratic in u at
# each delta of the estimation grid), validated in 47_protected_plugin_pilot.R against the
# bisection of 36. Calibration draws at delta0 are those of 36, so the nominal Wald and LR rows
# reproduce results/calib_enlarged.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/48_protected_plugin_calibration.R [cores=3] [--smoke]
# Output: results/protected_plugin/<config>.rds, results/protected_plugin.csv,
#         results/protected_plugin_thresholds.csv
# ==============================================================================
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
SMOKE  <- "--smoke" %in% args
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
B_DRAW <- 1000L
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))          # estimation grid of 36
AL    <- c(1e-3, 1e-5, 5e-8)
DG    <- c(0, 10^seq(-4, 4, by = 0.5), 1e6, 1e10)       # coarse grid of true ratios
UG    <- seq(-15, 15, by = 0.02)                        # u grid for the LR

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

## error contrasts at true ratio delta from standard normals Zk, zp and ssr (as calib_config of 36)
contrasts <- function(sp, delta, Zk, zp, ssr) {
  d <- 1 + delta * sp$ell; S0 <- sum(sp$xi^2 / d) + sp$xp^2
  s <- Zk * sqrt(d)
  b <- (colSums((sp$xi / d) * s) + sp$xp * zp) / S0
  list(Rk = s - outer(sp$xi, b), Rp = zp - sp$xp * b, ssr = ssr, S0 = S0)
}
## Wald coefficients (t = wa u + wg given r) and LR tabulated on UG
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
## rejection probabilities with u ~ N(lam, 1), c scalar or one threshold per draw;
## with se = TRUE, returns the mean and its Monte Carlo standard error
tail_w <- function(c, lam, P, se = FALSE) {
  p <- pnorm(lam - (c - P$wg) / P$wa) + pnorm((-c - P$wg) / P$wa - lam)
  if (se) c(mean(p), sd(p) / sqrt(length(p))) else mean(p)
}
tail_lr <- function(c, lam, LRt, se = FALSE) {
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
  if (se) c(mean(tot), sd(tot) / sqrt(length(tot))) else mean(tot)
}
## threshold where the (decreasing) tail estimate f crosses alpha
qsolve <- function(f, alpha, hi, it = 24L) {
  while (f(hi) > alpha && hi < 1e4) hi <- 2 * hi
  if (f(hi) > alpha) return(NA_real_)
  lo <- 0
  for (i in seq_len(it)) { mid <- (lo + hi) / 2; if (f(mid) > alpha) lo <- mid else hi <- mid }
  hi
}
interp_c <- function(dg, Q, dh) {                         # c_alpha(delta_hat), log scale in delta
  sapply(seq_len(ncol(Q)), function(k) {
    out <- rep(Q[1, k], length(dh)); pos <- dh > 0
    out[pos] <- approx(log10(dg[-1]), Q[-1, k], log10(dh[pos]), rule = 2)$y
    out })
}
## quantiles at one true ratio: no-margin (plug-in) and 2-se margin (protected)
quantiles_at <- function(sp, delta, Zk, zp, ssr, tq, cq, plain = TRUE) {
  P <- pieces(sp, contrasts(sp, delta, Zk, zp, ssr))
  out <- matrix(NA_real_, 4, length(AL), dimnames = list(c("W", "L", "Wp", "Lp"), NULL))
  for (k in seq_along(AL)) {
    out["Wp", k] <- qsolve(function(c) sum(tail_w(c, 0, P, TRUE) * c(1, 2)), AL[k], 3 * tq[k])
    out["Lp", k] <- qsolve(function(c) sum(tail_lr(c, 0, P$LRt, TRUE) * c(1, 2)), AL[k], 3 * cq[k])
    if (plain) {
      out["W", k] <- qsolve(function(c) tail_w(c, 0, P), AL[k], 3 * tq[k])
      out["L", k] <- qsolve(function(c) tail_lr(c, 0, P$LRt), AL[k], 3 * cq[k])
    }
  }
  out
}

calib_protected <- function(sp, delta0, seed_calib, tag) {
  B <- B_DRAW; nu <- sp$nu0
  tq <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  t0 <- Sys.time()
  set.seed(seed_calib + 1L)                               # common random numbers over true ratios
  Zk <- matrix(rnorm(sp$r * B), sp$r, B); zp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  Qc <- lapply(DG, function(d) quantiles_at(sp, d, Zk, zp, ssr, tq, cq))
  ## refinement by 0.1 decade within 0.5 decade of each maximizer of the protected quantiles
  arg <- c(sapply(seq_along(AL), function(k) DG[which.max(sapply(Qc, function(q) q["Wp", k]))]),
           sapply(seq_along(AL), function(k) DG[which.max(sapply(Qc, function(q) q["Lp", k]))]))
  ref <- unique(unlist(lapply(arg[arg > 0 & arg <= 1e4], function(a) 10^(log10(a) + seq(-0.5, 0.5, by = 0.1)))))
  ref <- ref[!sapply(ref, function(x) any(abs(log10(x) - log10(DG[DG > 0])) < 1e-9))]
  Qr <- lapply(ref, function(d) quantiles_at(sp, d, Zk, zp, ssr, tq, cq, plain = FALSE))
  allD <- c(DG, ref); allQ <- c(Qc, Qr)
  getq <- function(row) t(sapply(allQ, function(q) q[row, ]))
  QWp <- getq("Wp"); QLp <- getq("Lp")
  cW <- apply(QWp, 2, max, na.rm = TRUE); cL <- apply(QLp, 2, max, na.rm = TRUE)
  aW <- allD[apply(QWp, 2, which.max)]; aL <- allD[apply(QLp, 2, which.max)]
  QW <- t(sapply(Qc, function(q) q["W", ])); QL <- t(sapply(Qc, function(q) q["L", ]))
  ## calibration draws of 36 at delta0
  set.seed(seed_calib)
  Zk0 <- matrix(rnorm(sp$r * B), sp$r, B); zp0 <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr0 <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  P0 <- pieces(sp, contrasts(sp, delta0, Zk0, zp0, ssr0))
  pW <- interp_c(DG, QW, P0$dh1); pL <- interp_c(DG, QL, P0$dh1ml)
  res <- list()
  for (k in seq_along(AL)) {
    add <- function(m, thr, ct) {
      sz <- if (grepl("^lr", m)) tail_lr(ct, 0, P0$LRt) else tail_w(ct, 0, P0)
      pw <- if (grepl("^lr", m)) tail_lr(ct, tq[k], P0$LRt) else tail_w(ct, tq[k], P0)
      res[[length(res) + 1]] <<- data.frame(tag, method = m, alpha = AL[k], size_ratio = sz / AL[k],
                                            power = pw, lambda = tq[k], threshold = thr)
    }
    add("wald_t", tq[k], tq[k]); add("lrt_chi2", cq[k], cq[k])
    add("wald_protected", cW[k], cW[k]); add("lr_protected", cL[k], cL[k])
    add("wald_plugin", median(pW[, k]), pW[, k]); add("lr_plugin", median(pL[, k]), pL[, k])
  }
  out <- do.call(rbind, res)
  attr(out, "thresholds") <- data.frame(tag, alpha = AL, t_nominal = tq, chi2_nominal = cq,
    wald_protected = cW, wald_argmax = aW, lr_protected = cL, lr_argmax = aL,
    n_delta = length(allD), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  out
}

## ---- configurations of 36 --------------------------------------------------------------------
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
    role = if (k == 1) "lead" else sprintf("q%02d", round(100 * QS[k - 1])))
}
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2, unlinked = FALSE)
for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = 1000L, nb = 100L, h2 = h2, unlinked = FALSE)
for (n in c(200L, 1000L)) for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = 20L, h2 = h2, unlinked = TRUE)
if (SMOKE) cfgs <- cfgs[1]

OUT <- "calibration/results/protected_plugin"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "")
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    seed <- 4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L
  } else {
    g <- GENES[z$i, ]; id <- sprintf("geu_%s_%d_h%.2f", g$symbol, z$pos, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    seed <- 9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973)
  }
  out <- calib_protected(sp, z$h2 / (1 - z$h2), seed, tag)
  if (SMOKE) return(out)
  saveRDS(out, ck); ck
}
if (SMOKE) {
  t0 <- Sys.time(); out <- run_one(cfgs[[1]])
  cat(sprintf("%s  nu_delta = %.1f  (%.0f s)\n", out$config[1], out$nu_delta[1],
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  print(out[, c("method", "alpha", "size_ratio", "power", "threshold")], digits = 3, row.names = FALSE)
  print(attr(out, "thresholds")[, c("alpha", "wald_protected", "wald_argmax", "lr_protected", "lr_argmax", "n_delta")], digits = 4)
} else {
  cat("configurations:", length(cfgs), "\n")
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) paste("ERROR", conditionMessage(e))),
                 mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT, "rds$", full.names = TRUE); L <- lapply(fs, readRDS)
  write.csv(do.call(rbind, L), "calibration/results/protected_plugin.csv", row.names = FALSE)
  write.csv(do.call(rbind, lapply(L, attr, "thresholds")), "calibration/results/protected_plugin_thresholds.csv",
            row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
