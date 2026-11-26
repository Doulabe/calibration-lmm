# ==============================================================================
# 51_protected_refine.R
#
# Refinement of the protected thresholds of 49_protected_plugin_v2.R where the evaluation of
# 50_protected_plugin_eval.R (100,000 independent draws at delta0) exceeds the level, followed by
# a new evaluation. With 10,000 common draws the 2-se threshold at 5e-8 can rest on a few draws,
# and then lies below the exact one (50: Wald size 1.8 times the level for RPS26 56931934, h2 0.20).
#   flags      Wald (LR) when a wald_protected (lr_protected) row of 50 has a size ratio above 1
#   ratios     stage-2 ratios of 49 whose threshold is within 10% (Wald) or 5% (LR, at most 3) of the
#              stage-2 maximum at some level, and for the Wald the stage-1 neighbours of the maximizer
#              at 5e-8
#   stage 3    c+ (p_delta(c) + 2 se = alpha) recomputed at these ratios with common draws (seed +
#              5,000,000) by blocks of 100,000, adding blocks until the tail at the largest c+ of each
#              level rests on at least 100 effective draws or the cap is reached (Wald 1,000,000
#              draws, LR 300,000). Wald: coefficients of all draws kept, Brent's method as in 49.
#              LR: tails accumulated on a geometric grid of thresholds (0.8 to 2 times the stage-2
#              value, 41 points per level), c+ interpolated on the log of p + 2 se, u grid [-8, 8]
#              (size only, the measure beyond is below 1e-15).
#   evaluation size and power at delta0 with independent draws (seed + 6,000,000), 1,000,000 for
#              the Wald and 100,000 for the LR, at the stage-3 and stage-2 thresholds.
# --check compares the Wald coefficients, the LR table, the grid tails and the restricted u grid
# with the code of 49 on one chunk, and times them.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/51_protected_refine.R [cores=3] [--check] [--only=ID,...]
# Output: results/protected_refine/<config>.rds (wald_protected and/or lr_protected rows, which
#         replace those of 50 in the tables), results/protected_refine.csv
# ==============================================================================
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
CHECK  <- "--check" %in% args
ONLY   <- { a <- grep("^--only=", args, value = TRUE); if (length(a)) strsplit(sub("^--only=", "", a), ",")[[1]] else NULL }
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
CH2 <- 2500L; BLOCK <- 100000L; CAPW <- 1000000L; CAPL <- 300000L; EVW <- 1000000L; EVL <- 100000L; ESS_MIN <- 100
TINY <- "--tiny" %in% args                               # code path test with few draws, output in results/protected_refine_tiny
if (TINY) { BLOCK <- 5000L; CAPW <- 10000L; CAPL <- 5000L; EVW <- 5000L; EVL <- 5000L }
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))          # estimation grid of 36
AL    <- c(1e-3, 1e-5, 5e-8)
UG    <- seq(-10, 15, by = 0.025)                       # u grid for the LR (as 49)
UG3   <- seq(-8, 8, by = 0.025)                         # u grid of the stage-3 LR tails (size only)

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

## error contrasts at true ratio delta (as calib_config of 36)
contrasts <- function(sp, delta, Zk, zp, ssr) {
  d <- 1 + delta * sp$ell; S0 <- sum(sp$xi^2 / d) + sp$xp^2
  s <- Zk * sqrt(d)
  b <- (colSums((sp$xi / d) * s) + sp$xp * zp) / S0
  list(Rk = s - outer(sp$xi, b), Rp = zp - sp$xp * b, ssr = ssr, S0 = S0)
}
## Wald coefficients, REML and ML estimates, LR tabulated on UG (as 49; used by --check only)
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
## Wald coefficients only (pieces without the ML estimate and the LR table)
wald_pieces <- function(sp, C) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  Rk2 <- C$Rk * C$Rk
  A_r <- crossprod(Rk2, W2) + (C$Rp^2 + C$ssr); B_r <- crossprod(xi * C$Rk, W2) + xp * C$Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  dh1 <- pick_grid(-0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"), GRID2)
  iD <- 1 / (1 + outer(ell, dh1))
  A <- colSums(Rk2 * iD) + C$Rp^2 + C$ssr; Bx <- colSums(xi * C$Rk * iD) + xp * C$Rp
  S2 <- colSums(xi^2 * iD) + xp^2; sig_h <- sqrt((A - Bx^2 / S2) / (nu - 1))
  list(wa = sqrt(S2 / C$S0) / sig_h, wg = (Bx / S2) * sqrt(S2) / sig_h)
}
## LR tabulated on the u grid ug (pieces without the REML estimate and the Wald coefficients)
lr_table <- function(sp, C, ug) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- nu + 1L; B <- ncol(C$Rk)
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  Rk2 <- C$Rk * C$Rk
  A_r <- crossprod(Rk2, W2) + (C$Rp^2 + C$ssr); B_r <- crossprod(xi * C$Rk, W2) + xp * C$Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+")
  ml1 <- L1ml[cbind(seq_len(B), max.col(L1ml, "first"))]
  S2m <- matrix(S2g, B, length(GRID2), byrow = TRUE); ldm <- matrix(ld2, B, length(GRID2), byrow = TRUE)
  LRt <- matrix(0, B, length(ug))
  for (k in seq_along(ug)) {
    bv <- ug[k] / sqrt(C$S0)
    O <- -n * log(A_r + (2 * bv) * B_r + (bv * bv) * S2m) - ldm
    LRt[, k] <- ml1 - O[cbind(seq_len(B), max.col(O, "first"))]
  }
  LRt
}
## per-draw rejection probabilities with u ~ N(lam, 1) (as 49, with the u grid as an argument)
tw_draw <- function(c, lam, P) pnorm(lam - (c - P$wg) / P$wa) + pnorm((-c - P$wg) / P$wa - lam)
tl_draw <- function(c, lam, LRt, ug = UG) {
  K <- length(ug); z <- ug - lam
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
    zs <- ug[k] + (cb - v1) / (v2 - v1) * (ug[k + 1] - ug[k]) - lam
    m <- ifelse(L1[tr],
                ifelse(z[k] >= 0, PU[k] - pnorm(zs, lower.tail = FALSE), pnorm(zs) - PL[k]),
                ifelse(zs >= 0, pnorm(zs, lower.tail = FALSE) - PU[k + 1], PL[k + 1] - pnorm(zs)))
    rs <- rowsum(m, b); ix <- as.integer(rownames(rs)); tot[ix] <- tot[ix] + rs[, 1]
  }
  tot
}
## the same tails at every threshold of an increasing grid cg in one pass (draws x thresholds)
tl_grid <- function(cg, lam, LRt, ug) {
  K <- length(ug); M <- length(cg); B <- nrow(LRt); z <- ug - lam
  PU <- pnorm(z, lower.tail = FALSE); PL <- pnorm(z)
  cell <- ifelse(z[-K] >= 0, PU[-K] - PU[-1], PL[-1] - PL[-K])
  v1 <- LRt[, -K, drop = FALSE]; v2 <- LRt[, -1, drop = FALSE]
  jl <- matrix(findInterval(pmin(v1, v2), cg), B)       # thresholds 1..jl: whole cell in the set
  jh <- matrix(findInterval(pmax(v1, v2), cg), B)       # thresholds jl+1..jh: cell crossed
  Tm <- numeric(B * M)
  i <- which(jl > 0)                                    # whole cells, at index jl, then summed downwards
  if (length(i)) {
    b <- (i - 1L) %% B + 1L; k <- (i - 1L) %/% B + 1L
    rs <- rowsum(cell[k], b + B * (jl[i] - 1L)); Tm[as.integer(rownames(rs))] <- rs[, 1]
  }
  Tm <- matrix(Tm, B, M) %*% outer(seq_len(M), seq_len(M), ">=")
  i <- which(jh > jl)                                   # crossed cells, one term per threshold
  if (length(i)) {
    nr <- jh[i] - jl[i]
    b <- rep((i - 1L) %% B + 1L, nr); k <- rep((i - 1L) %/% B + 1L, nr); ii <- rep(jl[i], nr) + sequence(nr)
    c <- cg[ii]; a1 <- v1[cbind(b, k)]; a2 <- v2[cbind(b, k)]
    zs <- ug[k] + (c - a1) / (a2 - a1) * (ug[k + 1] - ug[k]) - lam
    m <- ifelse(a1 >= c,
                ifelse(z[k] >= 0, PU[k] - pnorm(zs, lower.tail = FALSE), pnorm(zs) - PL[k]),
                ifelse(zs >= 0, pnorm(zs, lower.tail = FALSE) - PU[k + 1], PL[k + 1] - pnorm(zs)))
    rs <- rowsum(m, b + B * (ii - 1L)); P <- numeric(B * M); P[as.integer(rownames(rs))] <- rs[, 1]
    Tm <- Tm + P
  }
  Tm + PL[1] * outer(findInterval(LRt[, 1], cg), seq_len(M), ">=") +
       PU[K] * outer(findInterval(LRt[, K], cg), seq_len(M), ">=")
}
## threshold where the (decreasing) tail estimate f crosses alpha: bracket from c0, then Brent (as 49)
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
mp  <- function(v) mean(v) + 2 * sd(v) / sqrt(length(v))
ess <- function(v) sum(v)^2 / sum(v^2)
## largest threshold where p + 2 se crosses alpha on the grid, log-linear interpolation
solve_grid <- function(cg, s1, s2, B, alpha) {
  f <- s1 / B + 2 * sqrt(pmax(s2 - s1^2 / B, 0) / (B - 1) / B)
  i <- max(c(0L, which(f >= alpha)))
  if (i == 0L || i == length(cg)) return(c(c = NA_real_, ess = NA_real_, edge = if (i == 0L) -1 else 1))
  w <- (log(alpha) - log(f[i])) / (log(f[i + 1]) - log(f[i]))
  c(c = cg[i] + w * (cg[i + 1] - cg[i]), ess = s1[i + 1]^2 / s2[i + 1], edge = 0)
}

## stage 3 and evaluation for one configuration and one statistic
refine_wald <- function(sp, delta0, seed, tq, Q2, S, D1, c2) {
  near <- function(Qk) S[Qk >= 0.9 * max(Qk, na.rm = TRUE) & !is.na(Qk)]
  R <- unique(unlist(lapply(1:3, function(k) near(Q2[, k]))))
  a3 <- S[which.max(Q2[, 3])]; j <- which.min(abs(log10(D1 + 1e-12) - log10(a3 + 1e-12)))
  R <- sort(unique(c(R, D1[pmax(1, j - 1):pmin(length(D1), j + 1)])))
  set.seed(seed + 5000000L)
  WA <- WG <- vector("list", length(R)); for (r in seq_along(R)) WA[[r]] <- WG[[r]] <- list()
  B <- 0L
  repeat {
    for (j in seq_len(BLOCK %/% CH2)) {
      Zk <- matrix(rnorm(sp$r * CH2), sp$r, CH2); zp <- if (sp$nperp) rnorm(CH2) else numeric(CH2)
      ssr <- if (sp$nrest > 0) rchisq(CH2, sp$nrest) else numeric(CH2)
      for (r in seq_along(R)) {
        P <- wald_pieces(sp, contrasts(sp, R[r], Zk, zp, ssr))
        WA[[r]][[length(WA[[r]]) + 1]] <- P$wa; WG[[r]][[length(WG[[r]]) + 1]] <- P$wg
      }
    }
    B <- B + BLOCK
    Pr <- lapply(seq_along(R), function(r) list(wa = unlist(WA[[r]]), wg = unlist(WG[[r]])))
    Cq <- sapply(seq_along(R), function(r) vapply(1:3, function(k)
      qsolve(function(c) mp(tw_draw(c, 0, Pr[[r]])), AL[k], c2[k]), 0))   # 3 x ratios
    c3 <- apply(Cq, 1, max, na.rm = TRUE); arg <- apply(Cq, 1, which.max)
    e3 <- vapply(1:3, function(k) ess(tw_draw(c3[k], 0, Pr[[arg[k]]])), 0)
    if (min(e3) >= ESS_MIN || B >= CAPW) break
  }
  rm(WA, WG, Pr)
  ## evaluation at delta0, independent draws
  set.seed(seed + 6000000L); acc <- array(0, c(2, 3, 5))
  for (j in seq_len(EVW %/% CH2)) {
    Zk <- matrix(rnorm(sp$r * CH2), sp$r, CH2); zp <- if (sp$nperp) rnorm(CH2) else numeric(CH2)
    ssr <- if (sp$nrest > 0) rchisq(CH2, sp$nrest) else numeric(CH2)
    P <- wald_pieces(sp, contrasts(sp, delta0, Zk, zp, ssr))
    for (k in 1:3) for (s in 1:2) {
      thr <- if (s == 1) c3[k] else c2[k]
      v0 <- tw_draw(thr, 0, P); v1 <- tw_draw(thr, tq[k], P)
      acc[s, k, ] <- acc[s, k, ] + c(length(v0), sum(v0), sum(v0^2), sum(v1), sum(v1^2))
    }
  }
  list(stat = "W", R = R, B = B, c3 = c3, arg = R[arg], ess3 = e3, Cq = Cq, acc = acc, c2 = c2)
}
refine_lr <- function(sp, delta0, seed, tq, Q2, S, c2) {
  near <- function(Qk) S[order(-Qk)][seq_len(min(3, sum(Qk >= 0.95 * max(Qk, na.rm = TRUE), na.rm = TRUE)))]
  R <- sort(unique(unlist(lapply(4:6, function(k) near(Q2[, k])))))
  R <- R[stats::complete.cases(Q2[match(R, S), 4:6, drop = FALSE])]
  if (length(R) > 3) R <- R[order(-Q2[match(R, S), 6])][1:3]           # at most 3, ranked at 5e-8
  R <- sort(R)
  G <- lapply(R, function(d) { q <- Q2[match(d, S), 4:6]
    lapply(1:3, function(k) q[k] * exp(seq(log(0.8), log(2), length.out = 41))) })
  cgs <- lapply(G, function(g) sort(unlist(g)))
  s1 <- lapply(cgs, function(g) numeric(length(g))); s2 <- s1
  set.seed(seed + 5000000L); B <- 0L
  repeat {
    for (j in seq_len(BLOCK %/% CH2)) {
      Zk <- matrix(rnorm(sp$r * CH2), sp$r, CH2); zp <- if (sp$nperp) rnorm(CH2) else numeric(CH2)
      ssr <- if (sp$nrest > 0) rchisq(CH2, sp$nrest) else numeric(CH2)
      for (r in seq_along(R)) {
        Tm <- tl_grid(cgs[[r]], 0, lr_table(sp, contrasts(sp, R[r], Zk, zp, ssr), UG3), UG3)
        s1[[r]] <- s1[[r]] + colSums(Tm); s2[[r]] <- s2[[r]] + colSums(Tm^2)
      }
    }
    B <- B + BLOCK
    sol <- lapply(seq_along(R), function(r) lapply(1:3, function(k) {
      ix <- match(G[[r]][[k]], cgs[[r]]); solve_grid(G[[r]][[k]], s1[[r]][ix], s2[[r]][ix], B, AL[k]) }))
    Cq <- sapply(seq_along(R), function(r) vapply(1:3, function(k) sol[[r]][[k]][["c"]], 0))
    Eq <- sapply(seq_along(R), function(r) vapply(1:3, function(k) sol[[r]][[k]][["ess"]], 0))
    Ed <- sapply(seq_along(R), function(r) vapply(1:3, function(k) sol[[r]][[k]][["edge"]], 0))
    Cq <- matrix(Cq, 3); Eq <- matrix(Eq, 3); Ed <- matrix(Ed, 3)
    c3 <- apply(Cq, 1, max, na.rm = TRUE); arg <- apply(Cq, 1, which.max)
    e3 <- Eq[cbind(1:3, arg)]
    if (min(e3) >= ESS_MIN || B >= CAPL) break
  }
  if (any(Ed != 0)) warning("LR threshold outside the grid at some ratio")
  ## evaluation at delta0, independent draws
  set.seed(seed + 6000000L); acc <- array(0, c(2, 3, 5))
  for (j in seq_len(EVL %/% CH2)) {
    Zk <- matrix(rnorm(sp$r * CH2), sp$r, CH2); zp <- if (sp$nperp) rnorm(CH2) else numeric(CH2)
    ssr <- if (sp$nrest > 0) rchisq(CH2, sp$nrest) else numeric(CH2)
    LRt <- lr_table(sp, contrasts(sp, delta0, Zk, zp, ssr), UG)
    for (k in 1:3) for (s in 1:2) {
      thr <- if (s == 1) c3[k] else c2[k]
      v0 <- tl_draw(thr, 0, LRt); v1 <- tl_draw(thr, tq[k], LRt)
      acc[s, k, ] <- acc[s, k, ] + c(length(v0), sum(v0), sum(v0^2), sum(v1), sum(v1^2))
    }
  }
  list(stat = "L", R = R, B = B, c3 = c3, arg = R[arg], ess3 = e3, Cq = Cq, edge = Ed, acc = acc, c2 = c2)
}
rows_of <- function(x, tag, tq) {
  m <- if (x$stat == "W") "wald_protected" else "lr_protected"
  do.call(rbind, lapply(1:3, function(k) {
    a <- x$acc[1, k, ]; b <- x$acc[2, k, ]; B <- a[1]
    data.frame(tag, method = m, alpha = AL[k],
      size_ratio = a[2] / B / AL[k], size_se = sqrt(max(a[3] - a[2]^2 / B, 0) / (B - 1) / B) / AL[k],
      size_ess = a[2]^2 / a[3], power = a[4] / B, power_se = sqrt(max(a[5] - a[4]^2 / B, 0) / (B - 1) / B),
      lambda = tq[k], threshold = x$c3[k], draws = B,
      threshold_stage2 = x$c2[k], size_ratio_stage2 = b[2] / B / AL[k],
      size_se_stage2 = sqrt(max(b[3] - b[2]^2 / B, 0) / (B - 1) / B) / AL[k],
      crn_draws = x$B, crn_ess = x$ess3[k], argmax = x$arg[k])
  }))
}

## ---- configurations of 36 (as 48, 49, 50) -----------------------------------------------------
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

IN2  <- "calibration/results/protected_plugin_v2"
IN3  <- "calibration/results/protected_plugin_eval"
OUT4 <- if (TINY) "calibration/results/protected_refine_tiny" else "calibration/results/protected_refine"
dir.create(OUT4, showWarnings = FALSE, recursive = TRUE)
## flags from the evaluation of 50 (or the configurations named by --only, both statistics)
flag <- lapply(cfgs, function(z) {
  if (!is.null(ONLY)) return(if (z$id %in% ONLY) c("W", "L") else character(0))
  f3 <- file.path(IN3, paste0(z$id, ".rds")); if (!file.exists(f3)) return(character(0))
  e <- readRDS(f3)
  c(if (any(e$size_ratio[e$method == "wald_protected"] > 1)) "W", if (any(e$size_ratio[e$method == "lr_protected"] > 1)) "L")
})
worst <- vapply(cfgs, function(z) { f3 <- file.path(IN3, paste0(z$id, ".rds")); if (!file.exists(f3)) return(0)
  e <- readRDS(f3); max(e$size_ratio[grepl("^(wald|lr)_protected$", e$method)]) }, 0)
keep <- lengths(flag) > 0
cfgs <- cfgs[keep][order(-worst[keep])]; flag <- flag[keep][order(-worst[keep])]

setup <- function(z) {
  if (z$kind == "syn") {
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    seed <- 4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L
  } else {
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    seed <- 9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973)
  }
  list(sp = sp, seed = seed)
}
run_one <- function(z, fl) {
  ck <- file.path(OUT4, paste0(z$id, ".rds")); if (file.exists(ck)) return(ck)
  prev <- readRDS(file.path(IN2, paste0(z$id, ".rds")))
  st <- setup(z); sp <- st$sp; delta0 <- z$h2 / (1 - z$h2)
  tag <- prev[1, c("config", "design", "gene", "pos", "role", "n", "blocks", "h2", "nu_delta")]
  stopifnot(abs(indices(sp, delta0)[["nu_delta"]] - tag$nu_delta) < 1e-8 * tag$nu_delta)
  tq <- qt(1 - AL / 2, df = sp$nu0 - 1)
  TH <- attr(prev, "thresholds"); PR <- attr(prev, "profiles")
  t0 <- Sys.time(); out <- list(); info <- list()
  if ("W" %in% fl) { x <- refine_wald(sp, delta0, st$seed, tq, PR$Q2, PR$S, PR$D1, TH$wald_protected)
                     out[[length(out) + 1]] <- rows_of(x, tag, tq); info$W <- x[c("R", "B", "c3", "arg", "ess3", "Cq", "c2")] }
  if ("L" %in% fl) { x <- refine_lr(sp, delta0, st$seed, tq, PR$Q2, PR$S, TH$lr_protected)
                     out[[length(out) + 1]] <- rows_of(x, tag, tq); info$L <- x[c("R", "B", "c3", "arg", "ess3", "Cq", "edge", "c2")] }
  out <- do.call(rbind, out); rownames(out) <- NULL
  attr(out, "refine") <- info; attr(out, "sec") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(out, ck); ck
}

if (CHECK) {
  z <- cfgs[[1]]; st <- setup(z); sp <- st$sp; delta0 <- z$h2 / (1 - z$h2)
  set.seed(1); Zk <- matrix(rnorm(sp$r * CH2), sp$r, CH2); zp <- if (sp$nperp) rnorm(CH2) else numeric(CH2)
  ssr <- if (sp$nrest > 0) rchisq(CH2, sp$nrest) else numeric(CH2)
  C <- contrasts(sp, delta0, Zk, zp, ssr)
  tm <- function(e) { t <- Sys.time(); v <- eval.parent(substitute(e)); list(v = v, s = as.numeric(difftime(Sys.time(), t, units = "secs"))) }
  P0 <- tm(pieces(sp, C)); PW <- tm(wald_pieces(sp, C)); L1 <- tm(lr_table(sp, C, UG)); L3 <- tm(lr_table(sp, C, UG3))
  cat(sprintf("%s  r = %d  seconds per chunk: pieces %.2f, wald_pieces %.2f, lr_table UG %.2f, UG3 %.2f\n",
              z$id, sp$r, P0$s, PW$s, L1$s, L3$s))
  cat("Wald coefficients, max abs difference:", max(abs(PW$v$wa - P0$v$wa)), max(abs(PW$v$wg - P0$v$wg)), "\n")
  cat("LR table, max abs difference:", max(abs(L1$v - P0$v$LRt)), "\n")
  cg <- sort(c(11, 22, 32) %o% exp(seq(log(0.8), log(2), length.out = 41)))
  for (lam in c(0, 2.5)) {
    TG <- tm(tl_grid(cg, lam, L1$v, UG)); TD <- tm(sapply(cg, function(c) tl_draw(c, lam, L1$v)))
    cat(sprintf("lambda %.1f: tl_grid vs tl_draw, max relative difference %.2e (%.2f s vs %.2f s for %d thresholds)\n",
                lam, max(abs(TG$v - TD$v) / pmax(TD$v, 1e-300)), TG$s, TD$s, length(cg)))
  }
  T3 <- tl_grid(cg, 0, L3$v, UG3); T1 <- tl_grid(cg, 0, L1$v, UG); m1 <- colMeans(T1)
  cat("u grid [-8, 8] vs [-10, 15] at lambda 0: max absolute difference per draw", signif(max(abs(T3 - T1)), 3),
      "; max relative difference of the mean tail where it exceeds 1e-10:",
      signif(max(abs(colMeans(T3)[m1 > 1e-10] / m1[m1 > 1e-10] - 1)), 3), "\n")
} else {
  cat("configurations flagged:", length(cfgs), "  Wald:", sum(vapply(flag, function(f) "W" %in% f, NA)),
      "  LR:", sum(vapply(flag, function(f) "L" %in% f, NA)), "\n")
  st <- mclapply(seq_along(cfgs), function(i) tryCatch(run_one(cfgs[[i]], flag[[i]]),
                 error = function(e) paste("ERROR", cfgs[[i]]$id, conditionMessage(e))),
                 mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT4, "rds$", full.names = TRUE)
  if (length(fs)) write.csv(do.call(rbind, lapply(fs, readRDS)), paste0(OUT4, ".csv"), row.names = FALSE)
  cat("configs done:", length(fs), "\n")
  if (TINY) for (f in fs) { x <- readRDS(f); print(x[, c("config", "method", "alpha", "size_ratio", "size_ess", "threshold",
                                                           "threshold_stage2", "crn_draws", "crn_ess", "argmax")], digits = 3)
                            cat("seconds:", attr(x, "sec"), "\n"); str(attr(x, "refine"), max.level = 2) }
}
