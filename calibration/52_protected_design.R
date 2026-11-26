# ==============================================================================
# 52_protected_design.R
#
# Protected critical values computed once per design, without the true ratio, then frozen and
# evaluated at several ratios with independent draws (in
# 49-51 the stage-2 set contained delta0, the refinement of 51 followed an observed excess, and the
# same variant received different thresholds in its two heritability configurations).
# A design is a candidate in its window (GEUVADIS) or a synthetic design; the 174 configurations
# of 36_calibration_enlarged.R contain 149 designs (18 GEUVADIS variants appear at both h2).
#   stage 1  quantiles on the grid of 49 (half decades, refined by 0.1 decade near each maximum),
#            B1 = 1,000 common draws
#   stage 2  2-se thresholds c+ (p_delta(c) + 2 se = alpha) with B2 = 10,000 common draws at the three
#            largest stage-1 values of each statistic and level (no delta0)
#   stage 3  for a statistic whose tail at the largest c+ of some level rests on fewer than 100
#            effective draws, the procedure of 51: common draws by blocks of 100,000 at the ratios
#            near the maximum until 100 effective draws carry the tail (caps 1,000,000 Wald, 300,000 LR);
#            then (28 October) any stage-2 ratio left out whose value exceeds the refined maximum of a level
#            is accepted if it rests on 100 effective draws and refined in a second pass otherwise (the
#            first version kept the 3 LR ratios of the 5e-8 order and missed the 1e-3 and 1e-5 maximizers)
#   evaluation  thresholds frozen; size (and power at the configurations) with independent draws at
#            the true ratio of each configuration of the design and at the maximizing ratios of each
#            statistic at 1e-5 and 5e-8 (worst case), by blocks (Wald 100,000, LR 25,000) until every
#            level has 100 effective draws or a size ratio below 0.1 (caps 1,000,000 and 200,000)
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/52_protected_design.R [cores=3] [--tiny] [--only=DESIGN,...]
# Output: results/protected_design/<design>.rds (list: thresholds, stages, evaluation rows),
#         results/protected_design_eval.csv, results/protected_design_thresholds.csv
# The published values come from a run on x86-64 Linux (R 4.4.3, reference BLAS, 30 processes,
# 10 h). The same seeds give other draws on arm64 macOS, since the signs of the eigenvectors
# returned by LAPACK differ between the platforms; the results then differ within Monte Carlo error.
# ==============================================================================
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
TINY   <- "--tiny" %in% args                             # code path test with few draws
ONLY   <- { a <- grep("^--only=", args, value = TRUE); if (length(a)) strsplit(sub("^--only=", "", a), ",")[[1]] else NULL }
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
B1 <- 1000L; B2 <- 10000L; CH2 <- 2500L; BLOCK <- 100000L; CAPW <- 1000000L; CAPL <- 300000L; ESS_MIN <- 100
BLW <- 100000L; BLL <- 25000L; EVCAPW <- 1000000L; EVCAPL <- 200000L; NEGL <- 0.1
if (TINY) { B1 <- 200L; B2 <- 1000L; BLOCK <- 5000L; CAPW <- 10000L; CAPL <- 5000L; BLW <- 5000L; BLL <- 2500L
            EVCAPW <- 10000L; EVCAPL <- 5000L }
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))          # estimation grid of 36
AL    <- c(1e-3, 1e-5, 5e-8)
DG    <- c(0, 10^seq(-4, 4, by = 0.5), 1e6, 1e10)       # coarse grid of true ratios (49)
UG    <- seq(-10, 15, by = 0.025)                       # u grid for the LR (49)
UG3   <- seq(-8, 8, by = 0.025)                         # u grid for LR sizes (51)

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

## draws (as 36, 48-51)
base_draws <- function(sp, B, seed) {
  set.seed(seed)
  Zk <- matrix(rnorm(sp$r * B), sp$r, B); zp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  list(Zk = Zk, zp = zp, ssr = ssr)
}
new_draws <- function(sp, B) {                            # continues the current stream (51)
  Zk <- matrix(rnorm(sp$r * B), sp$r, B); zp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  list(Zk = Zk, zp = zp, ssr = ssr)
}
contrasts <- function(sp, delta, Zk, zp, ssr) {
  d <- 1 + delta * sp$ell; S0 <- sum(sp$xi^2 / d) + sp$xp^2
  s <- Zk * sqrt(d)
  b <- (colSums((sp$xi / d) * s) + sp$xp * zp) / S0
  list(Rk = s - outer(sp$xi, b), Rp = zp - sp$xp * b, ssr = ssr, S0 = S0)
}
## Wald coefficients, REML and ML estimates, LR tabulated on UG (as 49)
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
pieces_all <- function(sp, delta, D) {                    # by chunks of CH2 draws (49)
  B <- ncol(D$Zk); idx <- split(seq_len(B), ceiling(seq_len(B) / CH2))
  lapply(idx, function(i) pieces(sp, contrasts(sp, delta, D$Zk[, i, drop = FALSE], D$zp[i], D$ssr[i])))
}
## Wald coefficients only and LR table only (51)
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
## per-draw rejection probabilities with u ~ N(lam, 1) (49, 51)
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
tl_grid <- function(cg, lam, LRt, ug) {                  # tl_draw at every threshold of cg in one pass (51)
  K <- length(ug); M <- length(cg); B <- nrow(LRt); z <- ug - lam
  PU <- pnorm(z, lower.tail = FALSE); PL <- pnorm(z)
  cell <- ifelse(z[-K] >= 0, PU[-K] - PU[-1], PL[-1] - PL[-K])
  v1 <- LRt[, -K, drop = FALSE]; v2 <- LRt[, -1, drop = FALSE]
  jl <- matrix(findInterval(pmin(v1, v2), cg), B)
  jh <- matrix(findInterval(pmax(v1, v2), cg), B)
  Tm <- numeric(B * M)
  i <- which(jl > 0)
  if (length(i)) {
    b <- (i - 1L) %% B + 1L; k <- (i - 1L) %/% B + 1L
    rs <- rowsum(cell[k], b + B * (jl[i] - 1L)); Tm[as.integer(rownames(rs))] <- rs[, 1]
  }
  Tm <- matrix(Tm, B, M) %*% outer(seq_len(M), seq_len(M), ">=")
  i <- which(jh > jl)
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
tail_v <- function(stat, c, lam, PL) {                    # per-draw values over the chunks (49)
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
qsolve <- function(f, alpha, c0) {                         # bracket from c0, then Brent (49)
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
ess <- function(v) if (sum(v) > 0) sum(v)^2 / sum(v^2) else 0
solve_grid <- function(cg, s1, s2, B, alpha) {             # largest crossing of p + 2 se on the grid (51)
  f <- s1 / B + 2 * sqrt(pmax(s2 - s1^2 / B, 0) / (B - 1) / B)
  i <- max(c(0L, which(f >= alpha)))
  if (i == 0L || i == length(cg)) return(c(c = NA_real_, ess = NA_real_, edge = if (i == 0L) -1 else 1))
  w <- (log(alpha) - log(f[i])) / (log(f[i + 1]) - log(f[i]))
  c(c = cg[i] + w * (cg[i + 1] - cg[i]), ess = s1[i + 1]^2 / s2[i + 1], edge = 0)
}
getq <- function(Q, row) t(sapply(Q, function(q) q[row, ]))
wmax <- function(Q) { Q[is.na(Q)] <- -Inf; apply(Q, 2, which.max) }     # row of the maximum, by column

## stage 1: plain and 2-se quantiles on the coarse grid and its refinement, B1 common draws (49)
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
## stage 2: 2-se thresholds and their effective numbers of draws at the three largest stage-1 values
stage2 <- function(sp, s1, seed) {
  QWp1 <- getq(s1$Q, "Wp"); QLp1 <- getq(s1$Q, "Lp")
  top3 <- function(Q) unlist(lapply(seq_len(ncol(Q)), function(k) s1$D[order(-Q[, k])[1:3]]))
  S <- sort(unique(c(top3(QWp1), top3(QLp1))))
  near <- function(d) which.min(abs(log10(s1$D + 1e-12) - log10(d + 1e-12)))
  D2 <- base_draws(sp, B2, seed)
  Q2 <- E2 <- matrix(NA_real_, length(S), 6, dimnames = list(NULL, c(paste0("Wp", 1:3), paste0("Lp", 1:3))))
  for (i in seq_along(S)) {
    PL <- pieces_all(sp, S[i], D2); j <- near(S[i])
    for (k in seq_along(AL)) {
      Q2[i, k] <- qsolve(function(c) tail_mp("W", c, 0, PL), AL[k], QWp1[j, k])
      Q2[i, 3 + k] <- qsolve(function(c) tail_mp("L", c, 0, PL), AL[k], QLp1[j, k])
      if (is.finite(Q2[i, k])) E2[i, k] <- ess(tail_v("W", Q2[i, k], 0, PL))
      if (is.finite(Q2[i, 3 + k])) E2[i, 3 + k] <- ess(tail_v("L", Q2[i, 3 + k], 0, PL))
    }
    rm(PL)
  }
  list(S = S, Q2 = Q2, E2 = E2)
}
## stage 3 (51): common draws by blocks at the ratios R until the tail at the largest c+ of each level
## rests on ESS_MIN effective draws or the cap is reached
s3_wald <- function(sp, R, seed, c0) {
  set.seed(seed); WA <- WG <- replicate(length(R), list(), simplify = FALSE); B <- 0L
  repeat {
    for (j in seq_len(BLOCK %/% CH2)) {
      dr <- new_draws(sp, CH2)
      for (r in seq_along(R)) {
        P <- wald_pieces(sp, contrasts(sp, R[r], dr$Zk, dr$zp, dr$ssr))
        WA[[r]][[length(WA[[r]]) + 1]] <- P$wa; WG[[r]][[length(WG[[r]]) + 1]] <- P$wg
      }
    }
    B <- B + BLOCK
    Pr <- lapply(seq_along(R), function(r) list(wa = unlist(WA[[r]]), wg = unlist(WG[[r]])))
    Cq <- matrix(sapply(seq_along(R), function(r) vapply(1:3, function(k)
      qsolve(function(c) mp(tw_draw(c, 0, Pr[[r]])), AL[k], c0[k]), 0)), 3)
    arg <- wmax(t(Cq)); c3 <- Cq[cbind(1:3, arg)]
    e3 <- vapply(1:3, function(k) if (is.finite(c3[k])) ess(tw_draw(c3[k], 0, Pr[[arg[k]]])) else 0, 0)
    if (min(e3) >= ESS_MIN || B >= CAPW) break
  }
  list(c = c3, arg = R[arg], ess = e3, B = B, R = R, Cq = Cq)
}
s3_lr <- function(sp, R, seed, Q2R) {
  G <- lapply(seq_along(R), function(r) lapply(1:3, function(k) Q2R[r, k] * exp(seq(log(0.8), log(2), length.out = 41))))
  cgs <- lapply(G, function(g) sort(unlist(g)))
  s1 <- lapply(cgs, function(g) numeric(length(g))); s2 <- s1
  set.seed(seed); B <- 0L
  repeat {
    for (j in seq_len(BLOCK %/% CH2)) {
      dr <- new_draws(sp, CH2)
      for (r in seq_along(R)) {
        Tm <- tl_grid(cgs[[r]], 0, lr_table(sp, contrasts(sp, R[r], dr$Zk, dr$zp, dr$ssr), UG3), UG3)
        s1[[r]] <- s1[[r]] + colSums(Tm); s2[[r]] <- s2[[r]] + colSums(Tm^2)
      }
    }
    B <- B + BLOCK
    sol <- lapply(seq_along(R), function(r) lapply(1:3, function(k) {
      ix <- match(G[[r]][[k]], cgs[[r]]); solve_grid(G[[r]][[k]], s1[[r]][ix], s2[[r]][ix], B, AL[k]) }))
    get <- function(f) matrix(sapply(seq_along(R), function(r) vapply(1:3, function(k) sol[[r]][[k]][[f]], 0)), 3)
    Cq <- get("c"); Eq <- get("ess"); Ed <- get("edge")
    arg <- wmax(t(Cq)); c3 <- Cq[cbind(1:3, arg)]; e3 <- Eq[cbind(1:3, arg)]; e3[!is.finite(e3)] <- 0
    if (min(e3) >= ESS_MIN || B >= CAPL) break
  }
  list(c = c3, arg = R[arg], ess = e3, B = B, R = R, Cq = Cq, edge = Ed)
}
wald_ratios <- function(Q, S, D1) {                         # within 10% of the maximum, and neighbours (51)
  R <- unique(unlist(lapply(1:3, function(k) S[!is.na(Q[, k]) & Q[, k] >= 0.9 * max(Q[, k], na.rm = TRUE)])))
  a3 <- S[wmax(Q)[3]]; j <- which.min(abs(log10(D1 + 1e-12) - log10(a3 + 1e-12)))
  sort(unique(c(R, D1[pmax(1, j - 1):pmin(length(D1), j + 1)])))
}
lr_ratios <- function(Q, S) {                               # within 5% of the maximum, at most 3 (51)
  near <- function(Qk) S[order(-Qk)][seq_len(min(3, sum(Qk >= 0.95 * max(Qk, na.rm = TRUE), na.rm = TRUE)))]
  R <- sort(unique(unlist(lapply(1:3, function(k) near(Q[, k])))))
  R <- R[stats::complete.cases(Q[match(R, S), , drop = FALSE])]
  if (length(R) > 3) R <- R[order(-Q[match(R, S), 3])][1:3]
  sort(R)
}
## completion (28 October): a stage-2 ratio left out of stage 3 whose value exceeds the refined maximum of a
## level is accepted when that value rests on ESS_MIN effective draws, and refined in a second pass otherwise
## (the 3 LR ratios of lr_ratios follow the 5e-8 order and could leave out the maximizers at 1e-3 and 1e-5)
complete_s3 <- function(sp, st, s2, x, fin, seed) {
  Q <- s2$Q2[, if (st == "W") 1:3 else 4:6, drop = FALSE]; E <- s2$E2[, if (st == "W") 1:3 else 4:6, drop = FALSE]
  out <- which(!s2$S %in% x$R)
  up <- Q[out, , drop = FALSE] > matrix(fin$c, length(out), 3, byrow = TRUE); up[is.na(up)] <- FALSE
  rel <- E[out, , drop = FALSE] >= ESS_MIN; rel[is.na(rel)] <- FALSE
  R2 <- s2$S[out[rowSums(up & !rel) > 0]]; x2 <- NULL
  if (length(R2)) {
    Q2R <- Q[match(R2, s2$S), , drop = FALSE]
    for (k in 1:3) Q2R[is.na(Q2R[, k]), k] <- fin$c[k]
    x2 <- if (st == "W") s3_wald(sp, R2, seed, fin$c) else s3_lr(sp, R2, seed, Q2R)
    for (k in 1:3) if (is.finite(x2$c[k]) && x2$c[k] > fin$c[k]) {
      fin$c[k] <- x2$c[k]; fin$arg[k] <- x2$arg[k]; fin$ess[k] <- x2$ess[k] }
  }
  for (k in 1:3) {
    a <- out[up[, k] & rel[, k]]
    if (length(a)) { m <- a[which.max(Q[a, k])]
      if (Q[m, k] > fin$c[k]) { fin$c[k] <- Q[m, k]; fin$arg[k] <- s2$S[m]; fin$ess[k] <- E[m, k] } }
  }
  fin$n2 <- length(R2); fin$B2 <- if (is.null(x2)) 0L else x2$B
  list(fin = fin, x2 = x2)
}
## evaluation with independent draws by blocks; stops when every level has ESS_MIN effective draws
## or a size ratio whose upper 2-se value is below NEGL
acc_rows <- function(acc, B, cst, tq, power, method) {
  data.frame(method = method, alpha = AL,
    size_ratio = acc[, 2] / B / AL, size_se = sqrt(pmax(acc[, 3] - acc[, 2]^2 / B, 0) / (B - 1) / B) / AL,
    size_ess = ifelse(acc[, 3] > 0, acc[, 2]^2 / acc[, 3], 0),
    power = if (power) acc[, 4] / B else NA_real_,
    power_se = if (power) sqrt(pmax(acc[, 5] - acc[, 4]^2 / B, 0) / (B - 1) / B) else NA_real_,
    lambda = tq, threshold = cst, draws = B)
}
enough <- function(acc, B) {
  r <- acc[, 2] / B / AL; se <- sqrt(pmax(acc[, 3] - acc[, 2]^2 / B, 0) / (B - 1) / B) / AL
  e <- ifelse(acc[, 3] > 0, acc[, 2]^2 / acc[, 3], 0)
  all(e >= ESS_MIN | r + 2 * se < NEGL)
}
eval_wald <- function(sp, d, seed, cst, tq, power) {
  set.seed(seed); acc <- matrix(0, 3, 5); nom <- matrix(0, 3, 5); B <- 0L
  repeat {
    for (j in seq_len(BLW %/% CH2)) {
      dr <- new_draws(sp, CH2); P <- wald_pieces(sp, contrasts(sp, d, dr$Zk, dr$zp, dr$ssr))
      for (k in 1:3) {
        v0 <- tw_draw(cst[k], 0, P); v1 <- if (power) tw_draw(cst[k], tq[k], P) else 0
        acc[k, ] <- acc[k, ] + c(length(v0), sum(v0), sum(v0^2), sum(v1), sum(v1^2))
        if (power) { w0 <- tw_draw(tq[k], 0, P); w1 <- tw_draw(tq[k], tq[k], P)
                     nom[k, ] <- nom[k, ] + c(length(w0), sum(w0), sum(w0^2), sum(w1), sum(w1^2)) }
      }
    }
    B <- B + BLW
    if (enough(acc, B) || B >= EVCAPW) break
  }
  rbind(acc_rows(acc, B, cst, tq, power, "wald_protected"), if (power) acc_rows(nom, B, tq, tq, TRUE, "wald_t"))
}
eval_lr <- function(sp, d, seed, cst, tq, power) {
  ug <- if (power) UG else UG3; cq <- qchisq(1 - AL, 1)
  set.seed(seed); acc <- matrix(0, 3, 5); nom <- matrix(0, 3, 5); B <- 0L
  repeat {
    for (j in seq_len(BLL %/% CH2)) {
      dr <- new_draws(sp, CH2); LRt <- lr_table(sp, contrasts(sp, d, dr$Zk, dr$zp, dr$ssr), ug)
      for (k in 1:3) {
        v0 <- tl_draw(cst[k], 0, LRt, ug); v1 <- if (power) tl_draw(cst[k], tq[k], LRt, ug) else 0
        acc[k, ] <- acc[k, ] + c(length(v0), sum(v0), sum(v0^2), sum(v1), sum(v1^2))
        if (power) { w0 <- tl_draw(cq[k], 0, LRt, ug); w1 <- tl_draw(cq[k], tq[k], LRt, ug)
                     nom[k, ] <- nom[k, ] + c(length(w0), sum(w0), sum(w0^2), sum(w1), sum(w1^2)) }
      }
    }
    B <- B + BLL
    if (enough(acc, B) || B >= EVCAPL) break
  }
  rbind(acc_rows(acc, B, cst, tq, power, "lr_protected"), if (power) acc_rows(nom, B, cq, tq, TRUE, "lrt_chi2"))
}

## ---- configurations of 36 (as 48-51) grouped into designs -------------------------------------
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
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2, unlinked = FALSE, role = "block1",
                                   id = sprintf("syn_n%d_b%d_h%.2f", n, nb, h2))
for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = 1000L, nb = 100L, h2 = h2, unlinked = FALSE,
                                                           role = "block1", id = sprintf("syn_n1000_b100_h%.2f", h2))
for (n in c(200L, 1000L)) for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = 20L, h2 = h2,
  unlinked = TRUE, role = "unlinked", id = sprintf("syn_n%d_b20_h%.2f_unlinked", n, h2))
dkey <- vapply(cfgs, function(z) if (z$kind == "geu") sprintf("geu_%s_%d", GENES$symbol[z$i], z$pos)
               else sprintf("syn_n%d_b%d%s", z$n, z$nb, if (z$unlinked) "_unlinked" else ""), "")
designs <- lapply(split(seq_along(cfgs), factor(dkey, unique(dkey))), function(ix) {
  z <- cfgs[[ix[1]]]
  seed <- if (z$kind == "geu") 100L * (z$pos %% 99991L) + z$i else 10000000L + 1000L * z$n + 10L * z$nb + as.integer(z$unlinked)
  c(z[intersect(c("kind", "i", "pos", "n", "nb", "unlinked"), names(z))], list(id = dkey[ix[1]], seed = seed,
    configs = lapply(ix, function(j) cfgs[[j]][c("id", "h2", "role")])))
})
if (!is.null(ONLY)) designs <- designs[vapply(designs, `[[`, "", "id") %in% ONLY]
## smallest nu_delta first (the costly designs), from the configurations of 49
f49 <- "calibration/results/protected_plugin_v2_thresholds.csv"
if (file.exists(f49) && is.null(ONLY)) {
  T49 <- read.csv(f49)
  nmin <- vapply(designs, function(dz) min(T49$nu_delta[T49$config %in% vapply(dz$configs, `[[`, "", "id")], Inf), 0)
  designs <- designs[order(nmin)]
}

OUT5 <- if (TINY) "calibration/results/protected_design_tiny" else "calibration/results/protected_design"
dir.create(OUT5, showWarnings = FALSE, recursive = TRUE)
run_design <- function(dz) {
  ck <- file.path(OUT5, paste0(dz$id, ".rds")); if (file.exists(ck)) return(ck)
  if (dz$kind == "syn") {
    set.seed(1000 * dz$n + dz$nb)
    Z <- gen_blocks(dz$n, dz$nb, 50L, 0.95)
    if (dz$unlinked) { set.seed(77 + dz$n); sp <- spectral_setup(Z, x = rnorm(dz$n)) } else sp <- spectral_setup(Z, j = 25L)
    base <- data.frame(design = "AR1_rho0.95", gene = NA, pos = NA, n = dz$n, blocks = dz$nb)
  } else {
    W <- load_geu(dz$i); j <- match(dz$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    base <- data.frame(design = "GEUVADIS_cis", gene = GENES$symbol[dz$i], pos = dz$pos, n = nrow(W$Z), blocks = NA)
  }
  nu <- sp$nu0; tq <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  t0 <- Sys.time()
  s1 <- stage1(sp, base_draws(sp, B1, dz$seed + 100000000L), tq, cq)
  s2 <- stage2(sp, s1, dz$seed + 200000000L)
  fin <- list(); s3 <- list()
  for (st in c("W", "L")) {
    Q <- s2$Q2[, if (st == "W") 1:3 else 4:6, drop = FALSE]; E <- s2$E2[, if (st == "W") 1:3 else 4:6, drop = FALSE]
    i2 <- wmax(Q); c2 <- Q[cbind(i2, 1:3)]; e2 <- E[cbind(i2, 1:3)]
    if (all(is.finite(e2)) && min(e2) >= ESS_MIN) {
      fin[[st]] <- list(c = c2, arg = s2$S[i2], ess = e2, B = B2, stage = 2L, n2 = 0L, B2 = 0L)
    } else {
      x <- if (st == "W") s3_wald(sp, wald_ratios(Q, s2$S, s1$D), dz$seed + 300000000L, c2)
           else { R <- lr_ratios(Q, s2$S); s3_lr(sp, R, dz$seed + 350000000L, Q[match(R, s2$S), , drop = FALSE]) }
      bad <- !is.finite(x$c); x$c[bad] <- c2[bad]; x$arg[bad] <- s2$S[i2][bad]
      fin[[st]] <- list(c = x$c, arg = x$arg, ess = x$ess, B = x$B, stage = 3L); s3[[st]] <- x
      y <- complete_s3(sp, st, s2, x, fin[[st]], dz$seed + if (st == "W") 370000000L else 380000000L)
      fin[[st]] <- y$fin; if (!is.null(y$x2)) s3[[paste0(st, "2")]] <- y$x2
    }
  }
  t1 <- Sys.time()
  ## evaluation at the true ratio of each configuration (size and power) and at the maximizing ratios (size)
  ev <- list(); d0s <- vapply(dz$configs, function(cf) cf$h2 / (1 - cf$h2), 0)
  for (cf in dz$configs) {
    d0 <- cf$h2 / (1 - cf$h2); k <- if (abs(cf$h2 - 0.05) < 1e-9) 1L else 2L
    tag <- cbind(data.frame(config = cf$id), base, data.frame(role = cf$role, h2 = cf$h2,
                 nu_delta = indices(sp, d0)[["nu_delta"]], ratio = d0, kind = "config"))
    ev[[length(ev) + 1]] <- cbind(tag, eval_wald(sp, d0, dz$seed + 400000000L + 10000000L * k, fin$W$c, tq, TRUE))
    ev[[length(ev) + 1]] <- cbind(tag, eval_lr(sp, d0, dz$seed + 600000000L + 10000000L * k, fin$L$c, tq, TRUE))
  }
  for (st in c("W", "L")) {
    a <- unique(fin[[st]]$arg[2:3])
    a <- a[vapply(a, function(x) all(abs(log10(x + 1e-12) - log10(d0s + 1e-12)) >= 0.05), NA)]
    a <- a[!duplicated(round(log10(a + 1e-12), 1))]
    for (m in seq_along(a)) {
      tag <- cbind(data.frame(config = NA_character_), base, data.frame(role = "maximizing", h2 = NA_real_,
                   nu_delta = indices(sp, a[m])[["nu_delta"]], ratio = a[m], kind = "argmax"))
      ev[[length(ev) + 1]] <- cbind(tag, if (st == "W") eval_wald(sp, a[m], dz$seed + 400000000L + 10000000L * (2L + m), fin$W$c, tq, FALSE)
                                         else eval_lr(sp, a[m], dz$seed + 600000000L + 10000000L * (2L + m), fin$L$c, tq, FALSE))
    }
  }
  t2 <- Sys.time()
  sec <- function(a, b) as.numeric(difftime(b, a, units = "secs"))
  thr <- data.frame(design = dz$id, alpha = AL, wald = fin$W$c, wald_arg = fin$W$arg, wald_ess = fin$W$ess, wald_B = fin$W$B,
                    wald_stage = fin$W$stage, wald_n2 = fin$W$n2, wald_B2 = fin$W$B2, lr = fin$L$c, lr_arg = fin$L$arg,
                    lr_ess = fin$L$ess, lr_B = fin$L$B, lr_stage = fin$L$stage, lr_n2 = fin$L$n2, lr_B2 = fin$L$B2,
                    sec_thresholds = sec(t0, t1), sec_eval = sec(t1, t2))
  out <- list(design = dz$id, thresholds = thr, eval = do.call(rbind, ev), stage1 = s1, stage2 = s2, stage3 = s3)
  saveRDS(out, ck)
  cat(sprintf("%s  %-28s %5.1f min\n", format(Sys.time(), "%d %H:%M"), dz$id, sec(t0, t2) / 60))
  ck
}

cat("designs:", length(designs), " configurations:", sum(lengths(lapply(designs, `[[`, "configs"))), "\n")
st <- mclapply(designs, function(dz) tryCatch(run_design(dz), error = function(e) paste("ERROR", dz$id, conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st)[grepl("ERROR", unlist(st))])
fs <- list.files(OUT5, "rds$", full.names = TRUE); L <- lapply(fs, readRDS)
write.csv(do.call(rbind, lapply(L, `[[`, "eval")), paste0(OUT5, "_eval.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(L, `[[`, "thresholds")), paste0(OUT5, "_thresholds.csv"), row.names = FALSE)
cat("designs done:", length(fs), "\n")
if (TINY) for (x in L) { print(x$thresholds, digits = 4, row.names = FALSE)
  print(x$eval[, c("config", "kind", "ratio", "method", "alpha", "size_ratio", "size_se", "size_ess", "power", "draws")],
        digits = 3, row.names = FALSE) }
