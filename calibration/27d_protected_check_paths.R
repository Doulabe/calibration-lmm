# ==============================================================================
# 27d_protected_check_paths.R
# Numerical control of the protected exact-law p-values of the likelihood-ratio path of
# the article on the candidate-specific kernel (27c_geuvadis_application_cov.R), for the
# two second-step selections of the Yoruba sample whose protected value was above the
# per-test level on the coarse grid of 27c (PEX6-region 43011766 after 43186801;
# FLVCR1-AS1 213326573 after 213015222). As in 33g_protected_check.R:
#   - the statistic is the one of 27c (ML likelihood ratio, delta_hat on the grid of 27c,
#     [0, 1e4], since the ML profile is unbounded, see 33h);
#   - stage 1: p(delta0) at 0 and 10^-3, ..., 10^10 by 0.25 decade (B1 draws);
#   - stage 2: 0.05-decade grid within 0.3 decade of the stage-1 maximum (B2 draws);
#   - Monte Carlo standard errors, controlled value = largest p_hat + 2 se.
# Functions copied from 27c (spectral_step, pick_grid, stats_from_coords); the
# likelihood-ratio part of exact_logp rewritten to return standard errors.
# Usage (from the repository root): nice -n 10 Rscript calibration/27d_protected_check_paths.R [cores=2] [--domain=1e3]
# --domain=1e3: statistic and distribution on the decision domain of the paths, [0, 1e3] with 41 values.
# Output: results/protected_check/<case>.rds (same format as 33g)
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(TRUE)
NCORES <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
suppressPackageStartupMessages(library(parallel))
source("calibration/00_functions.R")              # indices()
B1 <- 4000L; B2 <- 10000L; NBIS <- 40L
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))              # grid of delta_hat of 27c (definition of the statistic)
D1E3 <- any(grepl("^--domain=1e3$", commandArgs(TRUE)))      # decision domain of the paths of 27c: [0, 1e3], 41 values
if (D1E3) GRID2 <- c(0, 10^seq(-4, 3, length.out = 41))
D1 <- c(0, 10^seq(-3, 10, by = 0.25))
OUT <- "calibration/results/protected_check"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CASES <- data.frame(
  id = c("PEX6_YRI_lrpath_step2", "FLVCR1_YRI_lrpath_step2"),
  symbol = c("PEX6-region", "FLVCR1-AS1"), gene = c("ENSG00000124587.9", "ENSG00000198468.2"),
  pos_before = c(43186801L, 213015222L), pos = c(43011766L, 213326573L), stringsAsFactors = FALSE)
CASES$tag <- sub("\\..*$", "", CASES$gene)
if (D1E3) CASES$id <- paste0(CASES$id, "_D1e3")
V27 <- read.csv("calibration/results/geuvadis_app_cov_variants.csv", stringsAsFactors = FALSE)

## ---- data (as in 27c) ------------------------------------------------------------------------
D <- "data/geuvadis"
EXPR <- { f <- file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")
  hdr <- strsplit(readLines(gzfile(f), n = 1L), "\t")[[1]]
  rows <- list(); con <- gzfile(f); open(con); readLines(con, n = 1L)
  while (length(l <- readLines(con, n = 4000L))) { for (g in CASES$gene) { h <- grep(paste0("^", g, "\t"), l)
    if (length(h)) rows[[g]] <- strsplit(l[h[1]], "\t")[[1]] }
    if (length(rows) == nrow(CASES)) break }
  close(con)
  lapply(rows, function(r) setNames(as.numeric(r[-(1:4)]), hdr[-(1:4)])) }
POPS <- list(YRI = list(dir = file.path(D, "loci_yri"), samples = readLines(file.path(D, "geuvadis_yri_overlap.txt"))))
COV <- read.csv(file.path(D, "geuvadis_covariates.csv"), stringsAsFactors = FALSE)
covariates <- function(samples) { d <- COV[match(samples, COV$sample), ]; stopifnot(!anyNA(d$sex))
  cbind(male = as.numeric(d$sex == "male")) }
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_window <- function(tag) {
  f <- file.path(POPS$YRI$dir, paste0(tag, ".gt"))
  L <- strsplit(readLines(f), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  sf <- file.path(POPS$YRI$dir, paste0(tag, ".samples"))
  samples <- if (file.exists(sf)) readLines(sf) else POPS$YRI$samples
  stopifnot(nrow(G) == length(samples))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf); keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  list(G = G[, keep, drop = FALSE], pos = pos[keep], samples = samples)
}

## ---- reduced model on the candidate-specific kernel (copied from 27c) ---------------------------
spectral_step <- function(X, y, C, sel, j) {
  n <- nrow(X); m <- ncol(X)
  F <- cbind(1, C, X[, sel, drop = FALSE]); p <- ncol(F)
  H <- qr.Q(qr(F), complete = TRUE)[, -(1:p), drop = FALSE]
  Rk <- setdiff(seq_len(m), c(sel, j)); Kjk <- tcrossprod(X[, Rk, drop = FALSE]) / length(Rk)
  M <- crossprod(H, Kjk %*% H); e <- eigen((M + t(M)) / 2, symmetric = TRUE); keep <- e$values > 1e-8 * max(e$values)
  U <- e$vectors[, keep, drop = FALSE]; ell <- e$values[keep]
  xt <- as.numeric(crossprod(H, X[, j])); yt <- as.numeric(crossprod(H, y))
  xi <- as.numeric(crossprod(U, xt)); xperp <- xt - as.numeric(U %*% xi); xp <- sqrt(sum(xperp^2))
  sk <- as.numeric(crossprod(U, yt))
  nperp <- as.integer(xp > 1e-6 * sqrt(sum(xt^2)))
  sp_ <- if (nperp) sum(xperp * yt) / xp else 0
  ssr <- max(sum(yt^2) - sum(sk^2) - sp_^2, 0)
  ek <- eigen(Kjk, symmetric = TRUE); UF <- crossprod(ek$vectors, F); lk <- pmax(ek$values, 0)
  ldF <- vapply(GRID2, function(dl) as.numeric(determinant(crossprod(UF, UF / (1 + dl * lk)))$modulus), 0)
  nu <- ncol(H)
  list(ell = ell, xi = xi, xp = if (nperp) xp else 0, nu0 = nu, r = length(ell), nperp = nperp,
       nrest = nu - length(ell) - nperp, n = n, ldF = ldF - ldF[1], obs = list(sk = sk, sp = sp_, ssr = ssr))
}
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
stats_from_coords <- function(sp, Sk, Sp, Ssr) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- sp$n
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
  RSS2 <- A2 - sweep(B2^2, 2, S2g, "/")
  l1 <- -0.5 * sweep((nu - 1) * log(RSS2), 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, GRID2)
  l0 <- -0.5 * sweep(nu * log(A2), 2, ld2, "+"); dh0 <- pick_grid(l0, GRID2)
  ldm <- ld2 - sp$ldF
  at <- function(dl) { iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Sk^2 * iD) + Sp^2 + Ssr; Bx <- colSums(xi * Sk * iD) + xp * Sp; S <- colSums(xi^2 * iD) + xp^2
    list(A = A, Bx = Bx, S = S, RSS = A - Bx^2 / S) }
  Q1 <- at(dh1); Q0 <- at(dh0)
  list(t = (Q1$Bx / Q1$S) * sqrt((nu - 1) * Q1$S / Q1$RSS), Q = Q0$Bx / sqrt(Q0$S * Q0$A / nu),
       lrt = apply(-sweep(n * log(RSS2), 2, ldm, "+"), 1, max) - apply(-sweep(n * log(A2), 2, ldm, "+"), 1, max),
       dh1 = dh1, dh0 = dh0, l1 = l1)
}
mc_summary <- function(lp) {
  m <- max(lp); w <- exp(lp - m); p <- mean(w); se <- sd(w) / sqrt(length(w))
  c(log10p = (m + log(p)) / log(10), log10up = (m + log(p + 2 * se)) / log(10), rel_se = se / p)
}
exact_lrt <- function(sp, delta0, obs_lrt, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); stats_from_coords(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  roots <- function(f, level, lo0, hi0) { lo <- rep(lo0, B); hi <- rep(hi0, B)
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- f(mid) < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    (lo + hi) / 2 }
  split <- roots(function(u) ev(u)$Q, 0, -60, 60)                    # as in 27c
  fL <- function(u) ev(u)$lrt
  ur <- { lo <- split; hi <- rep(200, B)
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- fL(mid) < obs_lrt; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    (lo + hi) / 2 }
  ul <- { lo <- rep(-200, B); hi <- split
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; dn <- fL(mid) < obs_lrt; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
    (lo + hi) / 2 }
  a1 <- pnorm(-ur, log.p = TRUE); a2 <- pnorm(ul, log.p = TRUE)
  mc_summary(pmax(a1, a2) + log1p(exp(-abs(a1 - a2))))
}

run_case <- function(k) {
  cs <- CASES[k, ]; ck <- file.path(OUT, paste0(cs$id, ".rds"))
  if (file.exists(ck)) return(paste("cached", cs$id))
  set.seed(3390 + k); t0 <- Sys.time()
  W <- load_window(cs$tag); X <- scale(W$G); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[cs$gene]][W$samples])); C <- covariates(W$samples)
  jb <- match(cs$pos_before, W$pos); j <- match(cs$pos, W$pos); stopifnot(!is.na(jb), !is.na(j))
  sp <- spectral_step(X, y, C, jb, j)
  o <- stats_from_coords(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
  n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m)
  v <- V27[V27$pop == "YRI" & V27$gene == cs$symbol & V27$pos == cs$pos & V27$step == 2, ][1, ]
  nominal <- pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  R27 <- readRDS(file.path("calibration/results/geuvadis_app_cov", sprintf("YRI_%s.rds", cs$tag)))
  rc <- R27$rec$lrt; lrt_decision <- rc$lrt[rc$step == 2 & rc$j == j][1]
  prof1 <- t(vapply(D1, function(d0) exact_lrt(sp, d0, o$lrt, B1), numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, function(d0) exact_lrt(sp, d0, o$lrt, B2), numeric(3)))
  prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2)); prof$id <- cs$id
  i <- which.max(prof$log10p)
  info <- data.frame(id = cs$id, pop = "YRI", gene = cs$symbol, pos = cs$pos, stat = "lrt (path of this study)", step = 2L,
                     n = n, m = m, nu = sp$nu0, rank = sp$r, nrest = sp$nrest, ell_min = min(sp$ell), ell_max = max(sp$ell),
                     delta_hat_obs = o$dh1, level = pchisq(pen, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                     nominal = nominal, stored_protected = v$log10p_sup_lrt, stored_argmax = v$argsup_lrt,
                     t_check = nominal, t_stored = v$log10p_nom_lrt, score_check = o$lrt, score_gemma = lrt_decision,
                     new_max = prof$log10p[i], new_argmax = prof$delta0[i], controlled = max(prof$log10up),
                     at_1e4 = prof1[D1 == 1e4, "log10p"], at_1e10 = prof1[D1 == 1e10, "log10p"],
                     elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(info = info, profile = prof), ck)
  sprintf("done %s (%.1f min)", cs$id, info$elapsed_min)
}
st <- mclapply(seq_len(nrow(CASES)), function(k) tryCatch(run_case(k), error = function(e) paste("FAIL", CASES$id[k], conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st))
