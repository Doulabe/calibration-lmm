# ==============================================================================
# 27c_geuvadis_application_cov.R
#
# GEUVADIS application with sample covariates in the fixed-effect block:
# F = (intercept, covariates, variants already selected), the covariates being sex and,
# for EUR, three ancestry indicators (FIN, TSI, CEU; 33d_geuvadis_covariates.R).
# Same 8 cis windows, genotypes and expression (GD462 resk10, standardized) as
# 27_geuvadis_application.R, which had the intercept only.
#   1. stepwise paths of the article on the candidate-specific kernel K_jk (wald_h1: t at
#      delta_hat; score_h0: t at the null estimate delta_hat_0, as EMMAX and MLMM; lrt:
#      maximum-likelihood ratio with delta under each model), eBIC level log n + 2 log m;
#      conditional regression (no random effect, eBIC) and SuSiE (L = 10, top PIP per 95%
#      credible set) on expression and genotypes residualized on the covariates.
#   2. for every variant selected by the three paths, at its step: nu_j, nominal p-values
#      and p-values from the exact laws (Rao-Blackwellized, 1,000 draws of r) of the three
#      statistics, evaluated
#        - at delta0 = delta_hat (plug-in, not a valid p-value),
#        - maximized over the 95 per cent profile-REML interval of delta (descriptive),
#        - maximized over delta0 in {0, 10^-3, 10^-2.5, ..., 10^4} and the interval points:
#          the supremum over the parameter space, a valid (conservative) p-value.
#      The maximum-likelihood statistic includes the term log|F'V^-1F| (33c_fix_ml_term.R):
#      log|V| = log|H'VH| - log|F'V^-1F| + const, which depends on delta at every step once
#      F contains covariates.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/27c_geuvadis_application_cov.R [cores=3] [--pop=EUR|YRI|ALL] [--only=GENE,...]
# Output: results/geuvadis_app_cov/<pop>_<tag>.rds, results/geuvadis_app_cov_{paths,variants}.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
POP    <- { a <- grep("^--pop=", args, value = TRUE); if (length(a)) sub("^--pop=", "", a) else "ALL" }
suppressPackageStartupMessages({ library(parallel); library(susieR) })
source("calibration/lib/stepwise_baseline.R")    # stepwise_baseline() (conditional regression)
source("calibration/00_functions.R")              # indices()
B_EX  <- 1000L
K_MAX <- 10L
GRIDD <- c(0, 10^seq(-4, 3, length.out = 41))
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
DSUP  <- c(0, 10^seq(-3, 4, by = 0.5))                     # delta0 values of the supremum p-value

## ---- engine algebra and paths (14_selection_demo.R, without the marginal Bayes factors) ---------
pieces_x <- function(Uty, UtF, UtX, s_rank, delta, a, n, k) {
  d_inv <- 1 / (1 + delta * s_rank); WU_y <- d_inv * Uty
  q0 <- sum(Uty * WU_y); p0 <- as.numeric(crossprod(UtF, WU_y))
  M0 <- crossprod(UtF, d_inv * UtF); L0 <- chol(M0 + 1e-12 * diag(k))
  w0 <- backsolve(L0, forwardsolve(t(L0), p0)); S00 <- sum(p0 * w0)
  logdetM0 <- 2 * sum(log(diag(L0))); logdetA <- sum(log1p(delta * s_rank))
  alpha <- colSums(d_inv * UtX^2); beta <- as.numeric(crossprod(UtX, WU_y))
  FWX <- crossprod(UtF, d_inv * UtX); H <- backsolve(L0, forwardsolve(t(L0), FWX))
  phi <- colSums(FWX * H); psi <- as.numeric(crossprod(FWX, w0))
  t_j <- pmax(1 - a * alpha, 1e-12); g_j <- a / t_j; denom <- 1 + g_j * phi
  yZy <- q0 + g_j * beta^2; pMp <- S00 + 2 * g_j * beta * psi + g_j^2 * beta^2 * phi
  hp <- psi + g_j * beta * phi; S0 <- pMp - g_j * hp^2 / denom
  RSS0 <- pmax(yZy - S0, 1e-300)
  D <- pmax(alpha / t_j - phi / (denom * t_j^2), 1e-300); u <- (beta - hp / denom) / t_j
  list(D = D, u = u, RSS0 = RSS0, ldZ = logdetA + log(t_j), ldM = logdetM0 + log(denom))
}
step_stats <- function(y, Xc, F_k, U, s_raw, m_k) {
  n <- length(y); k <- ncol(F_k); nu <- n - k; mc <- ncol(Xc); G <- length(GRIDD)
  s_rank <- pmax(s_raw, 0) / m_k
  Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F_k); UtX <- crossprod(U, Xc)
  Dm <- um <- R0 <- lZ <- lM <- matrix(0, mc, G)
  for (g in seq_len(G)) {
    pc <- pieces_x(Uty, UtF, UtX, s_rank, GRIDD[g], GRIDD[g] / m_k, n, k)
    Dm[, g] <- pc$D; um[, g] <- pc$u; R0[, g] <- pc$RSS0; lZ[, g] <- pc$ldZ; lM[, g] <- pc$ldM
  }
  R1 <- pmax(R0 - um^2 / Dm, 1e-300)
  l0 <- -0.5 * (lZ + lM) - 0.5 * nu * log(R0)
  l1 <- -0.5 * (lZ + lM + log(Dm)) - 0.5 * (nu - 1) * log(R1)
  g1 <- max.col(l1, ties.method = "first"); g0 <- max.col(l0, ties.method = "first")
  ix1 <- cbind(seq_len(mc), g1); ix0 <- cbind(seq_len(mc), g0)
  t2w <- (nu - 1) * (um[ix1]^2 / Dm[ix1]) / R1[ix1]
  Q2 <- nu * um[ix0]^2 / (Dm[ix0] * R0[ix0]); t2s <- (nu - 1) * Q2 / pmax(nu - Q2, 1e-12)
  lrt <- apply(-lZ - n * log(R1), 1, max) - apply(-lZ - n * log(R0), 1, max)   # lZ = log|V|: full ML
  list(t2w = t2w, t2s = t2s, lrt = lrt, nu = nu)
}
TREAT <- c("wald_h1", "score_h0", "lrt")
run_treatments <- function(y, X, C) {
  n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m); a_n <- cbind(1, C)
  cache <- new.env()
  get_step <- function(sel) {
    key <- paste0("k:", paste(sort(sel), collapse = ","))
    if (!is.null(cache[[key]])) return(cache[[key]])
    rem <- setdiff(seq_len(m), sel)
    raw <- tcrossprod(X[, rem, drop = FALSE]); eg <- eigen(raw, symmetric = TRUE)
    F_k <- cbind(a_n, X[, sel, drop = FALSE])
    st <- step_stats(y, X[, rem, drop = FALSE], F_k, eg$vectors, pmax(eg$values, 0), length(rem) - 1L)
    st$rem <- rem; cache[[key]] <- st; st
  }
  out <- list()
  for (tr in TREAT) {
    sel <- integer(0); rec <- list()
    for (k in seq_len(K_MAX)) {
      st <- get_step(sel); ct <- (st$nu - 1) * expm1(pen / n)
      rank <- switch(tr, wald_h1 = st$t2w, score_h0 = st$t2s, lrt = st$lrt)
      b <- which.max(rank)
      ok <- switch(tr, wald_h1 = st$t2w[b] > ct, score_h0 = st$t2s[b] > ct, lrt = st$lrt[b] > pen)
      rec[[k]] <- data.frame(step = k, j = st$rem[b], t2w = st$t2w[b], t2s = st$t2s[b], lrt = st$lrt[b],
                             level_t2 = ct, pen = pen, accepted = isTRUE(ok))
      if (!isTRUE(ok)) break
      sel <- c(sel, st$rem[b])
    }
    out[[tr]] <- list(sel = sel, rec = do.call(rbind, rec))
  }
  out$pen <- pen
  out
}

## ---- spectral setup of a candidate at a step (candidate-specific kernel), with the ML term -----
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
## statistics from spectral coordinates (columns = data sets): t at the REML delta_hat, score Q at the
## null REML delta_hat_0, ML likelihood ratio with log|V| = sum log(1 + delta l_i) - log|F'V^-1F| + const
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
## p-values (log10) of the observed statistics from their exact laws at a given delta0 (27_geuvadis_application.R)
exact_logp <- function(sp, delta0, obs_t, obs_Q, obs_lrt, B = B_EX, NBIS = 40L) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- sp$n
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); stats_from_coords(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  z0 <- ev(rep(0, B))
  zz <- ev(rep(1, B)); wa <- zz$t - z0$t; wg <- z0$t                 # t(delta_hat) is affine in u given r
  lg_lin <- function(c0) { up <- (c0 - wg) / wa; dn <- (-c0 - wg) / wa
    lp <- pmax(pnorm(-up, log.p = TRUE), pnorm(dn, log.p = TRUE)) +
          log1p(exp(-abs(pnorm(-up, log.p = TRUE) - pnorm(dn, log.p = TRUE))))
    (max(lp) + log(mean(exp(lp - max(lp))))) / log(10) }
  roots <- function(f, level, lo0, hi0) { lo <- rep(lo0, B); hi <- rep(hi0, B)
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- f(mid) < level
      lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    (lo + hi) / 2 }
  fQ <- function(u) ev(u)$Q                                          # |Q| and |t(delta_hat_0)| are equivalent
  qu <- roots(fQ, abs(obs_Q), -200, 200); ql <- roots(fQ, -abs(obs_Q), -200, 200)
  lpQ <- { a1 <- pnorm(-qu, log.p = TRUE); a2 <- pnorm(ql, log.p = TRUE)
    lp <- pmax(a1, a2) + log1p(exp(-abs(a1 - a2))); (max(lp) + log(mean(exp(lp - max(lp))))) / log(10) }
  split <- roots(fQ, 0, -60, 60)
  fL <- function(u) ev(u)$lrt
  ur <- { lo <- split; hi <- rep(200, B)
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- fL(mid) < obs_lrt; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    (lo + hi) / 2 }
  ul <- { lo <- rep(-200, B); hi <- split
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; dn <- fL(mid) < obs_lrt; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
    (lo + hi) / 2 }
  lpL <- { a1 <- pnorm(-ur, log.p = TRUE); a2 <- pnorm(ul, log.p = TRUE)
    lp <- pmax(a1, a2) + log1p(exp(-abs(a1 - a2))); (max(lp) + log(mean(exp(lp - max(lp))))) / log(10) }
  c(wald = lg_lin(abs(obs_t)), null = lpQ, lrt = lpL)
}

variant_report <- function(X, y, C, sel_before, j, pos) {
  sp <- spectral_step(X, y, C, sel_before, j)
  o <- stats_from_coords(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
  nu <- sp$nu0; dh <- o$dh1; ix <- indices(sp, max(dh, 1e-6))
  prof <- 2 * (max(o$l1) - o$l1[1, ]); ci <- range(GRID2[prof <= qchisq(0.95, 1)])
  dci <- unique(c(ci[1], dh, ci[2], exp(mean(log(pmax(ci, 1e-4))))))
  dall <- unique(c(dci, DSUP))
  lp <- t(vapply(dall, function(d0) exact_logp(sp, d0, o$t, o$Q, o$lrt), numeric(3)))
  inci <- dall %in% dci
  t0 <- o$Q * sqrt((nu - 1) / pmax(nu - o$Q^2, 1e-12))               # t at the null estimate
  data.frame(pos = pos[j], j = j, step = length(sel_before) + 1L, nu = nu, delta_hat = dh, delta_hat0 = o$dh0,
             delta_ci_lo = ci[1], delta_ci_hi = ci[2], nu_j = ix[["nu_delta"]], nubar_j = ix[["nu_total"]],
             kappa_j = ix[["kappa"]], t = o$t, Q = o$Q, t_null = t0, lrt = o$lrt,
             log10p_nom_wald = (log(2) + pt(-abs(o$t), nu - 1, log.p = TRUE)) / log(10),
             log10p_nom_null = (log(2) + pt(-abs(t0), nu - 1, log.p = TRUE)) / log(10),
             log10p_nom_lrt = pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
             log10p_ex_wald = lp[dall == dh, "wald"][1], log10p_ex_null = lp[dall == dh, "null"][1],
             log10p_ex_lrt = lp[dall == dh, "lrt"][1],
             log10p_ci_wald = max(lp[inci, "wald"]), log10p_ci_null = max(lp[inci, "null"]), log10p_ci_lrt = max(lp[inci, "lrt"]),
             log10p_sup_wald = max(lp[, "wald"]), log10p_sup_null = max(lp[, "null"]), log10p_sup_lrt = max(lp[, "lrt"]),
             argsup_wald = dall[which.max(lp[, "wald"])], argsup_null = dall[which.max(lp[, "null"])],
             argsup_lrt = dall[which.max(lp[, "lrt"])])
}

## ---- data ------------------------------------------------------------------------------------
GENES <- data.frame(
  gene = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4", "ENSG00000203875.4",
           "ENSG00000198468.2", "ENSG00000124587.9", "ENSG00000230658.1", "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  chr = c("5", "12", "17", "6", "1", "6", "7", "19"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L, 213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)
GENES$tag <- sub("\\..*$", "", GENES$gene)
D <- "data/geuvadis"
EXPR <- { hdr <- strsplit(readLines(gzfile(file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")), n = 1L), "\t")[[1]]
  rows <- list(); con <- gzfile(file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")); open(con); readLines(con, n = 1L)
  while (length(l <- readLines(con, n = 4000L))) { for (g in GENES$gene) { h <- grep(paste0("^", g, "\t"), l)
    if (length(h)) rows[[g]] <- strsplit(l[h[1]], "\t")[[1]] }
    if (length(rows) == nrow(GENES)) break }
  close(con)
  lapply(rows, function(r) setNames(as.numeric(r[-(1:4)]), hdr[-(1:4)])) }
POPS <- list(EUR = list(dir = file.path(D, "loci"), samples = readLines(file.path(D, "geuvadis_eur_overlap.txt"))),
             YRI = list(dir = file.path(D, "loci_yri"), samples = readLines(file.path(D, "geuvadis_yri_overlap.txt"))))
COV <- read.csv(file.path(D, "geuvadis_covariates.csv"), stringsAsFactors = FALSE)
covariates <- function(pop, samples) {
  d <- COV[match(samples, COV$sample), ]; stopifnot(!anyNA(d$sex))
  C <- cbind(male = as.numeric(d$sex == "male"))
  if (pop == "EUR") C <- cbind(C, FIN = as.numeric(d$ancestry == "FIN"), TSI = as.numeric(d$ancestry == "TSI"),
                               CEU = as.numeric(d$ancestry == "CEU"))
  C
}
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_window <- function(pop, tag) {
  f <- file.path(POPS[[pop]]$dir, paste0(tag, ".gt")); if (!file.exists(f)) return(NULL)
  L <- strsplit(readLines(f), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  sf <- file.path(POPS[[pop]]$dir, paste0(tag, ".samples"))
  samples <- if (file.exists(sf)) readLines(sf) else POPS[[pop]]$samples
  stopifnot(nrow(G) == length(samples))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf); keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  list(G = G[, keep, drop = FALSE], pos = pos[keep], samples = samples)
}

OUT <- "calibration/results/geuvadis_app_cov"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_gene <- function(pop, i) {
  g <- GENES[i, ]; ck <- file.path(OUT, sprintf("%s_%s.rds", pop, g$tag))
  if (file.exists(ck)) return("cached")
  W <- load_window(pop, g$tag); if (is.null(W)) return(paste("no genotypes", pop, g$symbol))
  set.seed(2750 + i + if (pop == "YRI") 50L else 0L)
  X <- scale(W$G); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); stopifnot(!anyNA(y))
  C <- covariates(pop, W$samples)
  n <- nrow(X); m <- ncol(X); t0 <- Sys.time()
  tr <- run_treatments(y, X, C)
  ## conditional regression and SuSiE on expression and genotypes residualized on (1, C)
  qf <- qr(cbind(1, C)); yr <- as.numeric(scale(qr.resid(qf, y))); Xr <- scale(qr.resid(qf, X)); Xr[!is.finite(Xr)] <- 0
  ols <- tryCatch(stepwise_baseline(yr, Xr, kernel = "none", delta_mode = "once", stop = "ebic", K_max = K_MAX)$indices,
                  error = function(e) NA_integer_)
  sus <- tryCatch({ f <- susie(Xr, yr, L = 10, verbose = FALSE); cs <- f$sets$cs
    if (is.null(cs)) integer(0) else unique(vapply(cs, function(ix) ix[which.max(f$pip[ix])], 1L)) },
    error = function(e) NA_integer_)
  todo <- unique(do.call(rbind, lapply(TREAT, function(trn) {
    s <- tr[[trn]]$sel; if (!length(s)) return(NULL)
    data.frame(key = vapply(seq_along(s), function(k) paste(c(s[seq_len(k - 1)], "|", s[k]), collapse = ","), ""),
               k = seq_along(s), trn = trn) })))
  reps <- list()
  if (!is.null(todo)) for (key in unique(todo$key)) {
    parts <- strsplit(key, ",")[[1]]; bar <- which(parts == "|")
    before <- as.integer(parts[seq_len(bar - 1)]); j <- as.integer(parts[bar + 1])
    rr <- variant_report(X, y, C, before, j, W$pos)
    rr$treatments <- paste(sort(unique(todo$trn[todo$key == key])), collapse = "+"); reps[[key]] <- rr
  }
  lead_j <- match(g$lead_pos, W$pos)
  sels <- c(lapply(tr[TREAT], `[[`, "sel"), list(ols_ebic = ols, susie = sus))
  tagr2 <- function(js) if (!length(js) || all(is.na(js)) || is.na(lead_j)) NA else max(cor(X[, js, drop = FALSE], X[, lead_j])^2)
  res <- list(pop = pop, gene = g$symbol, tag = g$tag, n = n, m = m, n_cov = ncol(C), pen = tr$pen, pos = W$pos, lead_j = lead_j,
              sel = sels, rec = lapply(tr[TREAT], `[[`, "rec"), variants = if (length(reps)) do.call(rbind, reps) else NULL,
              lead_tag_r2 = vapply(sels, tagr2, 0), elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(res, ck); sprintf("done %s %s (%.0f s)", pop, g$symbol, res$elapsed)
}
pops <- if (POP == "ALL") names(POPS) else POP
jobs <- expand.grid(i = seq_len(nrow(GENES)), pop = pops, stringsAsFactors = FALSE)
only <- grep("^--only=", args, value = TRUE)
if (length(only)) jobs <- jobs[GENES$symbol[jobs$i] %in% strsplit(sub("^--only=", "", only), ",")[[1]], , drop = FALSE]
jobs <- jobs[order(jobs$pop == "EUR", decreasing = TRUE), ]           # EUR windows (largest) first
st <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(run_gene(jobs$pop[k], jobs$i[k]),
  error = function(e) paste("FAIL", jobs$pop[k], GENES$symbol[jobs$i[k]], conditionMessage(e))),
  mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st))
fs <- list.files(OUT, "rds$", full.names = TRUE)
if (length(fs)) {
  R <- lapply(fs, readRDS)
  paths <- do.call(rbind, lapply(R, function(r) do.call(rbind, lapply(names(r$sel), function(mn) {
    s <- r$sel[[mn]]; data.frame(pop = r$pop, gene = r$gene, n = r$n, m = r$m, method = mn,
      K = if (length(s) == 1 && is.na(s)) NA else length(s),
      positions = if (length(s) && !all(is.na(s))) paste(r$pos[s], collapse = ";") else "",
      lead_tag_r2 = r$lead_tag_r2[[mn]]) }))))
  write.csv(paths, "calibration/results/geuvadis_app_cov_paths.csv", row.names = FALSE)
  vars <- do.call(rbind, lapply(R, function(r) if (is.null(r$variants)) NULL else cbind(pop = r$pop, gene = r$gene, r$variants)))
  if (!is.null(vars)) write.csv(vars, "calibration/results/geuvadis_app_cov_variants.csv", row.names = FALSE)
  cat("collected", length(fs), "gene x population runs\n")
}
