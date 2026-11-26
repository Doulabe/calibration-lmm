# ==============================================================================
# 33g_protected_check.R
# Numerical control of the protected exact-law p-values of the GEUVADIS application
# (33e_geuvadis_comparators_cov.R) for the variant and test cases on which a conclusion of
# the article rests. 33e approximated the supremum over delta0 >= 0 of the exact-law
# p-value p(delta0) by its maximum over 16 values of delta0 in [0, 10^4], each estimated
# with 1,000 draws of the error contrasts, the search of delta_hat being bounded at 10^4.
# Here
#   - the search grid of delta_hat follows GEMMA's range, lambda in [1e-5, 1e5] (0.1
#     decade; 33e used [0, 1e4]); it cannot be extended further, because the
#     maximum-likelihood profile of the null model, used by GEMMA's LRT and score, grows
#     like (1/2) log(delta) without bound (33h_ml_profile_check.R), so that GEMMA's
#     statistics are defined on this bounded range;
#   - stage 1 evaluates p(delta0) at 0 and 10^-3, ..., 10^10 by 0.25 decade (B1 draws);
#   - stage 2 refines by 0.05 decade within 0.3 decade of the stage-1 maximum (B2 draws);
#   - every value comes with its Monte Carlo standard error (standard deviation of the
#     conditional probabilities divided by sqrt(B)), and the controlled protected value is
#     the largest upper bound p_hat + 2 se over both stages.
# Cases: the four cases of Table 5 of the article (PEX6-region EUR Wald, ERAP2 EUR
# Kenward-Roger, TRA2A-AS YRI score, ZNF266 YRI Wald), the score at ZNF266 YRI and the
# Wald test at the PEX6-region lead of the Yoruba sample, selected by GEMMA's Wald path.
# Functions copied from 33e (spectral_fixed, pick_grid, gemma_stats, kr_factors), the
# Monte Carlo part of exact_gemma rewritten to return standard errors.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/33g_protected_check.R [cores=3] [--time]
# Output: results/protected_check/<case>.rds, results/protected_check_{profile,summary}.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(TRUE)
NCORES <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
TIME_ONLY <- "--time" %in% args
suppressPackageStartupMessages(library(parallel))
source("calibration/00_functions.R")              # indices()
B1 <- 4000L; B2 <- 10000L; NBIS <- 40L
GRID2 <- c(0, 10^seq(-5, 5, by = 0.1))                     # search grid of delta_hat, GEMMA's range (33e: up to 10^4)
D1 <- c(0, 10^seq(-3, 10, by = 0.25))                      # stage-1 values of delta0
OUT <- "calibration/results/protected_check"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
COMP <- "calibration/results/geuvadis_comp_cov"

CASES <- data.frame(
  id   = c("PEX6_EUR_wald", "ERAP2_EUR_kr", "TRA2A_YRI_score", "ZNF266_YRI_wald", "ZNF266_YRI_score", "PEX6_YRI_wald"),
  pop  = c("EUR", "EUR", "YRI", "YRI", "YRI", "YRI"),
  symbol = c("PEX6-region", "ERAP2", "TRA2A-AS", "ZNF266", "ZNF266", "PEX6-region"),
  pos  = c(42944850L, 96252432L, 23134790L, 9477868L, 9477868L, 42915920L),
  stat = c("wald", "kr", "score", "wald", "score", "wald"), stringsAsFactors = FALSE)

## ---- data (as in 33e) ------------------------------------------------------------------------
GENES <- data.frame(                                       # all eight genes (33i sources this script)
  gene = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4", "ENSG00000203875.4",
           "ENSG00000198468.2", "ENSG00000124587.9", "ENSG00000230658.1", "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  stringsAsFactors = FALSE)
GENES$tag <- sub("\\..*$", "", GENES$gene)
D <- "data/geuvadis"
EXPR <- { f <- file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")
  hdr <- strsplit(readLines(gzfile(f), n = 1L), "\t")[[1]]
  rows <- list(); con <- gzfile(f); open(con); readLines(con, n = 1L)
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

## ---- reduced model with GEMMA's fixed kernel (copied from 33e) ----------------------------------
spectral_fixed <- function(X, y, C, sel, j, Kfull) {
  n <- nrow(X); F <- cbind(1, C, X[, sel, drop = FALSE]); p <- ncol(F)
  H <- qr.Q(qr(F), complete = TRUE)[, -(1:p), drop = FALSE]
  M <- crossprod(H, Kfull %*% H); e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values); U <- e$vectors[, keep, drop = FALSE]; ell <- e$values[keep]
  xt <- as.numeric(crossprod(H, X[, j])); yt <- as.numeric(crossprod(H, y))
  xi <- as.numeric(crossprod(U, xt)); xperp <- xt - as.numeric(U %*% xi); xp <- sqrt(sum(xperp^2))
  sk <- as.numeric(crossprod(U, yt)); nperp <- as.integer(xp > 1e-6 * sqrt(sum(xt^2)))
  sp_ <- if (nperp) sum(xperp * yt) / xp else 0
  ssr <- max(sum(yt^2) - sum(sk^2) - sp_^2, 0); nu <- ncol(H)
  ek <- eigen(Kfull, symmetric = TRUE); UF <- crossprod(ek$vectors, F); lk <- pmax(ek$values, 0)
  ldF <- vapply(GRID2, function(dl) as.numeric(determinant(crossprod(UF, UF / (1 + dl * lk)))$modulus), 0)
  list(ell = ell, xi = xi, xp = if (nperp) xp else 0, nu0 = nu, r = length(ell), nperp = nperp,
       nrest = nu - length(ell) - nperp, n = n, ldF = ldF - ldF[1], obs = list(sk = sk, sp = sp_, ssr = ssr),
       ell_min = min(ell), ell_max = max(ell))
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
gemma_stats <- function(sp, Sk, Sp, Ssr) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- sp$n
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
  RSS2 <- A2 - sweep(B2^2, 2, S2g, "/")
  l1 <- -0.5 * sweep((nu - 1) * log(RSS2), 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, GRID2)   # REML, H1
  ldm <- ld2 - sp$ldF
  ml0 <- -0.5 * sweep(n * log(A2), 2, ldm, "+"); ml1 <- -0.5 * sweep(n * log(RSS2), 2, ldm, "+")
  dm0 <- pick_grid(ml0, GRID2)                                                                   # ML, H0
  at <- function(dl) { iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Sk^2 * iD) + Sp^2 + Ssr; Bx <- colSums(xi * Sk * iD) + xp * Sp; S <- colSums(xi^2 * iD) + xp^2
    list(A = A, Bx = Bx, S = S, RSS = A - Bx^2 / S) }
  Q1 <- at(dh1); Q0 <- at(dm0)
  list(t = (Q1$Bx / Q1$S) * sqrt((nu - 1) * Q1$S / Q1$RSS),
       Fs = n * Q0$Bx^2 / (Q0$S * Q0$A),
       lrt = 2 * (apply(ml1, 1, max) - apply(ml0, 1, max)), dh1 = dh1, dm0 = dm0, l1 = l1)
}
kr_factors <- function(ell, xi2, xp2, nout, dl) {
  w <- 1 / (1 + outer(ell, dl)); w2 <- w * w; w3 <- w2 * w
  P <- colSums(xi2 * w) + xp2; Pg1 <- colSums(xi2 * ell * w2); Pe1 <- colSums(xi2 * w2) + xp2
  Ugg <- colSums(xi2 * ell^2 * w3); Uge <- colSums(xi2 * ell * w3); Uee <- colSums(xi2 * w3) + xp2
  Tgg <- colSums(ell^2 * w2); Tge <- colSums(ell * w2); Tee <- colSums(w2) + nout
  Igg <- 0.5 * (Tgg - 2 * Ugg / P + Pg1^2 / P^2); Ige <- 0.5 * (Tge - 2 * Uge / P + Pg1 * Pe1 / P^2)
  Iee <- 0.5 * (Tee - 2 * Uee / P + Pe1^2 / P^2); det <- Igg * Iee - Ige^2
  hgg <- Iee / det; hge <- -Ige / det; hee <- Igg / det
  A <- (hgg * Pg1^2 + 2 * hge * Pg1 * Pe1 + hee * Pe1^2) / P^2
  LoP <- (hgg * (Ugg - Pg1^2 / P) + 2 * hge * (Uge - Pg1 * Pe1 / P) + hee * (Uee - Pe1^2 / P)) / P
  list(ratio = 1 + 2 * LoP, df = 2 / A)
}

## ---- exact-law p-value of one statistic at delta0, with its Monte Carlo standard error ----------
mc_summary <- function(lp) {                  # lp: log conditional probabilities, one per draw of r
  m <- max(lp); w <- exp(lp - m); p <- mean(w); se <- sd(w) / sqrt(length(w))
  c(log10p = (m + log(p)) / log(10), log10up = (m + log(p + 2 * se)) / log(10), rel_se = se / p)
}
exact_one <- function(sp, delta0, obs, stat, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); gemma_stats(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  z0 <- ev(rep(0, B))
  lsum <- function(a1, a2) pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))
  if (stat %in% c("wald", "kr")) {
    zz <- ev(rep(1, B)); wa <- zz$t - z0$t; wg <- z0$t
    c0 <- if (stat == "wald") abs(obs$t) else {
      kf <- kr_factors(ell, xi^2, xp^2, sp$nperp + sp$nrest, z0$dh1)
      qt(obs$p_kr / 2, df = kf$df, lower.tail = FALSE) * sqrt(kf$ratio) }
    return(mc_summary(lsum(pnorm(-(c0 - wg) / wa, log.p = TRUE), pnorm((-c0 - wg) / wa, log.p = TRUE))))
  }
  f <- if (stat == "score") function(u) ev(u)$Fs else function(u) ev(u)$lrt
  level <- if (stat == "score") obs$Fs else obs$lrt
  split <- { lo <- rep(-60, B); hi <- rep(60, B)
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- ev(mid)$t < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
  ur <- { lo <- split; hi <- rep(250, B)
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- f(mid) < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
  ul <- { lo <- rep(-250, B); hi <- split
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; dn <- f(mid) < level; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }; (lo + hi) / 2 }
  mc_summary(lsum(pnorm(-ur, log.p = TRUE), pnorm(ul, log.p = TRUE)))
}

## ---- one case -------------------------------------------------------------------------------------
setup_case <- function(cs) {
  g <- GENES[GENES$symbol == cs$symbol, ]
  R <- readRDS(file.path(COMP, sprintf("%s_%s.rds", cs$pop, g$tag)))
  v <- R$variants[R$variants$pos == cs$pos, ][1, ]
  before <- NULL
  for (tn in strsplit(v$tests, "+", fixed = TRUE)[[1]]) { idx <- R$gemma_paths[[tn]]$idx
    if (length(idx) >= v$step && idx[v$step] == v$j) { before <- idx[seq_len(v$step - 1L)]; break } }
  stopifnot(!is.null(before))
  W <- load_window(cs$pop, g$tag)
  X <- scale(W$G); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); C <- covariates(cs$pop, W$samples)
  n <- nrow(X); m <- ncol(X); K <- tcrossprod(X) / m
  stopifnot(W$pos[v$j] == cs$pos)
  sp <- spectral_fixed(X, y, C, before, v$j, K)
  o <- gemma_stats(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
  pen <- log(n) + 2 * log(m); ncf <- 1L + ncol(C); nu_lev <- n - ncf - v$step + 1L
  level <- log10(2 * pt(-sqrt((nu_lev - 1) * expm1(pen / n)), nu_lev - 1))        # eBIC per-test level (33e)
  stored <- switch(cs$stat, wald = v$log10p_sup_wald, kr = v$log10p_sup_kr, score = v$log10p_sup_score)
  argst <- switch(cs$stat, wald = v$argsup_wald, kr = v$argsup_kr, score = v$argsup_score)
  nominal <- switch(cs$stat, wald = v$log10p_gemma_wald, kr = v$log10p_kr_full, score = v$log10p_gemma_score)
  list(sp = sp, obs = list(t = o$t, Fs = o$Fs, lrt = o$lrt, p_kr = 10^v$log10p_kr_full),
       info = data.frame(id = cs$id, pop = cs$pop, gene = cs$symbol, pos = cs$pos, stat = cs$stat, step = v$step,
                         n = n, m = m, nu = sp$nu0, rank = sp$r, nrest = sp$nrest, ell_min = sp$ell_min,
                         ell_max = sp$ell_max, delta_hat_obs = o$dh1, level = level, nominal = nominal,
                         stored_protected = stored, stored_argmax = argst,
                         t_check = (log(2) + pt(-abs(o$t), sp$nu0 - 1, log.p = TRUE)) / log(10),
                         t_stored = v$log10p_ours_wald,
                         score_check = pf(o$Fs, 1, sp$nu0 - 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                         score_gemma = v$log10p_gemma_score))
}
run_case <- function(k) {
  cs <- CASES[k, ]; ck <- file.path(OUT, paste0(cs$id, ".rds"))
  if (file.exists(ck)) return(paste("cached", cs$id))
  set.seed(3370 + k); S <- setup_case(cs); t0 <- Sys.time()
  prof1 <- t(vapply(D1, function(d0) exact_one(S$sp, d0, S$obs, cs$stat, B1), numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, function(d0) exact_one(S$sp, d0, S$obs, cs$stat, B2), numeric(3)))
  prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2))
  prof$id <- cs$id
  i <- which.max(prof$log10p)
  res <- list(info = cbind(S$info, new_max = prof$log10p[i], new_argmax = prof$delta0[i],
                           controlled = max(prof$log10up), at_1e4 = prof1[D1 == 1e4, "log10p"],
                           at_1e10 = prof1[D1 == 1e10, "log10p"],
                           elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))),
              profile = prof)
  saveRDS(res, ck)
  sprintf("done %s (%.1f min)", cs$id, res$info$elapsed_min)
}

if (identical(Sys.getenv("PC_SOURCE_ONLY"), "1")) stop("sourced for its functions only", call. = FALSE)
if (TIME_ONLY) {                                   # timing of one evaluation per case, then stop
  for (k in seq_len(nrow(CASES))) { cs <- CASES[k, ]; S <- setup_case(cs)
    tm <- system.time(e <- exact_one(S$sp, 1, S$obs, cs$stat, B1))[["elapsed"]]
    print(cbind(S$info[, c("id", "nu", "level", "nominal", "stored_protected", "t_check", "t_stored",
                           "score_check", "score_gemma")], sec_per_eval_B1 = tm, log10p_at_1 = e[["log10p"]])) }
  quit(save = "no")
}
st <- mclapply(seq_len(nrow(CASES)), function(k) tryCatch(run_case(k), error = function(e) paste("FAIL", CASES$id[k], conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st))
R <- lapply(list.files(OUT, "[.]rds$", full.names = TRUE), readRDS)
write.csv(do.call(rbind, lapply(R, `[[`, "profile")), "calibration/results/protected_check_profile.csv", row.names = FALSE)
S <- do.call(rbind, lapply(R, `[[`, "info"))
write.csv(S, "calibration/results/protected_check_summary.csv", row.names = FALSE)
print(S[, c("id", "level", "nominal", "stored_protected", "stored_argmax", "new_max", "new_argmax", "controlled", "at_1e4", "at_1e10")], digits = 3)
