# ==============================================================================
# 33_geuvadis_comparators.R
# Published software on the GEUVADIS application (8 cis windows, EUR n = 358 and
# YRI n = 87-89; same genotypes, expression and windows as 27_geuvadis_application.R).
#   A. GEMMA 0.98.5 -lmm 4 with the regional kernel K = ZZ'/m held fixed, forward
#      conditional paths of its Wald, LRT and score tests at the per-test level of the
#      eBIC rule (t^2 > (nu - 1)(exp(pen/n) - 1), pen = log n + 2 log m) and down to 1e-3.
#      For every variant entering a GEMMA path, at its step and with GEMMA's kernel:
#      nu_j, GEMMA's nominal p-values, the exact p-values of GEMMA's three statistics
#      (Wald with REML delta under H1; LRT with ML delta under both models, kernel
#      determinant and fixed-effect term included; score with the ML null delta and
#      GEMMA's n multiplier), the Kenward-Roger p-value (full n x n formulas at the REML fit
#      of MM4LMM 3.0.3; MM4LMM's own KR at the first step, where it runs; the reduced closed
#      form of 32_kr_exact.R) and its exact p-value.
#   B. MLMM 0.1.1 (extBIC and mBonf), QTLtools 1.3.1 (1,000 permutations, gene-level
#      5%), FarmCPU (rMVP 1.4.6, Bonferroni 0.05/m, LD clumping r2 > 0.5), BSLMM
#      (GEMMA -bslmm 1, 50,000 burn-in, 200,000 sampling; PIP > 0.5).
#   C. GridLMM (Runcie and Crawford 2019): first-step scan, REML Wald and Bayes factor
#      on an h2 grid of step 0.01, each marker downdated from its own kernel, prior
#      precision 1/tau2 = 25 on the standardized marker effect.
# SuSiE-inf and FINEMAP-inf on the same data: 33b_susieinf_geuvadis.py.
# Usage (from the repository root): nice -n 10 Rscript calibration/33_geuvadis_comparators.R [cores=2]
# Output: results/geuvadis_comp/<pop>_<tag>.rds, results/geuvadis_comp_{paths,variants,scan}.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(TRUE)
NCORES <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
suppressPackageStartupMessages({ library(parallel); library(bigmemory); library(rMVP); library(MM4LMM); library(GridLMM) })
.libPaths(c(Sys.getenv("COMPARATORS_RLIB", "comparators/Rlib"), .libPaths()))
suppressPackageStartupMessages(library(mlmm))
source("calibration/00_functions.R")              # indices()
GEMMA <- normalizePath(Sys.getenv("GEMMA_BIN", "comparators/GEMMA-0.98.5/bin/gemma"))
QTLT <- normalizePath(Sys.getenv("QTLTOOLS_BIN", "comparators/qtltools/bin/QTLtools"))
B_EX <- 1000L; K_MAX <- 10L; P_LIB <- 1e-3
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
OUT <- "calibration/results/geuvadis_comp"; TMP <- file.path(OUT, "tmp")
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

## ---- data (as in 27_geuvadis_application.R) -------------------------------------------------------
GENES <- data.frame(
  gene = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4", "ENSG00000203875.4",
           "ENSG00000198468.2", "ENSG00000124587.9", "ENSG00000230658.1", "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  chr = c(5L, 12L, 17L, 6L, 1L, 6L, 7L, 19L),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L, 213049214L, 42944850L, 23143113L, 9544276L),
  tss = c(96211643L, 56435637L, 33570086L, 86388451L, 213031430L, 42946958L, 23145322L, 9546234L),
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

## ---- reduced model with a FIXED kernel (GEMMA's): candidate and selected variants inside K ----------
spectral_fixed <- function(X, y, sel, j, Kfull) {
  n <- nrow(X); F <- cbind(1, X[, sel, drop = FALSE]); p <- ncol(F)
  H <- qr.Q(qr(F), complete = TRUE)[, -(1:p), drop = FALSE]
  M <- crossprod(H, Kfull %*% H); e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values); U <- e$vectors[, keep, drop = FALSE]; ell <- e$values[keep]
  xt <- as.numeric(crossprod(H, X[, j])); yt <- as.numeric(crossprod(H, y))
  xi <- as.numeric(crossprod(U, xt)); xperp <- xt - as.numeric(U %*% xi); xp <- sqrt(sum(xperp^2))
  sk <- as.numeric(crossprod(U, yt)); nperp <- as.integer(xp > 1e-6 * sqrt(sum(xt^2)))
  sp_ <- if (nperp) sum(xperp * yt) / xp else 0
  ssr <- max(sum(yt^2) - sum(sk^2) - sp_^2, 0); nu <- ncol(H)
  ## ML needs log|V| = log|H'VH| + log|F'V^{-1}F| + const: the second term on GRID2, eigen of the full kernel
  ek <- eigen(Kfull, symmetric = TRUE); UF <- crossprod(ek$vectors, F); lk <- pmax(ek$values, 0)
  ldF <- vapply(GRID2, function(dl) as.numeric(determinant(crossprod(UF, UF / (1 + dl * lk)))$modulus), 0)
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
## GEMMA's statistics from spectral coordinates (columns = datasets)
gemma_stats <- function(sp, Sk, Sp, Ssr) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; n <- sp$n
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2))); S2g <- colSums(xi^2 * W2) + xp^2
  A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
  RSS2 <- A2 - sweep(B2^2, 2, S2g, "/")
  l1 <- -0.5 * sweep((nu - 1) * log(RSS2), 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, GRID2)   # REML, H1
  ml0 <- -0.5 * sweep(n * log(A2), 2, ld2 + sp$ldF, "+"); ml1 <- -0.5 * sweep(n * log(RSS2), 2, ld2 + sp$ldF, "+")
  dm0 <- pick_grid(ml0, GRID2)                                                                   # ML, H0
  at <- function(dl) { iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Sk^2 * iD) + Sp^2 + Ssr; Bx <- colSums(xi * Sk * iD) + xp * Sp; S <- colSums(xi^2 * iD) + xp^2
    list(A = A, Bx = Bx, S = S, RSS = A - Bx^2 / S) }
  Q1 <- at(dh1); Q0 <- at(dm0)
  list(t = (Q1$Bx / Q1$S) * sqrt((nu - 1) * Q1$S / Q1$RSS),
       Fs = n * Q0$Bx^2 / (Q0$S * Q0$A),                                                        # GEMMA score
       lrt = 2 * (apply(ml1, 1, max) - apply(ml0, 1, max)), dh1 = dh1, dm0 = dm0, l1 = l1)
}
kr_factors <- function(ell, xi2, xp2, nout, dl) {                  # as in 32_kr_exact.R
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
## Kenward-Roger on the full n x n model (Kenward and Roger 1997, as coded in MM4LMM), at a given REML fit
kr_full <- function(y, Xf, K, s2) {
  n <- length(y); G <- list(K, diag(n))
  Si <- solve(s2[1] * K + s2[2] * diag(n)); SiX <- Si %*% Xf
  Phi <- solve(crossprod(Xf, SiX)); beta <- Phi %*% crossprod(SiX, y)
  P <- lapply(G, function(Gk) -crossprod(SiX, Gk %*% SiX))
  Pi <- Si - SiX %*% Phi %*% t(SiX); PG <- lapply(G, function(Gk) Pi %*% Gk)
  I <- matrix(0, 2, 2); for (i in 1:2) for (j in 1:2) I[i, j] <- 0.5 * sum(t(PG[[i]]) * PG[[j]])
  W <- solve(I); Lam <- 0 * Phi
  for (i in 1:2) for (j in 1:2)
    Lam <- Lam + W[i, j] * (crossprod(SiX, G[[i]] %*% Si %*% G[[j]] %*% SiX) - P[[i]] %*% Phi %*% P[[j]])
  PhiA <- Phi + 2 * Phi %*% Lam %*% Phi; L <- c(rep(0, ncol(Xf) - 1), 1)
  v <- vapply(P, function(Pk) as.numeric(t(L) %*% Phi %*% Pk %*% Phi %*% L), 0) / as.numeric(t(L) %*% Phi %*% L)
  df <- 2 / as.numeric(t(v) %*% W %*% v); Fst <- as.numeric((t(L) %*% beta)^2 / (t(L) %*% PhiA %*% L))
  c(F = Fst, df = df, p = pf(Fst, 1, df, lower.tail = FALSE))
}
lmean_exp <- function(lp) (max(lp) + log(mean(exp(lp - max(lp))))) / log(10)
## exact p-values (log10) of GEMMA's statistics and of the KR test, nuisance delta0 given
exact_gemma <- function(sp, delta0, obs, B = B_EX, NBIS = 40L) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r; nu <- sp$nu0
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); gemma_stats(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  z0 <- ev(rep(0, B)); zz <- ev(rep(1, B)); wa <- zz$t - z0$t; wg <- z0$t
  lin <- function(c0) { a1 <- pnorm(-(c0 - wg) / wa, log.p = TRUE); a2 <- pnorm((-c0 - wg) / wa, log.p = TRUE)
    lmean_exp(pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))) }
  kf <- kr_factors(ell, xi^2, xp^2, sp$nperp + sp$nrest, z0$dh1)       # REML fit with the candidate: a function of r
  c_kr <- qt(obs$p_kr / 2, df = kf$df, lower.tail = FALSE) * sqrt(kf$ratio)
  two_sided <- function(f, level) {                                     # f increasing-then-decreasing around split
    split <- { lo <- rep(-60, B); hi <- rep(60, B)
      for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- ev(mid)$t < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
    ur <- { lo <- split; hi <- rep(250, B)
      for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- f(mid) < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
    ul <- { lo <- rep(-250, B); hi <- split
      for (it in 1:NBIS) { mid <- (lo + hi) / 2; dn <- f(mid) < level; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }; (lo + hi) / 2 }
    a1 <- pnorm(-ur, log.p = TRUE); a2 <- pnorm(ul, log.p = TRUE)
    lmean_exp(pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))) }
  c(wald = lin(abs(obs$t)), kr = lin(c_kr), lrt = two_sided(function(u) ev(u)$lrt, obs$lrt),
    score = two_sided(function(u) ev(u)$Fs, obs$Fs))
}

## ---- GEMMA ---------------------------------------------------------------------------------------
write_gemma_inputs <- function(w, G, K) {
  m <- ncol(G)
  data.table::fwrite(data.frame(paste0("s", seq_len(m)), "A", "B", t(G)), file.path(w, "geno.txt"), col.names = FALSE, sep = ",")
  data.table::fwrite(as.data.frame(K), file.path(w, "K.txt"), col.names = FALSE, sep = "\t")
}
gemma_lmm4 <- function(w, y, G_sel, tag) {
  write.table(y, file.path(w, "y.txt"), row.names = FALSE, col.names = FALSE)
  write.table(cbind(1, G_sel), file.path(w, "cov.txt"), row.names = FALSE, col.names = FALSE)
  st <- system2(GEMMA, c("-g", file.path(w, "geno.txt"), "-p", file.path(w, "y.txt"), "-k", file.path(w, "K.txt"),
                         "-c", file.path(w, "cov.txt"), "-lmm", "4", "-notsnp", "-outdir", w, "-o", tag),
                stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("GEMMA failed")
  A <- data.table::fread(file.path(w, paste0(tag, ".assoc.txt")), data.table = FALSE)
  data.frame(idx = as.integer(sub("^s", "", A$rs)), p_wald = A$p_wald, p_lrt = A$p_lrt, p_score = A$p_score,
             beta = A$beta, se = A$se, l_remle = A$l_remle, l_mle = A$l_mle)
}

run_gene <- function(pop, i) {
  g <- GENES[i, ]; ck <- file.path(OUT, sprintf("%s_%s.rds", pop, g$tag))
  if (file.exists(ck)) return("cached")
  W <- load_window(pop, g$tag); if (is.null(W)) return(paste("no genotypes", pop, g$symbol))
  set.seed(3300 + i + if (pop == "YRI") 50L else 0L)
  X <- scale(W$G); X[!is.finite(X)] <- 0; G <- W$G
  y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); stopifnot(!anyNA(y))
  n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m); K <- tcrossprod(X) / m
  w <- file.path(TMP, sprintf("%s_%s", pop, g$tag)); dir.create(w, showWarnings = FALSE)
  write_gemma_inputs(w, G, K)
  t0 <- Sys.time()
  ## A. GEMMA paths
  cache <- list()
  run <- function(sel) { key <- paste0("k:", paste(sort(sel), collapse = ","))
    if (is.null(cache[[key]])) cache[[key]] <<- gemma_lmm4(w, y, G[, sel, drop = FALSE], paste0("r", length(cache)))
    cache[[key]] }
  paths <- list()
  for (test in c("wald", "lrt", "score")) {
    sel <- integer(0); pv <- numeric(0)
    for (k in seq_len(K_MAX)) { A <- run(sel); p <- A[[paste0("p_", test)]]; b <- which.min(p)
      if (!length(b) || p[b] > P_LIB) break
      sel <- c(sel, A$idx[b]); pv <- c(pv, p[b]) }
    paths[[test]] <- data.frame(idx = sel, p = pv)
  }
  lev_wald <- function(nu) 2 * pt(-sqrt((nu - 1) * expm1(pen / n)), nu - 1)      # eBIC-equivalent per-test level
  stop_at <- function(pp) { if (!nrow(pp)) return(integer(0))
    lev <- vapply(seq_len(nrow(pp)), function(k) lev_wald(n - k), 0)
    ok <- pp$p <= lev; if (!ok[1]) integer(0) else pp$idx[seq_len(min(c(which(!ok) - 1, nrow(pp))))] }
  gsel <- lapply(paths, stop_at)
  ## variant reports (each path prefix once)
  todo <- unique(do.call(rbind, lapply(names(paths), function(tn) {   # selected at the eBIC level + the next entry
    s <- paths[[tn]]$idx; s <- s[seq_len(min(length(s), length(gsel[[tn]]) + 1L))]; if (!length(s)) return(NULL)
    data.frame(key = vapply(seq_along(s), function(k) paste(c(s[seq_len(k - 1)], "|", s[k]), collapse = ","), ""),
               test = tn) })))
  reps <- list()
  if (!is.null(todo)) for (key in unique(todo$key)) {
    parts <- strsplit(key, ",")[[1]]; bar <- which(parts == "|")
    before <- as.integer(parts[seq_len(bar - 1)]); j <- as.integer(parts[bar + 1])
    sp <- spectral_fixed(X, y, before, j, K)
    o <- gemma_stats(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
    gem <- run(before); gj <- gem[gem$idx == j, ]
    ## MM4LMM Kenward-Roger with GEMMA's kernel and covariates
    Fc <- cbind(1, X[, before, drop = FALSE])
    R <- MMEst(Y = y, Cofactor = Fc, X = matrix(X[, j], n, 1, dimnames = list(NULL, "x")),
               VarList = list(K = K, Error = diag(n)), CritVar = 1e-9, CritLogLik = 1e-9, MaxIter = 5000, Verbose = FALSE)
    krf <- kr_full(y, cbind(Fc, X[, j]), K, R[[1]]$Sigma2)            # full-model KR at MM4LMM's REML fit
    kr <- if (!length(before)) {                                          # MM4LMM's own KR (it fails with cofactors)
      AnovaTest(R, TestedCombination = matrix(c(0, 1), nrow = 1), Type = "KR", Cofactor = Fc,
                X = matrix(X[, j], n, 1, dimnames = list(NULL, "x")), VarList = list(K = K, Error = diag(n)))[[1]][, 1]
    } else c(Stat = NA, pval = NA, ddf = NA)
    dlmm <- R[[1]]$Sigma2[1] / R[[1]]$Sigma2[2]
    kfo <- kr_factors(sp$ell, sp$xi^2, sp$xp^2, sp$nperp + sp$nrest, o$dh1)
    t2 <- o$t^2; p_kr_ours <- pf(t2 / kfo$ratio, 1, kfo$df, lower.tail = FALSE)
    dh <- o$dh1; ix <- indices(sp, max(dh, 1e-6))
    prof <- 2 * (max(o$l1) - o$l1[1, ]); ci <- range(GRID2[prof <= qchisq(0.95, 1)])
    dgrid <- unique(c(ci[1], dh, ci[2]))
    obs <- list(t = o$t, lrt = o$lrt, Fs = o$Fs, p_kr = krf[["p"]])
    lp <- t(vapply(dgrid, function(d0) exact_gemma(sp, d0, obs), numeric(4)))
    reps[[key]] <- data.frame(pos = W$pos[j], j = j, step = length(before) + 1L, tests = paste(sort(unique(todo$test[todo$key == key])), collapse = "+"),
      nu = sp$nu0, delta_hat = dh, delta_ci_lo = ci[1], delta_ci_hi = ci[2], nu_j = ix[["nu_delta"]],
      gemma_l_remle = gj$l_remle, ours_dh1 = dh, mm4lmm_delta = dlmm,
      log10p_gemma_wald = log10(gj$p_wald), log10p_gemma_lrt = log10(gj$p_lrt), log10p_gemma_score = log10(gj$p_score),
      log10p_ours_wald = (log(2) + pt(-abs(o$t), sp$nu0 - 1, log.p = TRUE)) / log(10),
      log10p_ours_lrt = pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
      log10p_ours_score = pf(o$Fs, 1, sp$nu0 - 1, lower.tail = FALSE, log.p = TRUE) / log(10),
      log10p_kr_mm4lmm = log10(kr[["pval"]]), df_kr_mm4lmm = kr[["ddf"]], log10p_kr_full = log10(krf[["p"]]), df_kr_full = krf[["df"]],
      log10p_kr_ours = log10(p_kr_ours), df_kr_ours = kfo$df,
      log10p_ex_wald = lp[dgrid == dh, "wald"][1], log10p_ex_kr = lp[dgrid == dh, "kr"][1],
      log10p_ex_lrt = lp[dgrid == dh, "lrt"][1], log10p_ex_score = lp[dgrid == dh, "score"][1],
      log10p_exmax_wald = max(lp[, "wald"]), log10p_exmax_kr = max(lp[, "kr"]),
      log10p_exmax_lrt = max(lp[, "lrt"]), log10p_exmax_score = max(lp[, "score"]))
  }
  ## B. MLMM, QTLtools, FarmCPU, BSLMM
  smp <- sprintf("i%03d", seq_len(n))
  ml <- tryCatch({ Xn <- X; colnames(Xn) <- paste0("s", seq_len(m)); rownames(Xn) <- smp
    Kn <- K; dimnames(Kn) <- list(smp, smp); sink(tempfile())
    fit <- mlmm(setNames(y, smp), Xn, Kn, nbchunks = 2, maxsteps = 11); sink()
    cof <- function(v) { v <- v[!is.na(v)]; as.integer(sub("^s", "", v)) }
    list(extBIC = cof(fit$opt_extBIC$cof), mbonf = cof(fit$opt_mbonf$cof)) },
    error = function(e) { if (sink.number() > 0) sink(); conditionMessage(e) })
  qt <- tryCatch({
    bed <- data.frame(g$chr, W$pos - 1L, W$pos, paste0("s", seq_len(m)), paste0("s", seq_len(m)), "+", t(G), check.names = FALSE)
    colnames(bed) <- c("#chr", "start", "end", "id", "gid", "strand", smp); bed <- bed[order(W$pos), ]
    data.table::fwrite(bed, file.path(w, "genotypes.bed"), sep = "\t", quote = FALSE)
    ph <- data.frame(g$chr, g$tss - 1L, g$tss, g$symbol, g$symbol, "+", t(y), check.names = FALSE)
    colnames(ph) <- colnames(bed); data.table::fwrite(ph, file.path(w, "expression.bed"), sep = "\t", quote = FALSE)
    for (nm in c("genotypes", "expression")) { fn <- file.path(w, paste0(nm, ".bed"))
      stopifnot(system2("bgzip", c("-f", shQuote(fn))) == 0, system2("tabix", c("-f", "-p", "bed", shQuote(paste0(fn, ".gz")))) == 0) }
    win <- max(abs(W$pos - g$tss)) + 1L
    common <- c("cis", "--vcf", file.path(w, "genotypes.bed.gz"), "--bed", file.path(w, "expression.bed.gz"), "--window", win,
                "--seed", 25000 + i)
    per <- file.path(w, "perm.txt")
    stopifnot(system2(QTLT, c(common, "--permute", "1000", "--out", per), stdout = FALSE, stderr = FALSE) == 0)
    p <- read.table(per); pb <- qbeta(0.05, p$V14, p$V15); tt <- qf(pb, 1, p$V13, lower.tail = FALSE)
    corr2 <- tt / (p$V13 + tt); cutoff <- pf(p$V12 * corr2 / (1 - corr2), 1, p$V12, lower.tail = FALSE)
    thr <- file.path(w, "thr.txt"); write.table(data.frame(g$symbol, cutoff), thr, row.names = FALSE, col.names = FALSE, quote = FALSE)
    mp <- file.path(w, "map.txt")
    stopifnot(system2(QTLT, c(common, "--mapping", thr, "--out", mp), stdout = FALSE, stderr = FALSE) == 0)
    tab <- if (file.exists(mp) && file.info(mp)$size > 0) read.table(mp) else data.frame()
    hits <- if (nrow(tab)) tab[tab$V21 == 1, , drop = FALSE] else tab
    list(sel = if (nrow(hits)) as.integer(sub("^s", "", hits$V8)) else integer(0), cutoff = cutoff, perm_p = p$V20) },
    error = function(e) conditionMessage(e))
  fc <- tryCatch({ geno <- bigmemory::as.big.matrix(t(G))
    map <- data.frame(SNP = paste0("s", seq_len(m)), CHROM = g$chr, POS = W$pos)
    r <- MVP.FarmCPU(phe = data.frame(id = smp, y = y), geno = geno, map = map, CV = NULL, maxLoop = 10,
                     method.bin = "static", bin.size = c(1e4, 5e4, 1e5), bin.selection = seq(5, 50, 5),
                     p.threshold = 0.05 / m, ncpus = 1, verbose = FALSE)
    p <- r[, ncol(r)]; o <- order(p, na.last = NA); sel <- integer(0)     # NaN: exact duplicates of pseudo-QTNs
    for (j in o[p[o] <= 0.05 / m]) if (!length(sel) || max(cor(X[, j], X[, sel, drop = FALSE])^2) <= 0.5) sel <- c(sel, j)
    list(sel = sel, minp = min(p, na.rm = TRUE), nan = which(is.na(p))) }, error = function(e) conditionMessage(e))
  bs <- tryCatch({ write.table(y, file.path(w, "yb.txt"), row.names = FALSE, col.names = FALSE)
    st <- system2(GEMMA, c("-g", file.path(w, "geno.txt"), "-p", file.path(w, "yb.txt"), "-bslmm", "1", "-notsnp",
                           "-w", "50000", "-s", "200000", "-seed", 3300 + i, "-outdir", w, "-o", "bs"), stdout = FALSE, stderr = FALSE)
    P <- data.table::fread(file.path(w, "bs.param.txt"), data.table = FALSE)
    idx <- as.integer(sub("^s", "", P$rs)); o <- order(-P$gamma)
    list(sel = idx[o][P$gamma[o] > 0.5], top = idx[o[1:5]], top_pip = P$gamma[o[1:5]]) },
    error = function(e) conditionMessage(e))
  ## C. GridLMM first-step scan, each marker downdated from its own kernel
  gr <- tryCatch({ Xn <- X; rownames(Xn) <- smp; colnames(Xn) <- paste0("s", seq_len(m))
    Kn <- K; dimnames(Kn) <- list(smp, smp); dat <- data.frame(y = y, ID = smp)
    common <- list(formula = y ~ 1 + (1 | ID), test_formula = ~1, reduced_formula = ~0, data = dat, X = Xn, X_ID = "ID",
                   relmat = list(ID = list(K = Kn, p = m)), h2_step = 0.01, proximal_markers = as.list(seq_len(m)),
                   proximal_Xs = list(ID = Xn), mc.cores = 1, verbose = FALSE)
    gre <- do.call(GridLMM_GWAS, c(common, list(method = "REML")))$results
    gbf <- do.call(GridLMM_GWAS, c(common, list(method = "BF", inv_prior_X = c(0, 1 / 0.04))))$results
    lbf <- gbf$log_posterior_factor - gbf$log_posterior_factor_reduced
    list(top_reml = which.min(gre$p_value_REML), p_reml = min(gre$p_value_REML),
         top_bf = which.max(lbf), log10bf = max(lbf) / log(10), lbf = lbf) },
    error = function(e) conditionMessage(e))
  lead_j <- match(g$lead_pos, W$pos)
  res <- list(pop = pop, gene = g$symbol, tag = g$tag, n = n, m = m, pos = W$pos, lead_j = lead_j, pen = pen,
              gemma_paths = paths, gemma_sel = gsel, variants = if (length(reps)) do.call(rbind, reps) else NULL,
              mlmm = ml, qtltools = qt, farmcpu = fc, bslmm = bs, gridlmm = gr,
              elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(res, ck); unlink(w, recursive = TRUE)
  sprintf("done %s %s (%.0f s)", pop, g$symbol, res$elapsed)
}

POPSEL <- { a <- grep("^--pop=", args, value = TRUE); if (length(a)) strsplit(sub("^--pop=", "", a), ",")[[1]] else c("EUR", "YRI") }
jobs <- expand.grid(i = seq_len(nrow(GENES)), pop = POPSEL, stringsAsFactors = FALSE)
only <- grep("^--only=", args, value = TRUE)
if (length(only)) jobs <- jobs[GENES$symbol[jobs$i] %in% strsplit(sub("^--only=", "", only), ",")[[1]], , drop = FALSE]
st <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(run_gene(jobs$pop[k], jobs$i[k]),
  error = function(e) paste("FAIL", jobs$pop[k], GENES$symbol[jobs$i[k]], conditionMessage(e))),
  mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st))
writeLines(format(Sys.time()), file.path(OUT, "DONE"))
