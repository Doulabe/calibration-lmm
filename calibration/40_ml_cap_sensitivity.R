# ==============================================================================
# 40_ml_cap_sensitivity.R
# Sensitivity of the statistics to the upper end delta_max of the bounded domain Omega on
# which every maximization over delta is taken. The
# maximum-likelihood profile is unbounded (33h), so delta_max is part of the definition of
# GEMMA's likelihood-ratio and score statistics. Two parts, the grid resolution being held
# at 0.1 decade so that only the upper end changes:
#  app    the 38 variants entering GEMMA's paths (33e, with GEMMA's kernel and fixed
#         effects): observed Wald (REML), likelihood-ratio and score (ML) statistics, their
#         nominal p-values and the decision of each path at the extended-BIC level, for
#         delta_max in {1e3, 1e4, 1e5 (GEMMA), 1e6, 1e7, 1e8}, grid from 1e-5 as in GEMMA;
#  calib  the seven configurations of Supplementary Table S4 (39_mc_error.R): exact size of
#         the Wald test (REML), the likelihood ratio (ML, as GEMMA), the score with the ML
#         null estimate referred to F(1, nu - 1) (GEMMA's score test) and the score with the
#         REML null estimate referred to t_{nu-1} (row "GEMMA score test" of Table S5), for
#         delta_max in {1e3, 1e4, 1e5, 1e6, 1e8}, grid from 1e-4 as in the calibration study,
#         with the SAME 1,000 draws for every delta_max (paired comparison), and the share
#         of null datasets whose estimates sit at the upper end of the domain.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/40_ml_cap_sensitivity.R app
#   nice -n 10 Rscript calibration/40_ml_cap_sensitivity.R calib [cores=1]
# Output: results/ml_cap_app.csv; results/ml_cap_calib/<config>.rds, results/ml_cap_calib.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
MODE <- commandArgs(TRUE)[1]
stopifnot(MODE %in% c("app", "calib"))
RES <- "calibration/results"

if (MODE == "app") {
  ## ---- part 1: observed statistics and path decisions of the application --------------------------
  Sys.setenv(PC_SOURCE_ONLY = "1")
  try(source("calibration/33g_protected_check.R"), silent = TRUE)   # data, spectral_fixed, gemma_stats
  stopifnot(exists("spectral_fixed"), exists("gemma_stats"))
  CAPS <- c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
  V <- read.csv(file.path(RES, "geuvadis_cov_gemma_variants.csv"), stringsAsFactors = FALSE)
  out <- list()
  for (i in seq_len(nrow(V))) {
    v <- V[i, ]; g <- GENES[GENES$symbol == v$gene, ]
    R <- readRDS(file.path(COMP, sprintf("%s_%s.rds", v$pop, g$tag)))
    w <- R$variants[R$variants$pos == v$pos & R$variants$step == v$step, ]; stopifnot(nrow(w) == 1)
    before <- NULL
    for (tn in strsplit(w$tests, "+", fixed = TRUE)[[1]]) { idx <- R$gemma_paths[[tn]]$idx
      if (length(idx) >= w$step && idx[w$step] == w$j) { before <- idx[seq_len(w$step - 1L)]; break } }
    stopifnot(!is.null(before))
    W <- load_window(v$pop, g$tag); X <- scale(W$G); X[!is.finite(X)] <- 0
    y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); C <- covariates(v$pop, W$samples)
    n <- nrow(X); m <- ncol(X); K <- tcrossprod(X) / m; stopifnot(W$pos[w$j] == v$pos)
    pen <- log(n) + 2 * log(m); nu_lev <- n - (1L + ncol(C)) - v$step + 1L
    lw <- log10(2 * pt(-sqrt((nu_lev - 1) * expm1(pen / n)), nu_lev - 1))   # eBIC-equivalent level of 33e, all three paths
    lev <- c(lrt = lw, wald = lw, score = lw)
    for (cap in CAPS) {
      GRID2 <- c(0, 10^seq(-5, log10(cap), by = 0.1))     # global: used by spectral_fixed and gemma_stats
      sp <- spectral_fixed(X, y, C, before, w$j, K); nu <- sp$nu0
      Sk <- matrix(sp$obs$sk); o <- gemma_stats(sp, Sk, sp$obs$sp, sp$obs$ssr)
      ## position of the maximizers on the grid (upper end = last grid value)
      W2 <- 1 / (1 + outer(sp$ell, GRID2)); ld2 <- colSums(log(1 + outer(sp$ell, GRID2)))
      S2g <- colSums(sp$xi^2 * W2) + sp$xp^2
      A2 <- crossprod(Sk^2, W2) + (sp$obs$sp^2 + sp$obs$ssr); B2 <- crossprod(sp$xi * Sk, W2) + sp$xp * sp$obs$sp
      RSS2 <- A2 - sweep(B2^2, 2, S2g, "/"); ldm <- ld2 - sp$ldF
      g_m0 <- max.col(-0.5 * sweep(sp$n * log(A2), 2, ldm, "+"), ties.method = "first")
      g_m1 <- max.col(-0.5 * sweep(sp$n * log(RSS2), 2, ldm, "+"), ties.method = "first")
      g_r1 <- max.col(-0.5 * sweep((nu - 1) * log(RSS2), 2, ld2 + log(S2g), "+"), ties.method = "first")
      G <- length(GRID2)
      lp <- c(wald = (log(2) + pt(-abs(o$t), nu - 1, log.p = TRUE)) / log(10),
              score = pf(o$Fs, 1, nu - 1, lower.tail = FALSE, log.p = TRUE) / log(10),
              lrt = pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10))
      for (tn in strsplit(v$tests, "+", fixed = TRUE)[[1]])
        out[[length(out) + 1]] <- data.frame(pop = v$pop, gene = v$gene, pos = v$pos, step = v$step, path = tn,
          selected_gemma = v$selected, cap = cap, log10p = lp[[tn]], level = lev[[tn]], select = lp[[tn]] < lev[[tn]],
          gemma_log10p = v[[paste0("log10p_gemma_", tn)]], reml_dh1 = o$dh1, ml_dh0 = o$dm0, ml_dh1 = GRID2[g_m1],
          at_end_ml0 = g_m0 == G, at_end_ml1 = g_m1 == G, at_end_reml1 = g_r1 == G)
    }
    cat(sprintf("%s %s %d step %d done\n", v$pop, v$gene, v$pos, v$step))
  }
  A <- do.call(rbind, out)
  write.csv(A, file.path(RES, "ml_cap_app.csv"), row.names = FALSE)
  ## summary: decisions that differ from those at GEMMA's range (1e5)
  ref <- A[A$cap == 1e5, c("pop", "gene", "pos", "step", "path", "select", "log10p")]
  M <- merge(A, ref, by = c("pop", "gene", "pos", "step", "path"), suffixes = c("", "_1e5"))
  M$changed <- M$select != M$select_1e5; M$shift <- M$log10p - M$log10p_1e5
  print(aggregate(cbind(changed, at_end_ml0, at_end_ml1) ~ cap + path, M, sum))
  print(aggregate(shift ~ cap + path, M, function(z) c(med = median(z), max = max(abs(z)))), digits = 3)
  print(M[M$changed, c("pop", "gene", "pos", "step", "path", "cap", "log10p", "level", "log10p_1e5")], digits = 3)
  cat("check at 1e5, max |ours - GEMMA| by path:\n")
  print(tapply(abs(ref$log10p - A$gemma_log10p[A$cap == 1e5]), ref$path, max), digits = 3)
  quit(save = "no")
}

## ---- part 2: exact sizes on the configurations of Table S4 -------------------------------------------
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # functions, genes and configurations of 36
NC <- { a <- commandArgs(TRUE); if (length(a) >= 2 && grepl("^[0-9]+$", a[2])) as.integer(a[2]) else 1L }
CAPS <- c(1e3, 1e4, 1e5, 1e6, 1e8); B_CAP <- 1000L; NB <- 34L
## the seven configurations of Table S4, chosen as in 39_mc_error.R
CHOSEN <- c("geu_SNHG5_85909462_h0.05", "geu_ERAP2_96252589_h0.05", "syn_n1000_b5_h0.05", "syn_n200_b20_h0.05_unlinked")
C1 <- read.csv(file.path(RES, "calib_enlarged.csv"))
U <- unique(C1[grepl("^geu_", C1$config), c("config", "nu_delta")])
for (tg in c(25, 60, 150)) CHOSEN <- c(CHOSEN, U$config[which.min(abs(log(U$nu_delta / tg)))])
cfg_id <- function(z) if (z$kind == "syn") sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "") else
  sprintf("geu_%s_%d_h%.2f", GENES$symbol[z$i], z$pos, z$h2)
ids <- vapply(cfgs, cfg_id, ""); stopifnot(!anyNA(match(CHOSEN, ids)))
sel <- cfgs[match(CHOSEN, ids)]
OUTC <- file.path(RES, "ml_cap_calib"); dir.create(OUTC, showWarnings = FALSE, recursive = TRUE)

cap_config <- function(sp, delta0, id) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L; B <- B_CAP
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)                # the same draws for every delta_max
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s); Rk2 <- Rk * Rk
  u0  <- rnorm(B)                                             # one null dataset per draw, for the boundary shares
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  cF  <- sqrt(qf(1 - AL, 1, nu - 1) * nu / n)                 # GEMMA's score: n Q^2 / nu > F(1, nu - 1)
  res <- list(); pd <- list()
  for (cap in CAPS) {
    grid <- c(0, 10^seq(-4, log10(cap), by = 0.1)); G <- length(grid)
    W2 <- 1 / (1 + outer(ell, grid)); ld2 <- colSums(log(1 + outer(ell, grid))); S2g <- colSums(xi^2 * W2) + xp^2
    A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
    RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
    l1 <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+"); dh1 <- pick_grid(l1, grid)
    L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+"); ml1 <- apply(L1ml, 1, max)
    end_reml1 <- max.col(l1, ties.method = "first") == G; end_ml1 <- max.col(L1ml, ties.method = "first") == G
    rm(A_r, B_r, RSS_r, l1)
    at_r <- function(dl) { iD <- 1 / (1 + outer(ell, dl))
      A <- colSums(Rk2 * iD) + Rp^2 + ssr; Bx <- colSums(xi * Rk * iD) + xp * Rp; S2 <- colSums(xi^2 * iD) + xp^2
      list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2) }
    Wd <- at_r(dh1); sig_h <- sqrt(Wd$RSS / (nu - 1))
    wa <- sqrt(Wd$S2 / S0) / sig_h; wg <- (Wd$Bx / Wd$S2) * sqrt(Wd$S2) / sig_h
    ev <- function(u, what) {
      bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
      A2 <- crossprod(sk2, W2) + (s2p^2 + ssr); out <- list()
      if ("lrt" %in% what) out$lrt <- ml1 - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max)
      qq <- function(O) { dh0 <- pick_grid(O, grid); iD <- 1 / (1 + outer(ell, dh0))
        Ah <- colSums(sk2 * iD) + s2p^2 + ssr; Bh <- colSums(xi * sk * iD) + xp * s2p; Sh <- colSums(xi^2 * iD) + xp^2
        Bh / sqrt(Sh * Ah / nu) }
      if ("Q" %in% what) out$Q <- qq(-0.5 * sweep(nu * log(A2), 2, ld2, "+"))        # REML null estimate (as in 36)
      if ("Qml" %in% what) out$Qml <- qq(-0.5 * sweep(n * log(A2), 2, ld2, "+"))     # ML null estimate (GEMMA)
      if ("end" %in% what) out$end0 <- max.col(-sweep(n * log(A2), 2, ld2, "+"), ties.method = "first") == G
      out
    }
    lo <- rep(-30, B); hi <- rep(30, B)
    for (it in 1:40) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    split <- (lo + hi) / 2
    rootsU <- function(c0) {                                  # as in 36
      f <- function(u) ev(u, "lrt")$lrt
      always <- f(split) >= c0
      lo <- split; hi <- rep(120, B); never_r <- f(hi) < c0
      for (it in seq_len(NB)) { mid <- (lo + hi) / 2; up <- f(mid) < c0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      ur <- (lo + hi) / 2; ur[never_r] <- Inf
      lo <- rep(-120, B); hi <- split; never_l <- f(lo) < c0
      for (it in seq_len(NB)) { mid <- (lo + hi) / 2; dn <- f(mid) < c0; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
      ul <- (lo + hi) / 2; ul[never_l] <- -Inf
      ur[always] <- split[always]; ul[always] <- split[always]
      list(up = ur, dn = ul)
    }
    rootsQ <- function(key, c0) {                             # as in 36b
      one <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
        for (it in seq_len(NB + 6)) { mid <- (lo + hi) / 2; up <- ev(mid, key)[[key]] < level
          lo[up] <- mid[up]; hi[!up] <- mid[!up] }
        (lo + hi) / 2 }
      list(up = one(c0), dn = one(-c0))
    }
    lin <- function(a, g, c0) list(up = (c0 - g) / a, dn = (-c0 - g) / a)
    perdraw <- function(rt) pnorm(-rt$up) + pnorm(rt$dn)
    end0 <- ev(u0, "end")$end0
    for (k in seq_along(AL)) {
      P <- list(wald = perdraw(lin(wa, wg, tq_[k])), lrt_ml = perdraw(rootsU(cq[k])),
                score_ml = perdraw(rootsQ("Qml", cF[k])), score_reml = perdraw(rootsQ("Q", tq_[k])))
      for (mt in names(P)) {
        pd[[sprintf("%s_%g_%g", mt, cap, AL[k])]] <- P[[mt]]
        res[[length(res) + 1]] <- data.frame(config = id, nu_delta = NA, cap = cap, method = mt, alpha = AL[k],
          size_ratio = mean(P[[mt]]) / AL[k], se_log10 = sd(P[[mt]]) / (sqrt(B) * mean(P[[mt]]) * log(10)),
          end_ml0 = mean(end0), end_ml1 = mean(end_ml1), end_reml1 = mean(end_reml1))
      }
    }
  }
  list(summary = do.call(rbind, res), perdraw = pd)
}
run_cfg <- function(z) {
  id <- cfg_id(z); ck <- file.path(OUTC, paste0(id, ".rds")); if (file.exists(ck)) return(paste("cached", id))
  if (z$kind == "syn") {                                      # data as in 36 and 39
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    sp <- if (z$unlinked) { set.seed(77 + z$n); spectral_setup(Z, x = rnorm(z$n)) } else spectral_setup(Z, j = 25L)
  } else { W <- load_geu(z$i); sp <- spectral_setup(W$Z, j = match(z$pos, W$pos)) }
  set.seed(40000 + sum(utf8ToInt(id)))                        # fresh draws, independent of 36 and 39
  t0 <- Sys.time(); r <- cap_config(sp, z$h2 / (1 - z$h2), id)
  r$summary$nu_delta <- C1$nu_delta[match(id, C1$config)]
  r$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs")); saveRDS(r, ck)
  sprintf("done %s (%.0f s)", id, r$secs)
}
suppressPackageStartupMessages(library(parallel))
print(unlist(mclapply(sel, function(z) tryCatch(run_cfg(z), error = function(e) paste("FAIL", cfg_id(z), conditionMessage(e))),
                      mc.cores = NC, mc.preschedule = FALSE)))
S <- do.call(rbind, lapply(list.files(OUTC, "rds$", full.names = TRUE), function(f) readRDS(f)$summary))
write.csv(S, file.path(RES, "ml_cap_calib.csv"), row.names = FALSE)
S$l10 <- log10(pmax(S$size_ratio, 1e-300))
base <- S[S$cap == 1e4, c("config", "method", "alpha", "l10")]
M <- merge(S, base, by = c("config", "method", "alpha"), suffixes = c("", "_1e4")); M$d <- M$l10 - M$l10_1e4
options(width = 200)
print(M[order(M$nu_delta, M$method, M$alpha, M$cap), c("config", "nu_delta", "method", "alpha", "cap", "size_ratio", "d",
        "se_log10", "end_ml0", "end_ml1", "end_reml1")], digits = 3, row.names = FALSE)
