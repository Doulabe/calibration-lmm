# ==============================================================================
# 33i_protected_check_all.R
# Numerical control of the protected exact-law p-values for all 38 variants entering
# GEMMA's paths in the GEUVADIS application (Supplementary Table of the GEMMA variants),
# for the four statistics (GEMMA's Wald, likelihood-ratio and score tests, and the
# Kenward-Roger test), with the procedure of 33g_protected_check.R:
#   - delta_hat searched on GEMMA's range, [0, 1e-5 ... 1e5] (definition of the statistic);
#   - stage 1: p(delta0) at 0 and 10^-3, ..., 10^10 by 0.25 decade (4,000 draws);
#   - stage 2: 0.05-decade grid within 0.3 decade of the stage-1 maximum (10,000 draws);
#   - controlled value = largest p_hat + 2 se over both stages.
# Data, reduced model and Wald/KR laws come from 33g (sourced). For the score and the
# likelihood ratio, whose rejection regions are found by bisection, the statistics are
# evaluated by lean versions of gemma_stats that compute only what the statistic needs
# (same formulas, the grid quantities computed once per delta0, row maxima by max.col),
# and the two crossings are separated at the sign change of the effect estimate at the
# grid maximizer of the null likelihood instead of the sign of t; the lean and full
# versions give the same p-values on the same draws (check with --check).
# The six (variant, statistic) pairs already controlled by 33g are reused.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/33i_protected_check_all.R [cores=4] [--check] [--time]
# Output: results/protected_check_all/<pop>_<gene>_<pos>_<step>_<stat>.rds
# ==============================================================================
Sys.setenv(PC_SOURCE_ONLY = "1")
args0 <- commandArgs(TRUE)
try(source("calibration/33g_protected_check.R"), silent = TRUE)   # data, GRID2, D1, B1, B2, functions
stopifnot(exists("exact_one"), exists("setup_case"), length(GRID2) == 102)
NCORES <- if (length(args0) && grepl("^[0-9]+$", args0[1])) as.integer(args0[1]) else 4L
OUT_ALL <- "calibration/results/protected_check_all"; dir.create(OUT_ALL, showWarnings = FALSE, recursive = TRUE)

## ---- lean statistics for the score and the likelihood ratio ----------------------------------------
lean_eval <- function(sp, stat, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; n <- sp$n
  W2 <- 1 / (1 + outer(ell, GRID2)); ldm <- colSums(log(1 + outer(ell, GRID2))) - sp$ldF
  S2g <- colSums(xi^2 * W2) + xp^2; ix <- seq_len(B)
  if (stat == "score") function(Sk, Sp, Ssr) {
    A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr)
    ml0 <- -0.5 * sweep(n * log(A2), 2, ldm, "+"); dm0 <- pick_grid(ml0, GRID2)
    iD <- 1 / (1 + outer(ell, dm0))
    A <- colSums(Sk^2 * iD) + Sp^2 + Ssr; Bx <- colSums(xi * Sk * iD) + xp * Sp; S <- colSums(xi^2 * iD) + xp^2
    list(stat = n * Bx^2 / (S * A), sgn = Bx)
  } else function(Sk, Sp, Ssr) {
    A2 <- crossprod(Sk^2, W2) + (Sp^2 + Ssr); B2 <- crossprod(xi * Sk, W2) + xp * Sp
    RSS2 <- A2 - sweep(B2^2, 2, S2g, "/")
    ml0 <- -0.5 * sweep(n * log(A2), 2, ldm, "+"); ml1 <- -0.5 * sweep(n * log(RSS2), 2, ldm, "+")
    g0 <- max.col(ml0, ties.method = "first"); g1 <- max.col(ml1, ties.method = "first")
    list(stat = 2 * (ml1[cbind(ix, g1)] - ml0[cbind(ix, g0)]), sgn = B2[cbind(ix, g0)])
  }
}
exact_lean <- function(sp, delta0, level, stat, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  f0 <- lean_eval(sp, stat, B)
  ev <- function(u) { bv <- u / sqrt(S0); f0(Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  split <- { lo <- rep(-60, B); hi <- rep(60, B)
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- ev(mid)$sgn < 0; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
  ur <- { lo <- split; hi <- rep(250, B)
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; up <- ev(mid)$stat < level; lo[up] <- mid[up]; hi[!up] <- mid[!up] }; (lo + hi) / 2 }
  ul <- { lo <- rep(-250, B); hi <- split
    for (it in 1:NBIS) { mid <- (lo + hi) / 2; dn <- ev(mid)$stat < level; hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }; (lo + hi) / 2 }
  a1 <- pnorm(-ur, log.p = TRUE); a2 <- pnorm(ul, log.p = TRUE)
  mc_summary(pmax(a1, a2) + log1p(exp(-abs(a1 - a2))))
}
exact_any <- function(S, d0, stat, B)
  if (stat %in% c("wald", "kr")) exact_one(S$sp, d0, S$obs, stat, B) else
    exact_lean(S$sp, d0, if (stat == "score") S$obs$Fs else S$obs$lrt, stat, B)

## ---- cases: the 38 variants x 4 statistics -----------------------------------------------------------
V <- read.csv("calibration/results/geuvadis_cov_gemma_variants.csv", stringsAsFactors = FALSE)
CA <- do.call(rbind, lapply(seq_len(nrow(V)), function(i) data.frame(pop = V$pop[i], symbol = V$gene[i], pos = V$pos[i],
  step = V$step[i], stat = c("lrt", "score", "wald", "kr"), stringsAsFactors = FALSE)))
CA$id <- sprintf("%s_%s_%d_%d_%s", CA$pop, gsub("[^A-Za-z0-9]", "", CA$symbol), CA$pos, CA$step, CA$stat)
cost <- with(CA, ifelse(stat %in% c("lrt", "score"), ifelse(pop == "EUR", 40, 10), ifelse(pop == "EUR", 2, 0.5)))
CA <- CA[order(-cost), ]
OLD <- do.call(rbind, lapply(list.files("calibration/results/protected_check", "_(wald|kr|score)[.]rds$", full.names = TRUE),
  function(f) { r <- readRDS(f); cbind(r$info[, c("pop", "gene", "pos", "step", "stat")], file = f) }))

setup_case2 <- function(pop, symbol, pos, step, stat) {    # setup_case of 33g, with the step
  g <- GENES[GENES$symbol == symbol, ]
  R <- readRDS(file.path(COMP, sprintf("%s_%s.rds", pop, g$tag)))
  v <- R$variants[R$variants$pos == pos & R$variants$step == step, ]; stopifnot(nrow(v) == 1)
  before <- NULL
  for (tn in strsplit(v$tests, "+", fixed = TRUE)[[1]]) { idx <- R$gemma_paths[[tn]]$idx
    if (length(idx) >= v$step && idx[v$step] == v$j) { before <- idx[seq_len(v$step - 1L)]; break } }
  stopifnot(!is.null(before))
  W <- load_window(pop, g$tag); X <- scale(W$G); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); C <- covariates(pop, W$samples)
  n <- nrow(X); m <- ncol(X); K <- tcrossprod(X) / m; stopifnot(W$pos[v$j] == pos)
  sp <- spectral_fixed(X, y, C, before, v$j, K)
  o <- gemma_stats(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
  pen <- log(n) + 2 * log(m); ncf <- 1L + ncol(C); nu_lev <- n - ncf - v$step + 1L
  level <- if (stat == "lrt") pchisq(pen, 1, lower.tail = FALSE, log.p = TRUE) / log(10) else
    log10(2 * pt(-sqrt((nu_lev - 1) * expm1(pen / n)), nu_lev - 1))
  stored <- v[[paste0("log10p_sup_", stat)]]; argst <- v[[paste0("argsup_", stat)]]
  nominal <- switch(stat, wald = v$log10p_gemma_wald, kr = v$log10p_kr_full, score = v$log10p_gemma_score, lrt = v$log10p_gemma_lrt)
  ours <- switch(stat, wald = (log(2) + pt(-abs(o$t), sp$nu0 - 1, log.p = TRUE)) / log(10), kr = NA,
                 score = pf(o$Fs, 1, sp$nu0 - 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                 lrt = pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10))
  list(sp = sp, obs = list(t = o$t, Fs = o$Fs, lrt = o$lrt, p_kr = 10^v$log10p_kr_full),
       info = data.frame(pop = pop, gene = symbol, pos = pos, step = step, stat = stat, n = n, m = m, nu = sp$nu0,
                         nu_j = v$nu_j, level = level, nominal = nominal, ours_nominal = ours,
                         stored_protected = stored, stored_argmax = argst))
}
run_one <- function(k) {
  cs <- CA[k, ]; ck <- file.path(OUT_ALL, paste0(cs$id, ".rds"))
  if (file.exists(ck)) return(paste("cached", cs$id))
  old <- if (!is.null(OLD)) OLD[OLD$pop == cs$pop & OLD$gene == cs$symbol & OLD$pos == cs$pos & OLD$step == cs$step &
                                OLD$stat == cs$stat, ] else NULL
  if (!is.null(old) && nrow(old) == 1) {                    # reuse the controlled value of 33g
    r <- readRDS(as.character(old$file)); S <- setup_case2(cs$pop, cs$symbol, cs$pos, cs$step, cs$stat)
    i <- which.max(r$profile$log10p)
    res <- list(info = cbind(S$info, new_max = r$profile$log10p[i], new_argmax = r$profile$delta0[i],
                             controlled = max(r$profile$log10up), at_1e4 = r$info$at_1e4, at_1e10 = r$info$at_1e10,
                             elapsed_min = r$info$elapsed_min, source = "33g"), profile = r$profile)
    saveRDS(res, ck); return(paste("reused", cs$id))
  }
  set.seed(4000 + k); t0 <- Sys.time()
  S <- setup_case2(cs$pop, cs$symbol, cs$pos, cs$step, cs$stat)
  prof1 <- t(vapply(D1, function(d0) exact_any(S, d0, cs$stat, B1), numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, function(d0) exact_any(S, d0, cs$stat, B2), numeric(3)))
  prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2)); prof$id <- cs$id
  i <- which.max(prof$log10p)
  res <- list(info = cbind(S$info, new_max = prof$log10p[i], new_argmax = prof$delta0[i], controlled = max(prof$log10up),
                           at_1e4 = prof1[D1 == 1e4, "log10p"], at_1e10 = prof1[D1 == 1e10, "log10p"],
                           elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")), source = "33i"),
              profile = prof)
  saveRDS(res, ck)
  sprintf("done %s (%.1f min)", cs$id, res$info$elapsed_min)
}

if ("--check" %in% args0) {                                 # lean versus full evaluator on the same draws
  for (st in c("score", "lrt")) {
    S <- setup_case2("YRI", "TRA2A-AS", 23134790L, 1L, st)
    for (d0 in c(1, 300)) {
      set.seed(9); e1 <- exact_lean(S$sp, d0, if (st == "score") S$obs$Fs else S$obs$lrt, st, 2000L)
      set.seed(9); e2 <- exact_one(S$sp, d0, S$obs, st, 2000L)
      cat(sprintf("%-5s delta0=%-4g lean %.4f  full %.4f\n", st, d0, e1[["log10p"]], e2[["log10p"]]))
    }
  }
  quit(save = "no")
}
if ("--time" %in% args0) {
  for (k in which(!duplicated(paste(CA$pop, CA$stat)))) { cs <- CA[k, ]
    S <- setup_case2(cs$pop, cs$symbol, cs$pos, cs$step, cs$stat)
    tm <- system.time(e <- exact_any(S, 1, cs$stat, B1))[["elapsed"]]
    cat(sprintf("%s %s %s: %.1f s per evaluation (B1), nominal %.2f ours %.2f\n", cs$pop, cs$stat, cs$id, tm,
                S$info$nominal, S$info$ours_nominal)) }
  quit(save = "no")
}
st <- parallel::mclapply(seq_len(nrow(CA)), function(k) tryCatch(run_one(k), error = function(e) paste("FAIL", CA$id[k], conditionMessage(e))),
                         mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st))
R <- lapply(list.files(OUT_ALL, "[.]rds$", full.names = TRUE), readRDS)
write.csv(do.call(rbind, lapply(R, `[[`, "info")), "calibration/results/protected_check_all_summary.csv", row.names = FALSE)
cat("cases:", length(R), "of", nrow(CA), "\n")
