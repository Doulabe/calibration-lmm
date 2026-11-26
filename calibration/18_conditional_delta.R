# ==============================================================================
# 18_conditional_delta.R
#
# Test of the conditional-delta lever proposed in the note (section sec:selection).
# At the start of the stepwise path the regional kernel absorbs the causal variants
# not yet selected, so delta_hat is inflated by (h2bg + K h2c) / h2bg and the first
# threshold is not crossed (17_delta_inflation_check.R). The lever estimates delta
# conditionally on a provisional set of signals. Its risk is the one of the Wald
# with delta_hat under H1: a delta estimated partly under the alternative.
#
# Datasets: the 1,000 datasets of 14_selection_demo.R (regenerated from the same
# seeds, checked against the stored causal sets) plus 800 extra global-null
# datasets (b = 26..125) for the calibration curves.
#
# Provisional sets: S_lib = the LRT path (rv below) down to the liberal level
# pen_lib = qchisq(1 - 1e-3, 1); S_sus = the SuSiE lead set (top PIP per 95% CS,
# L = 10), a sparse fit without background.
# Procedures (same stepwise path; each path is recorded down to pen_lib, so the
# selection at any stricter threshold is a prefix of the recorded path):
#   rv      LRT, delta re-estimated under each model
#   oracle  delta fixed at its design value delta* = h2bg / sigma2_e
#   c1      delta_tilde = REML with S_lib as fixed effects, then fixed along the path
#   c2      same with S_sus
#   c3      conditional LRT: delta re-estimated under each model, with the provisional
#           signals of S_lib that are neither in LD with the candidate (r2 < 0.1) nor
#           with the selected set (r2 < 0.5) as extra fixed effects
#   c4      same with S_sus
#   c5      fixed point: delta_tilde given the current selection, delta-fixed path,
#           iterated from S_lib until the selection is stable (run at each threshold)
# Thresholds (LRT scale): qchisq(1 - alpha, 1), alpha in 1e-3 .. 1e-7, and the eBIC
# penalty log n + 2 log m of the demonstration.
#
# Usage (from the repository root):
#   nice -n 15 Rscript calibration/18_conditional_delta.R [cores=5] [--smoke]
# Output: results/conditional_delta/<id>.rds and, once all datasets are done,
#   results/conditional_delta_{sel,pertest,meta}.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 5L
SMOKE  <- "--smoke" %in% args
suppressPackageStartupMessages({ library(parallel); library(susieR) })

B <- 25L; B_NULL <- 125L; K_MAX <- 10L
GRIDD   <- c(0, 10^seq(-4, 3, length.out = 41))         # delta nodes of the LRT (as in 14)
GR      <- c(0, 10^seq(-4, 3, length.out = 281))        # finer grid for delta_tilde (as in 17)
PEN_LIB <- qchisq(1 - 1e-3, 1)
ALPHAS  <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-7)
R2_LD <- 0.1; R2_SEL <- 0.5
PROCS <- c("rv", "oracle", "c1", "c2", "c3", "c4")

## exact per-candidate pieces at one delta, identical to pieces_x() of 14_selection_demo.R
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

## profile-ML LRT of each candidate, delta re-estimated under each model (lrt of 14)
lrt_grid <- function(y, Xc, F_k, U, s_raw, m_k) {
  n <- length(y); k <- ncol(F_k); G <- length(GRIDD); s_rank <- s_raw / m_k
  Uty <- as.numeric(crossprod(U, y)); UtF <- crossprod(U, F_k); UtX <- crossprod(U, Xc)
  A0 <- A1 <- matrix(0, ncol(Xc), G)
  for (g in seq_len(G)) {
    pc <- pieces_x(Uty, UtF, UtX, s_rank, GRIDD[g], GRIDD[g] / m_k, n, k)
    A0[, g] <- -pc$ldZ - n * log(pc$RSS0)
    A1[, g] <- -pc$ldZ - n * log(pmax(pc$RSS0 - pc$u^2 / pc$D, 1e-300))
  }
  apply(A1, 1, max) - apply(A0, 1, max)
}
## LRT of each candidate at a fixed delta (oracle arm of 16)
lrt_fixed <- function(y, Xc, F_k, U, s_raw, m_k, delta) {
  n <- length(y)
  pc <- pieces_x(as.numeric(crossprod(U, y)), crossprod(U, F_k), crossprod(U, Xc), s_raw / m_k,
                 delta, delta / m_k, n, ncol(F_k))
  -n * log(pmax(1 - pc$u^2 / (pc$D * pc$RSS0), 1e-300))
}

## all procedures on one dataset
analyse <- function(y, X, delta_star, S_sus, pertest = FALSE) {
  n <- nrow(X); m <- ncol(X); cache <- new.env()
  eig <- function(ex) {                                  # kernel of the window minus ex
    key <- paste0("x:", paste(sort(ex), collapse = ","))
    if (is.null(cache[[key]])) { keep <- setdiff(seq_len(m), ex)
      e <- eigen(tcrossprod(X[, keep, drop = FALSE]), symmetric = TRUE)
      cache[[key]] <- list(U = e$vectors, s = pmax(e$values, 0), keep = keep) }
    cache[[key]] }
  Fm <- function(ix) cbind(1, X[, ix, drop = FALSE])
  r2 <- function(a, b) (crossprod(X[, a, drop = FALSE], X[, b, drop = FALSE]) / (n - 1))^2
  rv_stat <- function(sel) { e <- eig(sel)
    list(idx = e$keep, stat = lrt_grid(y, X[, e$keep, drop = FALSE], Fm(sel), e$U, e$s, length(e$keep) - 1L),
         df2 = rep(n - length(sel) - 2, length(e$keep))) }
  fixed_stat <- function(delta) function(sel) { e <- eig(sel)
    list(idx = e$keep, stat = lrt_fixed(y, X[, e$keep, drop = FALSE], Fm(sel), e$U, e$s, length(e$keep) - 1L, delta),
         df2 = rep(n - length(sel) - 2, length(e$keep))) }
  cond_stat <- function(S) function(sel) {
    P <- setdiff(S, sel)
    if (length(P) && length(sel)) P <- P[apply(r2(P, sel), 1, max) < R2_SEL]
    if (!length(P)) return(rv_stat(sel))
    cand <- setdiff(seq_len(m), sel); ld <- r2(cand, P) >= R2_LD
    key <- apply(ld, 1, function(v) paste(which(v), collapse = ","))
    out <- lapply(split(seq_along(cand), key), function(ix) {
      Pg <- P[!ld[ix[1], ]]; cg <- cand[ix]; e <- eig(c(sel, Pg))
      list(idx = cg, stat = lrt_grid(y, X[, cg, drop = FALSE], Fm(c(sel, Pg)), e$U, e$s, length(e$keep) - 1L),
           df2 = rep(n - length(sel) - length(Pg) - 2, length(cg))) })
    list(idx = unlist(lapply(out, `[[`, "idx"), use.names = FALSE),
         stat = unlist(lapply(out, `[[`, "stat"), use.names = FALSE),
         df2 = unlist(lapply(out, `[[`, "df2"), use.names = FALSE))
  }
  run_path <- function(fun) {
    sel <- integer(0); rec <- numeric(0); first <- NULL
    for (k in seq_len(K_MAX)) {
      st <- fun(sel); if (k == 1L) first <- st
      b <- which.max(st$stat); if (!length(b) || !(st$stat[b] > PEN_LIB)) break
      sel <- c(sel, st$idx[b]); rec <- c(rec, st$stat[b])
    }
    list(sel = sel, stat = rec, first = first)
  }
  delta_reml <- function(S) {                            # REML, fixed effects [1, X_S], kernel minus S
    e <- eig(S); ev <- e$s / length(e$keep); F <- Fm(S); k <- ncol(F)
    Uty <- as.numeric(crossprod(e$U, y)); UtF <- crossprod(e$U, F)
    ll <- vapply(GR, function(d) { w <- 1 / (1 + d * ev); M <- crossprod(UtF, w * UtF)
      p <- crossprod(UtF, w * Uty); R <- sum(w * Uty^2) - sum(p * solve(M + 1e-10 * diag(k), p))
      -0.5 * sum(log1p(d * ev)) - 0.5 * as.numeric(determinant(M)$modulus) - 0.5 * (n - k) * log(R) }, 0)
    GR[which.max(ll)]
  }
  at <- function(p, thr) { k <- which(p$stat <= thr); p$sel[seq_len(if (length(k)) k[1] - 1L else length(p$sel))] }

  paths <- list()
  paths$rv <- run_path(rv_stat); S_lib <- paths$rv$sel
  paths$oracle <- run_path(fixed_stat(delta_star))
  d_null <- delta_reml(integer(0)); d_lib <- delta_reml(S_lib); d_sus <- delta_reml(S_sus)
  fcache <- list()
  fpath <- function(d) { key <- sprintf("%.8g", d)
    if (is.null(fcache[[key]])) fcache[[key]] <<- run_path(fixed_stat(d)); fcache[[key]] }
  paths$c1 <- fpath(d_lib); paths$c2 <- fpath(d_sus)
  paths$c3 <- if (length(S_lib)) run_path(cond_stat(S_lib)) else paths$rv
  paths$c4 <- if (length(S_sus)) run_path(cond_stat(S_sus)) else paths$rv

  thr <- c(setNames(qchisq(1 - ALPHAS, 1), paste0("a", format(ALPHAS, scientific = TRUE))),
           ebic = log(n) + 2 * log(m))
  sels <- lapply(paths, function(p) lapply(thr, function(t) at(p, t)))
  fix <- lapply(thr, function(t) {                       # c5, fixed point at threshold t
    S <- S_lib; seen <- character(0)
    for (it in 1:8) { d <- delta_reml(S); sel <- at(fpath(d), t)
      key <- paste(sort(sel), collapse = ",")
      if (setequal(sel, S) || key %in% seen) break
      seen <- c(seen, key); S <- sel }
    list(sel = sel, delta = d, iter = it) })
  sels$c5 <- lapply(fix, `[[`, "sel")

  max1 <- vapply(paths, function(p) max(p$first$stat), 0)
  pt <- NULL
  if (pertest) pt <- do.call(rbind, lapply(names(paths), function(nm) {
    f <- paths[[nm]]$first; t2 <- f$df2 * expm1(f$stat / n); p <- pf(t2, 1, f$df2, lower.tail = FALSE)
    data.frame(proc = nm, n_cand = length(p), rej_1e2 = sum(p < 1e-2), rej_1e3 = sum(p < 1e-3),
               rej_1e4 = sum(p < 1e-4), max1 = max(f$stat)) }))
  list(thr = thr, sels = sels, rec = lapply(paths, function(p) p[c("sel", "stat")]),
       S_lib = S_lib, S_sus = S_sus, delta = c(star = delta_star, null = d_null, lib = d_lib, sus = d_sus),
       fix_delta = vapply(fix, `[[`, 0, "delta"), fix_iter = vapply(fix, `[[`, 0L, "iter"),
       max1 = max1, pertest = pt)
}

metrics <- function(sel, causal, X) {                    # identical to 14_selection_demo.R
  kc <- length(causal); ks <- length(sel)
  if (ks == 0) return(c(K_hat = 0, tp = 0, fp08 = 0, fp025 = 0, recall = if (kc) 0 else NA,
                        precision = NA, f1 = if (kc) 0 else NA, any_off = 0, exact = 0))
  if (kc == 0) return(c(K_hat = ks, tp = 0, fp08 = ks, fp025 = ks, recall = NA, precision = 0,
                        f1 = NA, any_off = 1, exact = 0))
  r2 <- cor(X[, sel, drop = FALSE], X[, causal, drop = FALSE])^2
  best <- apply(r2, 1, max); tagged <- apply(r2 >= 0.8, 2, any)
  tp <- sum(best >= 0.8); prec <- tp / ks; rec <- mean(tagged)
  c(K_hat = ks, tp = tp, fp08 = sum(best < 0.8), fp025 = sum(best < 0.25), recall = rec,
    precision = prec, f1 = if (prec + rec > 0) 2 * prec * rec / (prec + rec) else 0,
    any_off = as.numeric(any(best < 0.25)), exact = sum(sel %in% causal))
}

## ---- data and design (as in 14_selection_demo.R, plus extra global-null datasets) ------
GENES <- data.frame(tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000124587"),
                    symbol = c("ERAP2", "RPS26", "SLFN5", "PEX6-region"), stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_Z <- function(tag) {
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(tag, ".gt"))), "\t", fixed = TRUE)
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  Z <- scale(G[, keep, drop = FALSE]); Z[!is.finite(Z)] <- 0; Z
}
ZS <- lapply(GENES$tag, load_Z)
grid <- rbind(
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = c(1L, 3L), h2c = c(0.08, 0.15),
              b = seq_len(B), stringsAsFactors = FALSE),
  expand.grid(gene = seq_len(nrow(GENES)), h2bg = c(0.10, 0.25), Kc = 0L, h2c = 0,
              b = (B + 1L):B_NULL, stringsAsFactors = FALSE))
grid$cell <- sprintf("%s_bg%.2f_K%d_c%.2f", GENES$symbol[grid$gene], grid$h2bg, grid$Kc, grid$h2c)
grid$id <- sprintf("%s_b%03d", grid$cell, grid$b)
IN  <- "calibration/results/selection_demo"
OUT <- "calibration/results/conditional_delta"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

one_rep <- function(i) {
  g <- grid[i, ]; ck <- file.path(OUT, paste0(g$id, ".rds"))
  if (!SMOKE && file.exists(ck) && file.size(ck) > 0) return("cached")
  Z <- ZS[[g$gene]]; n <- nrow(Z); m <- ncol(Z)
  set.seed(30000 + 1000 * g$gene + 100 * g$Kc + round(100 * g$h2bg) + 7919 * g$b + round(1000 * g$h2c))
  causal <- integer(0)
  while (length(causal) < g$Kc) {
    cand <- sample.int(m, 1)
    if (!length(causal) || max(cor(Z[, cand], Z[, causal, drop = FALSE])^2) < 0.5) causal <- c(causal, cand)
  }
  bvec <- rnorm(m); u <- as.numeric(Z %*% bvec); u <- u / sd(u) * sqrt(g$h2bg)
  sig <- as.numeric(Z[, causal, drop = FALSE] %*% rep(sqrt(g$h2c), g$Kc))
  s2e <- max(1 - g$h2bg - g$Kc * g$h2c, 0.05)
  y <- sig + u + rnorm(n) * sqrt(s2e)
  f14 <- file.path(IN, paste0(g$id, ".rds")); lrt14 <- NULL; S_sus <- NULL
  if (file.exists(f14)) {                                # dataset of the demonstration
    st <- readRDS(f14)
    if (!identical(as.integer(st$causal), as.integer(causal))) stop("causal set differs from 14")
    lrt14 <- st$selections$lrt; S_sus <- st$selections$susie
  }
  if (is.null(S_sus) || (length(S_sus) == 1 && is.na(S_sus))) {
    f <- susie(Z, y, L = 10, verbose = FALSE); cs <- f$sets$cs
    S_sus <- if (is.null(cs)) integer(0) else unique(vapply(cs, function(ix) ix[which.max(f$pip[ix])], 1L))
  }
  t0 <- Sys.time()
  a <- analyse(y, Z, g$h2bg / s2e, as.integer(S_sus), pertest = g$Kc == 0)
  res <- c(list(id = g$id, gene = GENES$symbol[g$gene], h2bg = g$h2bg, Kc = g$Kc, h2c = g$h2c, b = g$b,
                causal = causal, lrt14 = lrt14,
                check_rv = if (is.null(lrt14)) NA else identical(sort(as.integer(a$sels$rv$ebic)), sort(as.integer(lrt14))),
                elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))), a)
  if (SMOKE) return(res)
  saveRDS(res, ck); "done"
}

collect <- function() {
  fs <- list.files(OUT, "rds$", full.names = TRUE)
  rows <- lapply(fs, function(f) {
    r <- readRDS(f); Z <- ZS[[match(r$gene, GENES$symbol)]]
    sel <- do.call(rbind, lapply(names(r$sels), function(p) do.call(rbind, lapply(names(r$thr), function(t)
      data.frame(id = r$id, gene = r$gene, h2bg = r$h2bg, Kc = r$Kc, h2c = r$h2c, b = r$b, proc = p,
                 thr_name = t, thr = r$thr[[t]], t(metrics(r$sels[[p]][[t]], r$causal, Z)))))))
    meta <- data.frame(id = r$id, gene = r$gene, h2bg = r$h2bg, Kc = r$Kc, h2c = r$h2c, b = r$b,
                       delta_star = r$delta[["star"]], delta_null = r$delta[["null"]],
                       delta_lib = r$delta[["lib"]], delta_sus = r$delta[["sus"]],
                       n_lib = length(r$S_lib), n_sus = length(r$S_sus), check_rv = r$check_rv, elapsed = r$elapsed,
                       t(setNames(r$max1, paste0("max1_", names(r$max1)))),
                       t(setNames(r$fix_iter, paste0("fixiter_", names(r$fix_iter)))))
    pt <- if (is.null(r$pertest)) NULL else cbind(data.frame(id = r$id, gene = r$gene, h2bg = r$h2bg, b = r$b), r$pertest)
    list(sel = sel, meta = meta, pt = pt) })
  write.csv(do.call(rbind, lapply(rows, `[[`, "sel")), "calibration/results/conditional_delta_sel.csv", row.names = FALSE)
  write.csv(do.call(rbind, lapply(rows, `[[`, "meta")), "calibration/results/conditional_delta_meta.csv", row.names = FALSE)
  write.csv(do.call(rbind, lapply(rows, `[[`, "pt")), "calibration/results/conditional_delta_pertest.csv", row.names = FALSE)
  length(fs)
}

if (SMOKE) {
  ids <- sub("^--id=", "", grep("^--id=", args, value = TRUE))
  smoke_i <- if (length(ids)) match(ids, grid$id) else
    c(which(grid$gene == 1 & grid$h2bg == 0.10 & grid$Kc == 3 & grid$h2c == 0.15 & grid$b == 1),
      which(grid$gene == 2 & grid$h2bg == 0.25 & grid$Kc == 0 & grid$b == 30))
  for (i in smoke_i) {
    r <- one_rep(i)
    cat(sprintf("\n%s  causal: %s  (%.1f s)  check_rv vs 14: %s\n", r$id, paste(r$causal, collapse = ","),
                r$elapsed, r$check_rv))
    cat(sprintf("  delta* %.3f  delta_null %.3f  delta_lib %.3f (S_lib: %s)  delta_sus %.3f (S_sus: %s)\n",
                r$delta[["star"]], r$delta[["null"]], r$delta[["lib"]], paste(r$S_lib, collapse = ","),
                r$delta[["sus"]], paste(r$S_sus, collapse = ",")))
    for (p in names(r$sels)) cat(sprintf("  %-7s ebic: %-22s a1e-3: %s\n", p, paste(r$sels[[p]]$ebic, collapse = ","),
                                         paste(r$sels[[p]][[1]], collapse = ",")))
    if (!is.null(r$pertest)) print(r$pertest, row.names = FALSE)
  }
} else {
  ord <- order(grid$b, grid$gene, grid$h2bg, grid$Kc, grid$h2c)
  st <- unlist(mclapply(ord, function(i) tryCatch(one_rep(i), error = function(e)
    paste("FAIL", grid$id[i], conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE))
  cat("datasets:", sum(st == "done"), "new,", sum(st == "cached"), "cached,", sum(grepl("FAIL", st)), "failed\n")
  print(head(st[grepl("FAIL", st)]))
  if (!any(grepl("FAIL", st))) cat("collected:", collect(), "datasets\n")
}
