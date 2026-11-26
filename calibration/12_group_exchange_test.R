# ==============================================================================
# 12_group_exchange_test.R
#
# Lighter conditioning than the sign-flip CLR (10_clr_test.R). Under H0 the spectral
# coordinates in a group of directions with (nearly) equal eigenvalues are (nearly)
# iid N(0, sigma2 (1 + delta l)), hence exchangeable and symmetric. The reference
# distribution permutes the values within each group and flips their signs, i.e.
# conditions only on the multiset of |values| of each group, not on each |s_i|.
# The test is exact when the eigenvalues of a group are equal (always true for the
# group of zero eigenvalues: the x_perp direction and the directions outside the
# kernel), approximately similar otherwise.
#
# Groups: kernel eigenvalues sorted, consecutive runs with l_max / l_min <= ratio;
# ratio = 1 gives singletons (= the sign-flip CLR when there is no zero group).
# Statistic: the same profile-ML LRT as 10_clr_test.R (both profiles recomputed,
# since A(delta) changes when values move between unequal eigenvalues).
# Tails by importance sampling: sequential tilted assignment within each group,
# value k chosen at position i with probability prop. to cosh(theta |c_i| |v_k|),
# sign tilted by exp(+-theta c_i v), symmetric mixture of +-theta; c_i = w_i(d*) xi_i
# at the delta d* maximising the observed H1 profile; theta from the with-replacement
# saddlepoint of sum c_i v_i. Weights: uniform target / proposal.
#
# Same configurations, H1 simulation and comparators as 10_clr_test.R, but the
# datasets are shared across groupings (paired).
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/12_group_exchange_test.R [cores=5] [--smoke]
# Output: results/group_exchange/<config>.rds, results/group_exchange.csv
# ==============================================================================

args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 5L
SMOKE  <- "--smoke" %in% args
src <- readLines("calibration/10_clr_test.R")
cut <- grep("^## ---- configurations", src)
eval(parse(text = src[1:(cut - 1)]))           # helpers of 06, 08 and 10
RATIOS <- c(1, 1.25, 2, 4)
N0  <- if (SMOKE) 150L else 600L
N1  <- if (SMOKE) 80L else 300L
MIS <- if (SMOKE) 1000L else 1500L
REDUCED <- "--reduced" %in% args; TAIL <- "--tail_size" %in% args
if (REDUCED) { RATIOS <- c(2, 4); N0 <- 400L; N1 <- 200L; MIS <- 1000L }
if (TAIL) { RATIOS <- 4; N0 <- 3000L; N1 <- 0L; MIS <- 1000L }

logcosh <- function(x) { ax <- abs(x); ax + log1p(exp(-2 * ax)) - log(2) }

make_groups <- function(ell, ratio) {           # kernel directions only
  o <- order(ell, decreasing = TRUE); g <- integer(length(ell)); gid <- 0L; top <- NA_real_
  for (idx in o) {
    if (is.na(top) || top / ell[idx] > ratio * (1 + 1e-12)) { gid <- gid + 1L; top <- ell[idx] }
    g[idx] <- gid
  }
  split(seq_along(ell), g)
}

## profile-ML LRT for signed values V (M x K: kernel directions 1..r, then the x_perp
## slot if any). A(delta) uses the kernel directions plus the whole zero group, whose total
## sum of squares SSz is invariant under exchange; B(delta) uses every slot.
lrt_of <- function(V, ellK, xiK, SSz, n, r) {
  W <- 1 / (1 + outer(ellK, GRIDL)); ld <- colSums(log(1 + outer(ellK[seq_len(r)], GRIDL)))
  S <- colSums(xiK^2 * W)
  A <- (V[, seq_len(r), drop = FALSE]^2) %*% W[seq_len(r), , drop = FALSE] + SSz
  B <- V %*% (xiK * W)
  h1 <- sweep(-n * log(pmax(A - sweep(B^2, 2, S, "/"), 1e-300)), 2, ld, "-")
  h0 <- sweep(-n * log(A), 2, ld, "-")
  apply(h1, 1, max) - apply(h0, 1, max)
}
rowmax <- function(L) L[cbind(seq_len(nrow(L)), max.col(L, ties.method = "first"))]

## log p-value of the group-exchange test by importance sampling.
## v_obs: signed values of the K slots; zero_vals: values of the zero group (x_perp slot
## value first when present); groups: list of kernel-direction index sets.
ge_logp <- function(v_obs, zero_vals, ellK, xiK, groups, has_p, n, r, M = MIS, TMAX = 40L) {
  K <- length(v_obs); SSz <- sum(zero_vals^2)
  lam_obs <- lrt_of(matrix(v_obs, 1), ellK, xiK, SSz, n, r)
  W <- 1 / (1 + outer(ellK, GRIDL)); ld <- colSums(log(1 + outer(ellK[seq_len(r)], GRIDL)))
  S <- colSums(xiK^2 * W)
  Ao <- colSums(v_obs[seq_len(r)]^2 * W[seq_len(r), , drop = FALSE]) + SSz
  Bo <- colSums(v_obs * xiK * W)
  gs <- which.max(-ld - n * log(Ao - Bo^2 / S))
  cvec <- xiK * W[, gs]; tobs <- abs(Bo[gs])
  pools <- lapply(groups, function(ix) list(pos = ix, vals = abs(v_obs[ix]), zero = FALSE))
  if (has_p) pools[[length(pools) + 1]] <- list(pos = K, vals = abs(zero_vals), zero = TRUE)
  ## tilted slots: the TMAX largest |c_i| overall (and every zero-group slot)
  tilt_set <- head(order(-abs(cvec)), TMAX)
  ## theta: with-replacement saddlepoint of sum c_i v_i at |T_obs|
  Kp <- function(th) sum(unlist(lapply(pools, function(p) vapply(p$pos, function(i) {
    x <- th * abs(cvec[i]) * p$vals; lc <- logcosh(x); w <- exp(lc - max(lc))
    sum(abs(cvec[i]) * p$vals * tanh(x) * w) / sum(w) }, 0))))
  cmax <- max(abs(cvec)) * max(unlist(lapply(pools, `[[`, "vals")))
  lo <- 0; hi <- 30 / max(cmax, 1e-12)
  if (Kp(hi) < tobs) th <- hi else { for (it in 1:40) { mid <- (lo + hi) / 2
    if (Kp(mid) < tobs) lo <- mid else hi <- mid }; th <- (lo + hi) / 2 }
  comp <- rep(c(1, -1), length.out = M)                 # mixture component of each sample
  V <- matrix(0, M, K); lpval <- numeric(M); lsp <- numeric(M); lsm <- numeric(M); lp0 <- 0
  for (p in pools) {
    g <- length(p$vals)
    tpos <- p$pos[p$pos %in% tilt_set | p$zero]; tpos <- tpos[order(-abs(cvec[tpos]))]
    upos <- setdiff(p$pos, tpos)
    lp0 <- lp0 - (if (p$zero) log(g) else lgamma(g + 1))
    avail <- matrix(TRUE, M, g)
    for (i in tpos) {                                   # tilted sequential choices (Gumbel-max)
      L <- matrix(rep(logcosh(th * abs(cvec[i]) * p$vals), each = M), M, g); L[!avail] <- -Inf
      mx <- rowmax(L); lse <- mx + log(rowSums(exp(L - mx)))
      k <- max.col(L - log(-log(matrix(runif(M * g), M, g))), ties.method = "first")
      lpval <- lpval + L[cbind(seq_len(M), k)] - lse
      avail[cbind(seq_len(M), k)] <- FALSE
      mag <- p$vals[k]; x <- th * abs(cvec[i]) * mag
      s_al <- ifelse(runif(M) < plogis(2 * comp * x), 1, -1)     # +1 = aligned with c_i
      lsp <- lsp + plogis(2 * s_al * x, log.p = TRUE); lsm <- lsm + plogis(-2 * s_al * x, log.p = TRUE)
      V[, i] <- s_al * (if (cvec[i] >= 0) 1 else -1) * mag
    }
    nu_ <- length(upos)
    if (nu_ > 0 && !p$zero) {                            # uniform assignment of the rest
      keys <- matrix(runif(M * g), M, g); keys[!avail] <- Inf
      idx <- t(apply(keys, 1, order))[, seq_len(nu_), drop = FALSE]
      lpval <- lpval - lgamma(nu_ + 1)
      sg <- matrix(sample(c(-1, 1), M * nu_, replace = TRUE), M, nu_)
      V[, upos] <- sg * matrix(p$vals[idx], M, nu_)
      lsp <- lsp - nu_ * log(2); lsm <- lsm - nu_ * log(2)
    }
  }
  lp0 <- lp0 - K * log(2)
  mx <- pmax(lsp, lsm); lq <- lpval + log(0.5) + mx + log(exp(lsp - mx) + exp(lsm - mx))
  lw <- lp0 - lq
  lamV <- lrt_of(V, ellK, xiK, SSz, n, r)
  hit <- lamV >= lam_obs - 1e-9
  lp <- if (any(hit)) { v <- lw[hit]; mv <- max(v); mv + log(sum(exp(v - mv))) - log(M) } else -Inf
  c(lp = lp, lam = lam_obs, ess = if (any(hit)) { z <- lw[hit] - max(lw[hit]); sum(exp(z))^2 / sum(exp(2 * z)) } else 0)
}

ge_config <- function(sp, delta0, tag) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  has_p <- sp$nperp == 1L; nrest <- sp$nrest
  ellK <- c(ell, if (has_p) 0); xiK <- c(xi, if (has_p) xp)
  groups_by <- lapply(RATIOS, function(ra) make_groups(ell, ra))
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  tq_ <- qt(1 - AL / 2, df = nu - 1); cLR <- n * log1p(tq_^2 / (nu - 1))
  one <- function(lam) {
    beta <- lam / sqrt(S0)
    s <- rnorm(r) * sqrt(d0) + xi * beta
    zp <- if (has_p) rnorm(1) + xp * beta else numeric(0)
    zr <- if (nrest > 0) rnorm(nrest) else numeric(0)
    v <- c(s, zp); zero_vals <- c(zp, zr)
    out <- sapply(seq_along(RATIOS), function(k)
      ge_logp(v, zero_vals, ellK, xiK, groups_by[[k]], has_p, n, r)[c("lp", "ess")])
    lam_obs <- lrt_of(matrix(v, 1), ellK, xiK, sum(zero_vals^2), n, r)
    q0 <- { w <- 1 / (1 + delta0 * ell); A <- sum(s^2 * w) + sum(zero_vals^2); B <- sum(xi * s * w) + xp * sum(zp)
      Sx <- sum(xi^2 * w) + xp^2; (B / Sx) * sqrt(Sx * (nu - 1) / (A - B^2 / Sx)) }
    c(setNames(out["lp", ], paste0("lp_r", RATIOS)), setNames(out["ess", ], paste0("ess_r", RATIOS)),
      lrt = lam_obs, t_or = q0)
  }
  res <- list()
  H0 <- t(vapply(seq_len(N0), function(j) one(0), numeric(2 * length(RATIOS) + 2)))
  for (a0 in c(0.05, 0.01, if (TAIL) 1e-3)) {
    ct <- qt(1 - a0 / 2, nu - 1)
    res[[length(res) + 1]] <- data.frame(tag, lambda = 0, alpha = a0,
      method = c(paste0("ge_r", RATIOS), "lrt_ml", "oracle_t"),
      reject = c(colMeans(H0[, paste0("lp_r", RATIOS), drop = FALSE] <= log(a0)),
                 mean(H0[, "lrt"] > n * log1p(ct^2 / (nu - 1))), mean(abs(H0[, "t_or"]) > ct)),
      n_sim = N0)
  }
  if (N1 > 0) for (k in seq_along(AL)) {
    H1 <- t(vapply(seq_len(N1), function(j) one(tq_[k]), numeric(2 * length(RATIOS) + 2)))
    res[[length(res) + 1]] <- data.frame(tag, lambda = tq_[k], alpha = AL[k],
      method = c(paste0("ge_r", RATIOS), "lrt_ml", "oracle_t"),
      reject = c(colMeans(H1[, paste0("lp_r", RATIOS), drop = FALSE] <= log(AL[k])),
                 mean(H1[, "lrt"] > cLR[k]), mean(abs(H1[, "t_or"]) > tq_[k])),
      n_sim = N1)
  }
  out <- do.call(rbind, res)
  attr(out, "groups") <- vapply(groups_by, length, 1L)
  attr(out, "ess_h0") <- apply(H0[, paste0("ess_r", RATIOS), drop = FALSE], 2, median)
  out
}

## ---- configurations: 16 GEUVADIS + 2 AR(1) (one with a zero-eigenvalue group) --------
cfgs <- list()
for (i in seq_len(8)) for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "geu", i = i, h2 = h2)
cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = 1000L, nb = 5L, h2 = 0.05)   # m < n: zero group
cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = 200L, nb = 20L, h2 = 0.05)
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, j = match(g$lead_pos, pos))
}
if (REDUCED) cfgs <- Filter(function(z) z$kind == "geu" && z$i %in% c(1L, 2L, 3L, 6L), cfgs)
if (TAIL) cfgs <- Filter(function(z) z$kind == "geu" && z$i == 1L && z$h2 == 0.20, cfgs)
OUT <- paste0("calibration/results/group_exchange", if (REDUCED) "_reduced" else if (TAIL) "_tail" else "")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", n = z$n, blocks = z$nb,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(12100 + 1000 * z$n + z$nb + round(100 * z$h2))
  } else {
    id <- sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", n = nrow(L$Z), blocks = NA,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(12900 + z$i + round(100 * z$h2))
  }
  out <- ge_config(sp, z$h2 / (1 - z$h2), tag)
  if (!SMOKE) saveRDS(out, ck)
  cat(sprintf("[%7.1f s] %s done; groups %s; median IS ess under H0 %s\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), id,
              paste(attr(out, "groups"), collapse = "/"), paste(round(attr(out, "ess_h0")), collapse = "/")))
  if (SMOKE) out else ck
}
if (SMOKE) {
  out <- run_one(cfgs[[2]])                       # ERAP2 h2 = 0.20
  op <- options(width = 160); print(out[, c("lambda", "alpha", "method", "reject", "n_sim")], digits = 3)
  options(op)
} else {
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e)
    paste("ERROR", conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fx <- do.call(rbind, lapply(list.files(OUT, "rds$", full.names = TRUE), readRDS))
  write.csv(fx, paste0(OUT, ".csv"), row.names = FALSE)
  cat("configs done:", length(unique(fx$config)), "of", length(cfgs), "\n")
}
cat("done in", format(Sys.time() - t0), "\n")
