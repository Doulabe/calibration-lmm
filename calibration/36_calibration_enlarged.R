# ==============================================================================
# 36_calibration_enlarged.R
#
# Enlarged exact-calibration study for the article: many candidates
# per GEUVADIS window (not only the published lead), chosen to span the range of
# the index nu_j, plus AR(1) designs with few and many blocks and candidates
# unrelated to the region. For each configuration, exact size relative to alpha and
# exact power at lambda = t_{nu-1, alpha} of:
#   oracle_t        Student test with delta known
#   wald_t          Wald, delta_hat with the candidate (REML), Student nu-1 (GEMMA, MM4LMM default)
#   score_t         score, null delta_hat_0 (REML), Student nu-1 (EMMAX, P3D, GCTA, MLMM step)
#   lrt_chi2        likelihood ratio, ML under each model, chi2_1 (GEMMA LRT)
#   kr              Kenward-Roger (pbkrtest, MM4LMM, SAS)
#   satterthwaite   Satterthwaite degrees of freedom, expected information (lmerTest-type)
#   wald_edgeworth, lrt_edgeworth   Edgeworth thresholds (GEMMA -crt, Zhou & Stephens 2014)
# The computations are those of 22_edgeworth_check.R and 32_kr_exact.R (same formulas,
# same grid, B = 1000 draws of the error contrasts per configuration).
# Candidates of the GEUVADIS windows (n = 358): the published lead and nine variants at
# the quantiles 0.02, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.98 of log nu over the
# variants of the window (nu from results/nu_delta_per_snp_geuvadis.csv), for
# h2 in {0.05, 0.20}; the kernel is the window without the candidate.
# AR(1) designs (blocks of 50 markers, rho = 0.95): n in {200, 1000} x {5, 20} blocks
# (the 8 configurations of 06/22/32), n = 1000 with 100 blocks, and a candidate
# unrelated to the region for n in {200, 1000} with 20 blocks.
#
# Usage (from the repository root): nice -n 10 Rscript calibration/36_calibration_enlarged.R [cores=3] [--smoke]
# Output: results/calib_enlarged/<config>.rds, results/calib_enlarged.csv
# ==============================================================================
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
SMOKE  <- "--smoke" %in% args
B_DRAW <- if (SMOKE) 200L else 1000L
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
AL    <- c(1e-3, 1e-5, 5e-8)
NBIS  <- 28L

pick_grid <- function(O, grid) {                          # as in 06, 22, 32
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
## Rothenberg / Zhou-Stephens factors (22_edgeworth_check.R)
crt_factors <- function(sp, dl, reml = TRUE) {
  ell <- sp$ell; xi2 <- sp$xi^2; xp2 <- sp$xp^2; nout <- sp$nperp + sp$nrest
  w <- 1 / (1 + outer(ell, dl)); w2 <- w * w; w3 <- w2 * w
  P   <- colSums(xi2 * w) + xp2
  Pg1 <- colSums(xi2 * ell * w2); Pe1 <- colSums(xi2 * w2) + xp2
  Ugg <- colSums(xi2 * ell^2 * w3); Uge <- colSums(xi2 * ell * w3); Uee <- colSums(xi2 * w3) + xp2
  Tgg <- colSums(ell^2 * w2); Tge <- colSums(ell * w2); Tee <- colSums(w2) + nout
  if (reml) {
    Igg <- 0.5 * (Tgg - 2 * Ugg / P + Pg1^2 / P^2)
    Ige <- 0.5 * (Tge - 2 * Uge / P + Pg1 * Pe1 / P^2)
    Iee <- 0.5 * (Tee - 2 * Uee / P + Pe1^2 / P^2)
  } else { Igg <- 0.5 * Tgg; Ige <- 0.5 * Tge; Iee <- 0.5 * Tee }
  det <- Igg * Iee - Ige^2
  hgg <- Iee / det; hge <- -Ige / det; hee <- Igg / det
  cg <- -Pg1 / P; ce <- -Pe1 / P
  Lgg <- 2 * Ugg / P - cg^2; Lge <- 2 * Uge / P - cg * ce; Lee <- 2 * Uee / P - ce^2
  V <- hgg * cg^2 + 2 * hge * cg * ce + hee * ce^2
  H <- hgg * Lgg + 2 * hge * Lge + hee * Lee
  ## Kenward-Roger pieces (32_kr_exact.R): ratio Phi_A / Phi and Satterthwaite df 2 / V
  LoP <- (hgg * (Ugg - Pg1^2 / P) + 2 * hge * (Uge - Pg1 * Pe1 / P) + hee * (Uee - Pe1^2 / P)) / P
  list(a = H - V / 2, b = H - V, c = 1.5 * V, ratio = 1 + 2 * LoP, df = 2 / V)
}

calib_config <- function(sp, delta0, tag, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1  <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+")
  dh1 <- pick_grid(l1, GRID2)                                # REML with the candidate
  L1ml <- -sweep(n * log(RSS_r), 2, ld2, "+")
  ml1 <- apply(L1ml, 1, max); dh1ml <- pick_grid(L1ml, GRID2) # ML with the candidate
  at_r <- function(dl) {
    iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Rk2 * iD) + Rp^2 + ssr; Bx <- colSums(xi * Rk * iD) + xp * Rp
    S2 <- colSums(xi^2 * iD) + xp^2
    list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
  }
  Wd <- at_r(dh1); sig_h <- sqrt(Wd$RSS / (nu - 1))
  wa <- sqrt(Wd$S2 / S0) / sig_h; wg <- (Wd$Bx / Wd$S2) * sqrt(Wd$S2) / sig_h
  O0 <- at_r(rep(delta0, B)); sig0 <- sqrt(O0$RSS / (nu - 1))
  ev <- function(u, what) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    A2 <- crossprod(sk2, W2) + (s2p^2 + ssr); out <- list()
    if ("lrt" %in% what) out$lrt <- ml1 - apply(-sweep(n * log(A2), 2, ld2, "+"), 1, max)
    if ("Q" %in% what) {
      dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
      iD <- 1 / (1 + outer(ell, dh0))
      Ah <- colSums(sk2 * iD) + s2p^2 + ssr; Bh <- colSums(xi * sk * iD) + xp * s2p
      Sh <- colSums(xi^2 * iD) + xp^2
      out$Q <- Bh / sqrt(Sh * Ah / nu)
    }
    out
  }
  lo <- rep(-30, B); hi <- rep(30, B)                        # split point, Q(u) = 0
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < 0
    lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  rootsU <- function(key, c0) {
    f <- function(u) ev(u, key)[[key]]
    always <- f(split) >= c0
    lo <- split; hi <- rep(120, B); never_r <- f(hi) < c0
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; up <- f(mid) < c0
      lo[up] <- mid[up]; hi[!up] <- mid[!up] }
    ur <- (lo + hi) / 2; ur[never_r] <- Inf
    lo <- rep(-120, B); hi <- split; never_l <- f(lo) < c0
    for (it in seq_len(NBIS)) { mid <- (lo + hi) / 2; dn <- f(mid) < c0
      hi[dn] <- mid[dn]; lo[!dn] <- mid[!dn] }
    ul <- (lo + hi) / 2; ul[never_l] <- -Inf
    ur[always] <- split[always]; ul[always] <- split[always]
    list(up = ur, dn = ul)
  }
  rootsQ <- function(c0) {
    one <- function(level) { lo <- rep(-120, B); hi <- rep(120, B)
      for (it in seq_len(NBIS + 6)) { mid <- (lo + hi) / 2; up <- ev(mid, "Q")$Q < level
        lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      (lo + hi) / 2 }
    list(up = one(c0), dn = one(-c0))
  }
  lin <- function(a, g, c0) list(up = (c0 - g) / a, dn = (-c0 - g) / a)
  rej <- function(rt, lam) mean(pnorm(lam - rt$up) + pnorm(rt$dn - lam))
  fW <- crt_factors(sp, dh1, reml = TRUE); fL <- crt_factors(sp, dh1ml, reml = FALSE)
  tq_ <- qt(1 - AL / 2, df = nu - 1); cq <- qchisq(1 - AL, 1)
  res <- list()
  add <- function(method, k, rt) res[[length(res) + 1]] <<- data.frame(tag, method = method, alpha = AL[k],
      size_ratio = rej(rt, 0) / AL[k], power = rej(rt, tq_[k]), lambda = tq_[k])
  for (k in seq_along(AL)) {
    add("oracle_t", k, list(up = tq_[k] * sig0, dn = -tq_[k] * sig0))
    add("wald_t", k, lin(wa, wg, tq_[k]))
    add("kr", k, lin(wa, wg, qt(1 - AL[k] / 2, df = fW$df) * sqrt(fW$ratio)))
    add("satterthwaite", k, lin(wa, wg, qt(1 - AL[k] / 2, df = fW$df)))
    add("wald_edgeworth", k, lin(wa, wg, sqrt(pmax(cq[k] * (1 + (fW$a + fW$b) / 2 + fW$c * cq[k] / 6), 0))))
    add("lrt_chi2", k, rootsU("lrt", cq[k]))
    add("lrt_edgeworth", k, rootsU("lrt", cq[k] * (1 + fL$a / 2)))
    add("score_t", k, rootsQ(tq_[k]))
  }
  out <- do.call(rbind, res)
  attr(out, "extra") <- data.frame(tag, p_boundary = mean(dh1 < 1e-4), kr_df_med = median(fW$df),
                                   kr_df_q10 = unname(quantile(fW$df, 0.1)), kr_ratio_med = median(fW$ratio))
  out
}

## ---- configurations ---------------------------------------------------------------------------
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L, 213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)
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
  list(Z = Z, pos = pos)
}
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
    role = if (k == 1) "lead" else sprintf("q%02d", round(100 * QS[k - 1])))
}
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2, unlinked = FALSE)
for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = 1000L, nb = 100L, h2 = h2, unlinked = FALSE)
for (n in c(200L, 1000L)) for (h2 in c(0.05, 0.20)) cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = 20L, h2 = h2, unlinked = TRUE)
if (SMOKE) cfgs <- cfgs[c(1, 2, length(cfgs) - 6)]

OUT <- "calibration/results/calib_enlarged"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GEU <- list()
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "")
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L)
  } else {
    g <- GENES[z$i, ]; id <- sprintf("geu_%s_%d_h%.2f", g$symbol, z$pos, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973))
  }
  out <- calib_config(sp, z$h2 / (1 - z$h2), tag)
  if (SMOKE) return(out)
  saveRDS(out, ck); ck
}
if (SMOKE) {
  for (z in cfgs) { t0 <- Sys.time(); out <- run_one(z)
    cat(sprintf("\n%s  nu_delta = %.1f  (%.0f s)\n", out$config[1], out$nu_delta[1],
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(reshape(out[, c("method", "alpha", "size_ratio")], idvar = "method", timevar = "alpha", direction = "wide"),
          row.names = FALSE, digits = 3) }
} else {
  cat("configurations:", length(cfgs), "\n")
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) paste("ERROR", conditionMessage(e))),
                 mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT, "rds$", full.names = TRUE); L <- lapply(fs, readRDS)
  write.csv(do.call(rbind, L), "calibration/results/calib_enlarged.csv", row.names = FALSE)
  write.csv(do.call(rbind, lapply(L, attr, "extra")), "calibration/results/calib_enlarged_extra.csv", row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
