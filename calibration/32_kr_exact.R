# ==============================================================================
# 32_kr_exact.R
#
# Exact size and power of the standard small-sample corrections for the test of one
# fixed effect when the variance ratio is estimated: Kenward-Roger (1997) as in
# pbkrtest and MM4LMM, and Satterthwaite (1946) degrees of freedom as in lmerTest
# (here with the expected REML information). Both are functions of the REML fit with
# the candidate, hence of the error contrasts r only, and the Wald statistic is
# conditionally normal given r, so their exact laws are
#   P(reject) = E_r[ Phi((mu - c(r))/sqrt v) + Phi((-c(r) - mu)/sqrt v) ].
# Linear covariance parametrization theta = (V_g, V_e), V = V_g L + V_e I in the reduced
# coordinates (as pbkrtest / MM4LMM). With P = xi'V^{-1}xi, Phi = 1/P, W = inverse
# expected REML information, P_i = d P / d theta_i, Q_ij = xi'V^{-1}G_iV^{-1}G_jV^{-1}xi:
#   Lambda = Phi { sum_ij W_ij (Q_ij - P_i Phi P_j) } Phi,  Phi_A = Phi + 2 Lambda,
#   t_KR = beta_hat / sqrt(Phi_A) ~ t_m,  m = 2 / A,  A = sum_ij W_ij (P_i/P)(P_j/P) (Phi version)
#   or A_A = (Phi_A/Phi)^2 A (pbkrtest, Phi_A in Theta and U); lambda = 1 for one coefficient.
#   Satterthwaite: t(delta_hat) ~ t_{2/A}.
# Draws, seeds and configurations are those of 06_treatments_pilot.R and
# 22_edgeworth_check.R (8 AR(1) designs, 16 GEUVADIS windows, B = 1000), so the
# uncorrected row wald_t reproduces their wald_h1 / wald_t rows exactly.
# Part B validates the closed forms against MM4LMM::AnovaTest(Type = "KR") on GEUVADIS.
#
# Usage (from the repository root): nice -n 15 Rscript calibration/32_kr_exact.R [cores=2] [--smoke] [--validate]
# Output: results/kr_exact.csv (long format), results/kr_validation.csv
# ==============================================================================
args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 2L
SMOKE  <- "--smoke" %in% args; VALIDATE <- "--validate" %in% args
B_DRAW <- if (SMOKE) 200L else 1000L
source("calibration/00_functions.R")
suppressPackageStartupMessages(library(parallel))
GRID2 <- c(0, 10^seq(-4, 4, length.out = 81))
AL    <- c(1e-3, 1e-5, 5e-8)

pick_grid <- function(O, grid) {                          # as in 06 and 22
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

## Kenward-Roger and Satterthwaite pieces at (V_g, V_e) = (delta, 1); all ratios are scale-free
kr_factors <- function(ell, xi2, xp2, nout, dl) {
  w <- 1 / (1 + outer(ell, dl)); w2 <- w * w; w3 <- w2 * w
  P   <- colSums(xi2 * w) + xp2
  Pg1 <- colSums(xi2 * ell * w2); Pe1 <- colSums(xi2 * w2) + xp2          # -dP/dV_g, -dP/dV_e
  Ugg <- colSums(xi2 * ell^2 * w3); Uge <- colSums(xi2 * ell * w3); Uee <- colSums(xi2 * w3) + xp2
  Tgg <- colSums(ell^2 * w2); Tge <- colSums(ell * w2); Tee <- colSums(w2) + nout
  Igg <- 0.5 * (Tgg - 2 * Ugg / P + Pg1^2 / P^2)                        # REML information, candidate fixed
  Ige <- 0.5 * (Tge - 2 * Uge / P + Pg1 * Pe1 / P^2)
  Iee <- 0.5 * (Tee - 2 * Uee / P + Pe1^2 / P^2)
  det <- Igg * Iee - Ige^2
  hgg <- Iee / det; hge <- -Ige / det; hee <- Igg / det
  A <- (hgg * Pg1^2 + 2 * hge * Pg1 * Pe1 + hee * Pe1^2) / P^2            # g'Wg, g = grad log Phi
  LoP <- (hgg * (Ugg - Pg1^2 / P) + 2 * hge * (Uge - Pg1 * Pe1 / P) + hee * (Uee - Pe1^2 / P)) / P
  ratio <- 1 + 2 * LoP                                                    # Phi_A / Phi
  list(ratio = ratio, df_satt = 2 / A, df_krA = 2 / (ratio^2 * A), A = A)
}

kr_config <- function(sp, delta0, tag, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  S2g <- colSums(xi^2 * W2) + xp^2
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)              # same draws as 06 and 22
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  Rk2 <- Rk * Rk
  A_r <- crossprod(Rk2, W2) + (Rp^2 + ssr); B_r <- crossprod(xi * Rk, W2) + xp * Rp
  RSS_r <- A_r - sweep(B_r^2, 2, S2g, "/")
  l1  <- -0.5 * sweep((nu - 1) * log(RSS_r), 2, ld2 + log(S2g), "+")
  dh1 <- pick_grid(l1, GRID2)                                # REML with the candidate
  at_r <- function(dl) {
    iD <- 1 / (1 + outer(ell, dl))
    A <- colSums(Rk2 * iD) + Rp^2 + ssr; Bx <- colSums(xi * Rk * iD) + xp * Rp
    S2 <- colSums(xi^2 * iD) + xp^2
    list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
  }
  Wd <- at_r(dh1); sig_h <- sqrt(Wd$RSS / (nu - 1))
  wa <- sqrt(Wd$S2 / S0) / sig_h; wg <- (Wd$Bx / Wd$S2) * sqrt(Wd$S2) / sig_h
  lin <- function(a, g, c0) list(up = (c0 - g) / a, dn = (-c0 - g) / a)
  rej <- function(rt, lam) mean(pnorm(lam - rt$up) + pnorm(rt$dn - lam))
  kf <- kr_factors(ell, xi^2, xp^2, sp$nperp + sp$nrest, dh1)
  tq_ <- qt(1 - AL / 2, df = nu - 1)
  res <- list()
  add <- function(method, k, rt) res[[length(res) + 1]] <<- data.frame(tag, method = method, alpha = AL[k],
      size_ratio = rej(rt, 0) / AL[k], power = rej(rt, tq_[k]), lambda = tq_[k])
  for (k in seq_along(AL)) {
    add("wald_t", k, lin(wa, wg, tq_[k]))
    add("satterthwaite", k, lin(wa, wg, qt(1 - AL[k] / 2, df = kf$df_satt)))
    add("kr", k, lin(wa, wg, qt(1 - AL[k] / 2, df = kf$df_satt) * sqrt(kf$ratio)))
    add("kr_phiA", k, lin(wa, wg, qt(1 - AL[k] / 2, df = kf$df_krA) * sqrt(kf$ratio)))
    add("vcov_adj_only", k, lin(wa, wg, tq_[k] * sqrt(kf$ratio)))
  }
  out <- do.call(rbind, res)
  attr(out, "factors") <- data.frame(tag, ratio_med = median(kf$ratio), ratio_q90 = unname(quantile(kf$ratio, 0.9)),
    df_satt_med = median(kf$df_satt), df_satt_q10 = unname(quantile(kf$df_satt, 0.1)),
    df_krA_med = median(kf$df_krA), p_boundary = mean(dh1 < 1e-4))
  out
}

## ---- configurations, seeds and data as in 06_treatments_pilot.R and 22_edgeworth_check.R ----
cfgs <- list()
for (n in c(200L, 1000L)) for (nb in c(5L, 20L)) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "syn", n = n, nb = nb, h2 = h2)
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
for (i in seq_len(nrow(GENES))) for (h2 in c(0.05, 0.20))
  cfgs[[length(cfgs) + 1]] <- list(kind = "geu", i = i, h2 = h2)
if (SMOKE) cfgs <- cfgs[c(which(vapply(cfgs, function(z) z$kind == "geu" && z$i == 1 && z$h2 == 0.05, TRUE)),
                          which(vapply(cfgs, function(z) z$kind == "syn" && z$n == 1000 && z$nb == 20 && z$h2 == 0.20, TRUE)))]
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

## ---- Part B: validation against MM4LMM (Kenward-Roger) on the GEUVADIS leads --------------
if (VALIDATE) {
  suppressPackageStartupMessages(library(MM4LMM))
  rows <- list()
  for (i in c(1L, 2L, 6L)) {
    L <- load_geu(i); Z <- L$Z; j <- L$j; n <- nrow(Z)
    set.seed(91 + i)
    Kj <- tcrossprod(Z[, -j]) / (ncol(Z) - 1)
    for (h2 in c(0.05, 0.20)) {
      d0 <- h2 / (1 - h2); ev <- eigen(Kj, symmetric = TRUE)
      y <- as.numeric(ev$vectors %*% (rnorm(n) * sqrt(1 + d0 * pmax(ev$values, 0))))   # null phenotype
      x <- Z[, j]
      R <- MMEst(Y = y, Cofactor = matrix(1, n, 1), X = matrix(x, n, 1, dimnames = list(NULL, "x")),
                 VarList = list(K = Kj, Error = diag(n)), CritVar = 1e-10, CritLogLik = 1e-10, MaxIter = 5000,
                 Verbose = FALSE)
      TC <- matrix(c(0, 1), nrow = 1)                     # explicit contrast (MM4LMM 3.0.3 needs it for KR)
      kr <- AnovaTest(R, TestedCombination = TC, Type = "KR", Cofactor = matrix(1, n, 1),
                      X = matrix(x, n, 1, dimnames = list(NULL, "x")), VarList = list(K = Kj, Error = diag(n)))
      w3 <- AnovaTest(R, TestedCombination = TC)
      sig <- R[[1]]$Sigma2; dl <- sig[1] / sig[2]
      ## same quantities from the reduced model at the MM4LMM REML estimate
      sp <- spectral_setup(Z, j = j)
      kf <- kr_factors(sp$ell, sp$xi^2, sp$xp^2, sp$nperp + sp$nrest, dl)
      ## Wald t at that delta from the reduced coordinates
      s <- as.numeric(crossprod(sp$U, y - mean(y)))
      w <- 1 / (1 + dl * sp$ell); Sd <- sum(sp$xi^2 * w) + sp$xp^2
      xperp <- (x - mean(x)) - as.numeric(sp$U %*% sp$xi); yc <- y - mean(y)
      yperp <- yc - as.numeric(sp$U %*% s)
      Bx <- sum(sp$xi * s * w) + sum(xperp * yperp)
      A <- sum(s^2 * w) + sum(yperp^2); RSS <- A - Bx^2 / Sd
      t2 <- (Bx / Sd)^2 * Sd * (n - 2) / RSS
      kr_v <- kr[[1]][, 1]; w3_v <- unlist(w3[[1]])
      rows[[length(rows) + 1]] <- data.frame(gene = GENES$symbol[i], h2 = h2, delta_mm4lmm = dl,
        F_kr_mm4lmm = kr_v[["Stat"]], df_kr_mm4lmm = kr_v[["ddf"]], p_kr_mm4lmm = kr_v[["pval"]],
        F_kr_ours = t2 / kf$ratio, df_satt_ours = kf$df_satt, df_krA_ours = kf$df_krA, ratio_ours = kf$ratio,
        p_kr_ours = pf(t2 / kf$ratio, 1, kf$df_satt, lower.tail = FALSE),
        wald_mm4lmm = w3_v[1], t2_ours = t2, w3 = paste(names(w3_v), signif(w3_v, 6), collapse = " "))
    }
  }
  V <- do.call(rbind, rows); write.csv(V, "calibration/results/kr_validation.csv", row.names = FALSE)
  print(V, digits = 6)
  quit(save = "no")
}

OUT <- "calibration/results/kr_exact"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", n = z$n, blocks = z$nb, h2 = z$h2,
                      nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2))
  } else {
    id <- sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck) && !SMOKE) return(ck)
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", n = nrow(L$Z), blocks = NA, h2 = z$h2,
                      nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(7000 + z$i + round(100 * z$h2))
  }
  out <- kr_config(sp, z$h2 / (1 - z$h2), tag)
  if (SMOKE) return(out)
  saveRDS(out, ck); ck
}
if (SMOKE) {
  for (z in cfgs) { t0 <- Sys.time(); out <- run_one(z)
    cat(sprintf("\n%s  nu_delta = %.1f  (%.0f s)\n", out$config[1], out$nu_delta[1],
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    print(reshape(out[, c("method", "alpha", "size_ratio")], idvar = "method", timevar = "alpha",
                  direction = "wide"), row.names = FALSE, digits = 3)
    print(attr(out, "factors"), row.names = FALSE, digits = 3) }
} else {
  st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e) paste("ERROR", conditionMessage(e))),
                 mc.cores = NCORES, mc.preschedule = FALSE)
  print(unlist(st)[grepl("ERROR", unlist(st))])
  fs <- list.files(OUT, "rds$", full.names = TRUE); L <- lapply(fs, readRDS)
  write.csv(do.call(rbind, L), "calibration/results/kr_exact.csv", row.names = FALSE)
  write.csv(do.call(rbind, lapply(L, attr, "factors")), "calibration/results/kr_factors.csv", row.names = FALSE)
  cat("configs done:", length(fs), "of", length(cfgs), "\n")
}
