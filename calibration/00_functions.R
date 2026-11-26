# ==============================================================================
# 01_exact_size_pilot.R
#
# Question: does the sampling variability of the variance ratio delta matter for
# the conditional test of a candidate in a REGIONAL linear mixed model?
#
# Model (after projecting out the intercept):
#   y = x beta + u + e,  u ~ N(0, sigma^2 delta K),  e ~ N(0, sigma^2 I),
#   K = Z_{-j} Z_{-j}' / m_k  (regional kernel, candidate j excluded).
# Spectral reduction: s = Gamma' H' y ~ N(xi beta, sigma^2 diag(1 + delta l_i)).
#
# Plug-in Wald test: t(delta_hat) referred to Student t_{nu0-1}, delta_hat = REML
# with the candidate as fixed effect (per-candidate REML, GEMMA-type Wald test).
#
# Exact size by Rao-Blackwellised Monte Carlo (Kackar-Harville structure):
#   b = beta_hat(delta0) is independent of the error contrasts r, and
#   delta_hat, sigma_hat, d = beta_hat(delta_hat) - b are functions of r only, so
#   t(delta_hat) | r ~ N(m(r), v(r)),
#   m = d sqrt(S2(delta_hat)) / sigma_hat,  v = S2(delta_hat) / (S2(delta0) sigma_hat^2).
#   P(|t| > c) = E_r[ Phi((m - c)/sqrt v) + Phi((-c - m)/sqrt v) ]  (exact, any alpha).
#
# Closed-form indices at delta0 (a_i = l_i/(1+delta0 l_i), p = candidate loading):
#   nu_delta = nu0 Var_U(a) / (abar_p - abar_U)^2   (effective df of the delta part)
#   kappa    = Var(delta_hat) Var_p(a),  Var(delta_hat) ~ 2 / (nu0 Var_U(a))
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/01_exact_size_pilot.R [B=10000]
# Outputs: calibration/results/exact_size_{synthetic,geuvadis}.csv
# ==============================================================================

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1")
if (!exists("B_DRAW")) B_DRAW <- 10000L
CHUNK  <- 2500L
OUT    <- "calibration/results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ALPHA  <- c(5e-2, 1e-3, 1e-5, 5e-8)
GRID   <- c(0, 10^seq(-4, 4, length.out = 161))

## ---------------------------------------------------------------- spectral setup
spectral_setup <- function(Z, j = NULL, x = NULL) {
  n <- nrow(Z)
  if (is.null(x)) { x <- Z[, j]; Zk <- Z[, -j, drop = FALSE] } else Zk <- Z
  x  <- x - mean(x)
  mk <- ncol(Zk)
  ev <- eigen(tcrossprod(Zk) / mk, symmetric = TRUE)
  keep <- ev$values > 1e-8 * max(ev$values)
  U  <- ev$vectors[, keep, drop = FALSE]
  ell <- ev$values[keep]
  xi <- drop(crossprod(U, x))
  xp <- sqrt(max(sum(x^2) - sum(xi^2), 0))
  nu0 <- n - 1L; r <- length(ell)
  nperp <- as.integer(xp > 1e-6 * sqrt(sum(x^2)))
  if (!nperp) xp <- 0
  list(ell = ell, xi = xi, xp = xp, nu0 = nu0, r = r, nperp = nperp,
       nrest = nu0 - r - nperp, U = U)
}

## ------------------------------------------------------------- closed-form index
indices <- function(sp, delta0) {
  w <- 1 / (1 + delta0 * sp$ell); a <- sp$ell * w
  S2 <- sum(sp$xi^2 * w) + sp$xp^2
  aU <- sum(a) / sp$nu0; vU <- sum(a^2) / sp$nu0 - aU^2
  p  <- sp$xi^2 * w / S2
  ap <- sum(p * a); vp <- sum(p * a^2) - ap^2
  nu_d <- sp$nu0 * vU / (ap - aU)^2
  c(nu_delta = nu_d, nu_total = 1 / (1 / (sp$nu0 - 1) + 1 / nu_d),
    kappa = 2 * vp / (sp$nu0 * vU), sd_delta_asy = sqrt(2 / (sp$nu0 * vU)),
    abar_p = ap, abar_U = aU, lbar_p = sum(p * sp$ell), frac_kernel = sum(p))
}

## --------------------------------------------------------------- exact size by RB-MC
exact_size <- function(sp, delta0, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu0 <- sp$nu0; r <- sp$r
  cc  <- qt(1 - ALPHA / 2, df = nu0 - 1)
  Wg  <- 1 / (1 + outer(ell, GRID))                          # r x G
  ldg <- colSums(log(1 + outer(ell, GRID)))                  # G
  S2g <- colSums(xi^2 * Wg) + xp^2                           # G
  d0  <- 1 + delta0 * ell
  acc <- list(); nch <- ceiling(B / CHUNK)
  for (ch in seq_len(nch)) {
    b  <- CHUNK
    s  <- matrix(rnorm(r * b), r, b) * sqrt(d0)              # kernel directions
    sp_ <- if (sp$nperp) rnorm(b) else numeric(b)            # x_perp direction
    ssr <- if (sp$nrest > 0) rchisq(b, sp$nrest) else numeric(b)
    s2 <- s * s
    ## REML with x as fixed effect, grid + parabolic refinement on log(delta)
    Ag  <- crossprod(s2, Wg) + (sp_^2 + ssr)                 # b x G
    Bg  <- crossprod(xi * s, Wg) + xp * sp_                  # b x G
    RSSg <- Ag - sweep(Bg^2, 2, S2g, "/")
    obj <- -0.5 * (sweep((nu0 - 1) * log(RSSg), 2, ldg + log(S2g), "+"))
    ## null-model REML (no fixed effect) for the score-type variant
    obj0 <- -0.5 * (sweep(nu0 * log(Ag), 2, ldg, "+"))
    pick <- function(O) {
      g <- max.col(O, ties.method = "first"); dh <- GRID[g]
      int <- g > 2 & g < length(GRID)                        # interior on log grid
      if (any(int)) {
        ii <- which(int); lg <- log(GRID)
        x1 <- lg[g[ii] - 1]; x2 <- lg[g[ii]]; x3 <- lg[g[ii] + 1]
        f1 <- O[cbind(ii, g[ii] - 1)]; f2 <- O[cbind(ii, g[ii])]; f3 <- O[cbind(ii, g[ii] + 1)]
        den <- (x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1)
        num <- (x2 - x1)^2 * (f2 - f3) - (x2 - x3)^2 * (f2 - f1)
        xs <- ifelse(abs(den) > 1e-14, x2 - 0.5 * num / den, x2)
        xs <- pmin(pmax(xs, x1), x3)
        dh[ii] <- exp(xs)
      }
      dh
    }
    dh  <- pick(obj); dh0 <- pick(obj0)
    ## quantities at delta_hat (per replicate)
    at_delta <- function(dl) {
      D  <- 1 + outer(ell, dl); iD <- 1 / D
      A  <- colSums(s2 * iD) + sp_^2 + ssr
      Bx <- colSums(xi * s * iD) + xp * sp_
      S2 <- colSums(xi^2 * iD) + xp^2
      list(A = A, Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
    }
    Q  <- at_delta(dh)
    bh <- Q$Bx / Q$S2; sig2 <- Q$RSS / (nu0 - 1)
    tt <- bh * sqrt(Q$S2 / sig2)
    ## truth
    Q0 <- at_delta(rep(delta0, b))
    b0 <- Q0$Bx / Q0$S2
    dd <- bh - b0
    m  <- dd * sqrt(Q$S2 / sig2); v <- Q$S2 / (Q0$S2 * sig2)
    rb <- sapply(cc, function(c0) pnorm((m - c0) / sqrt(v)) + pnorm((-c0 - m) / sqrt(v)))
    ## oracle check (delta known): must reproduce alpha
    sig0 <- Q0$RSS / (nu0 - 1)
    orc <- sapply(cc, function(c0) 2 * pnorm(-c0 * sqrt(sig0)))
    ## score-type statistic with null-model delta_hat (plain MC, moderate alpha only)
    Qn <- at_delta(dh0)
    tsc <- (Qn$Bx / sqrt(Qn$S2)) / sqrt(Qn$A / nu0)
    acc[[ch]] <- list(rb = rb, orc = orc,
                      emp = sapply(cc[1:2], function(c0) as.numeric(abs(tt) > c0)),
                      emp0 = sapply(cc[1:2], function(c0) as.numeric(abs(tsc) > c0)),
                      dh = dh, dh0 = dh0, v = v)
  }
  rb  <- do.call(rbind, lapply(acc, `[[`, "rb"))
  orc <- do.call(rbind, lapply(acc, `[[`, "orc"))
  emp <- do.call(rbind, lapply(acc, `[[`, "emp"))
  emp0 <- do.call(rbind, lapply(acc, `[[`, "emp0"))
  dh  <- unlist(lapply(acc, `[[`, "dh")); dh0 <- unlist(lapply(acc, `[[`, "dh0"))
  size <- colMeans(rb); se <- apply(rb, 2, sd) / sqrt(nrow(rb))
  out <- c(p_dhat0 = mean(dh < 1e-4), med_dhat_rel = median(dh) / delta0,
           q10_dhat_rel = unname(quantile(dh, .1)) / delta0,
           q90_dhat_rel = unname(quantile(dh, .9)) / delta0,
           p_dhat0_null = mean(dh0 < 1e-4))
  for (k in seq_along(ALPHA)) {
    nm <- format(ALPHA[k], scientific = TRUE)
    out[paste0("ratio_", nm)]  <- size[k] / ALPHA[k]
    out[paste0("rse_", nm)]    <- se[k] / size[k]
    out[paste0("oracle_", nm)] <- mean(orc[, k]) / ALPHA[k]
  }
  out["emp_ratio_5e-02"]  <- mean(emp[, 1]) / ALPHA[1]
  out["emp_ratio_1e-03"]  <- mean(emp[, 2]) / ALPHA[2]
  out["score_null_ratio_5e-02"] <- mean(emp0[, 1]) / ALPHA[1]
  out["score_null_ratio_1e-03"] <- mean(emp0[, 2]) / ALPHA[2]
  out
}

pred_ratio <- function(nu_tot, nu0) {
  cc <- qt(1 - ALPHA / 2, df = nu0 - 1)
  setNames(2 * pt(-cc, df = nu_tot) / ALPHA,
           paste0("satt_ratio_", format(ALPHA, scientific = TRUE)))
}

run_config <- function(sp, delta0, tag) {
  ix <- indices(sp, delta0)
  ex <- exact_size(sp, delta0)
  pr <- pred_ratio(ix[["nu_total"]], sp$nu0)
  c(tag, list(r = sp$r, nu0 = sp$nu0, delta0 = delta0), as.list(ix), as.list(ex),
    as.list(pr))
}

## ----------------------------------------------------------------- synthetic designs
gen_blocks <- function(n, nb, bs, rho) {
  Z <- matrix(0, n, nb * bs)
  for (b in seq_len(nb)) {
    E <- matrix(rnorm(n * bs), n, bs); X <- E
    for (k in 2:bs) X[, k] <- rho * X[, k - 1] + sqrt(1 - rho^2) * E[, k]
    Z[, (b - 1) * bs + seq_len(bs)] <- X
  }
  scale(Z)
}

