# ==============================================================================
# 06b_nonmonotone_fix.R
#
# 06_treatments_pilot.R located rejection regions assuming each statistic is
# monotone in |u - split| on each side. A diagnostic on ERAP2 showed that the
# marginal Bayes factor (and the plug-in BF at delta_hat0) is NOT monotone: it
# rises, peaks near u ~ 20-40 and falls back towards 1 for very large effects,
# because the null model absorbs the effect into an ever larger delta and the
# slab variance tau2 is fixed (information paradox of fixed g-priors, Liang et
# al. 2008). A bounded rejection interval was then reported as empty.
#
# This script recomputes, on the SAME error-contrast draws (same seeds, same
# draw order as 06), the four statistics whose shape is not guaranteed:
#   score_h0 (|Q|), signflip (-log p), bf_h0, mbf.
# On each side of the split point the function is maximised by golden-section
# search on g in [0, 25] (larger |u| has negligible probability for lambda <= 6);
# the rejection interval is [rising crossing, falling crossing], each found by
# bisection, the falling crossing being set to +Inf when the statistic is still
# above the threshold at g = 30. P(reject | r) sums the two side intervals.
# Also reports the share of draws whose side maximum is interior (non-monotone).
#
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/06b_nonmonotone_fix.R [cores=3]
# Output: calibration/results/treatments_pilot_fix/<config>.rds and
#         calibration/results/treatments_pilot_final.csv (06 rows for the
#         other treatments + corrected rows for these four)
# ==============================================================================

args   <- commandArgs(TRUE)
NCORES <- if (length(args) >= 1 && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
src <- readLines("calibration/06_treatments_pilot.R")
cut <- grep("^## ---- configurations", src)
eval(parse(text = src[1:(cut - 1)]))           # definitions of 06 (B_DRAW, grids, helpers)
B_DRAW <- 1000L
GMAX <- 25; GFALL <- 30; NGOLD <- 32L; NB2 <- 26L

fix_config <- function(sp, delta0, tag, B = B_DRAW) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu <- sp$nu0; r <- sp$r; n <- nu + 1L
  W2 <- 1 / (1 + outer(ell, GRID2)); ld2 <- colSums(log(1 + outer(ell, GRID2)))
  WQ <- 1 / (1 + outer(ell, QG)); ldQ <- colSums(log(1 + outer(ell, QG)))
  SQg <- colSums(xi^2 * WQ) + xp^2; vQ <- TAU2 * SQg
  tq <- log(QG); dq <- tq[2] - tq[1]
  lwq <- log(c(dq / 2, rep(dq, length(tq) - 2), dq / 2)) + 0.5 * tq - log(pi) - log1p(QG)
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  ## identical draws to 06 (same order: s, spp, ssr)
  s   <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B)
  ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b   <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk  <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u, key) {
    bv <- u / sqrt(S0)
    sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv; sk2 <- sk * sk
    if (key == "mbf") {
      AQ <- crossprod(sk2, WQ) + (s2p^2 + ssr); BQ <- crossprod(xi * sk, WQ) + xp * s2p
      l0 <- -0.5 * sweep(nu * log(AQ), 2, ldQ, "+")
      lbf <- sweep(-(nu / 2) * log1p(-sweep(BQ^2 / AQ, 2, TAU2 / (1 + vQ), "*")), 2,
                   0.5 * log1p(vQ), "-")
      lw <- sweep(l0, 2, lwq, "+")
      return(rowLSE(lw + lbf) - rowLSE(lw))
    }
    A2 <- crossprod(sk2, W2) + (s2p^2 + ssr)
    dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
    iD <- 1 / (1 + outer(ell, dh0))
    Ah <- colSums(sk2 * iD) + s2p^2 + ssr
    Bh <- colSums(xi * sk * iD) + xp * s2p
    Sh <- colSums(xi^2 * iD) + xp^2
    Q <- Bh / sqrt(Sh * Ah / nu)
    if (key == "absQ") return(abs(Q))
    if (key == "bf0") { v <- TAU2 * Sh
      return(-0.5 * log1p(v) - (nu / 2) * log1p(-v * Q^2 / (nu * (1 + v)))) }
    if (key == "sf") return(-log_signflip_p(rbind(abs(xi * iD * sk), abs(xp * s2p)), abs(Bh)))
    stop("unknown key")
  }
  ## split point: signed score numerator at delta_hat0 equals 0 (as in 06)
  numQ <- function(u) {
    bv <- u / sqrt(S0); sk <- Rk + outer(xi, bv); s2p <- Rp + xp * bv
    A2 <- crossprod(sk * sk, W2) + (s2p^2 + ssr)
    dh0 <- pick_grid(-0.5 * sweep(nu * log(A2), 2, ld2, "+"), GRID2)
    colSums(xi * sk / (1 + outer(ell, dh0))) + xp * s2p
  }
  lo <- rep(-30, B); hi <- rep(30, B)
  for (it in 1:40) { mid <- (lo + hi) / 2; up <- numQ(mid) < 0
    lo[up] <- mid[up]; hi[!up] <- mid[!up] }
  split <- (lo + hi) / 2
  gr <- (sqrt(5) - 1) / 2
  side <- function(key, dir, cvec) {
    f <- function(g) ev(split + dir * g, key)
    ## golden-section maximisation of f on g in [0, GMAX], one evaluation per step
    a <- rep(0, B); z <- rep(GMAX, B)
    x1 <- z - gr * (z - a); x2 <- a + gr * (z - a); f1 <- f(x1); f2 <- f(x2)
    for (it in seq_len(NGOLD)) {
      mv <- f1 < f2                                  # maximum lies in [x1, z]
      an <- ifelse(mv, x1, a); zn <- ifelse(mv, z, x2)
      x1n <- ifelse(mv, x2, zn - gr * (zn - an)); x2n <- ifelse(mv, an + gr * (zn - an), x1)
      fe <- f(ifelse(mv, x2n, x1n))
      f1n <- ifelse(mv, f2, fe); f2n <- ifelse(mv, fe, f1)
      a <- an; z <- zn; x1 <- x1n; x2 <- x2n; f1 <- f1n; f2 <- f2n
    }
    gpk <- (a + z) / 2; fpk <- f(gpk)
    f0 <- f(rep(0, B)); fend <- f(rep(GMAX, B))
    ## guard: if an endpoint beats the interior search, use it
    use0 <- f0 > fpk; gpk[use0] <- 0; fpk[use0] <- f0[use0]
    useE <- fend > fpk; gpk[useE] <- GMAX; fpk[useE] <- fend[useE]
    ffall <- f(rep(GFALL, B))
    interior <- gpk > 0.05 & gpk < GMAX - 0.05
    lapply(cvec, function(c0) {
      has <- fpk >= c0
      ## rising crossing in [0, gpk]
      lo <- rep(0, B); hi <- gpk; start_in <- f0 >= c0
      for (it in seq_len(NB2)) { mid <- (lo + hi) / 2; up <- f(mid) < c0
        lo[up] <- mid[up]; hi[!up] <- mid[!up] }
      g1 <- ifelse(start_in, 0, (lo + hi) / 2)
      ## falling crossing in [gpk, GFALL]
      lo <- gpk; hi <- rep(GFALL, B); stays <- ffall >= c0
      for (it in seq_len(NB2)) { mid <- (lo + hi) / 2; above <- f(mid) >= c0
        lo[above] <- mid[above]; hi[!above] <- mid[!above] }
      g2 <- ifelse(stays, Inf, (lo + hi) / 2)
      g1[!has] <- NA; g2[!has] <- NA
      list(g1 = g1, g2 = g2, interior = interior)
    })
  }
  prob <- function(R, L, lam) {                      # R/L: side results for one threshold
    pr <- ifelse(is.na(R$g1), 0, pnorm(split + R$g2 - lam) - pnorm(split + R$g1 - lam))
    pl <- ifelse(is.na(L$g1), 0, pnorm(split - L$g1 - lam) - pnorm(split - L$g2 - lam))
    mean(pr + pl)
  }
  tq_ <- qt(1 - AL / 2, df = nu - 1)
  qa  <- tq_ * sqrt(nu / (nu - 1 + tq_^2)); v0 <- TAU2 * S0
  cBF <- -0.5 * log1p(v0) - (nu / 2) * log1p(-v0 * qa^2 / (nu * (1 + v0)))
  specs <- list(score_h0 = list(key = "absQ", c = tq_), signflip = list(key = "sf", c = -log(AL)),
                bf_h0 = list(key = "bf0", c = cBF), mbf = list(key = "mbf", c = cBF))
  res <- list()
  for (mm in names(specs)) {
    sp_ <- specs[[mm]]
    R <- side(sp_$key, +1, sp_$c); L <- side(sp_$key, -1, sp_$c)
    for (k in seq_along(AL)) {
      res[[length(res) + 1]] <- data.frame(tag, method = mm, alpha = AL[k],
        size_ratio = prob(R[[k]], L[[k]], 0) / AL[k], power = prob(R[[k]], L[[k]], tq_[k]),
        lambda = tq_[k], nonmonotone_share = mean(c(R[[k]]$interior, L[[k]]$interior)))
    }
  }
  do.call(rbind, res)
}

## ---- same configurations and seeds as 06 ------------------------------------------
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
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
load_geu <- function(i) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t",
                fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  list(Z = Z, j = match(g$lead_pos, pos))
}
OUT <- "calibration/results/treatments_pilot_fix"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
run_one <- function(z) {
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f", z$n, z$nb, z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    set.seed(1000 * z$n + z$nb)
    Z <- gen_blocks(z$n, z$nb, 50L, 0.95); sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", n = z$n, blocks = z$nb,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(4242 + 1000 * z$n + z$nb + round(100 * z$h2))
    out <- fix_config(sp, z$h2 / (1 - z$h2), tag)
  } else {
    id <- sprintf("geu_%s_h%.2f", GENES$symbol[z$i], z$h2)
    ck <- file.path(OUT, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    L <- load_geu(z$i); sp <- spectral_setup(L$Z, j = L$j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", n = nrow(L$Z), blocks = NA,
                      h2 = z$h2, nu_delta = indices(sp, z$h2 / (1 - z$h2))[["nu_delta"]])
    set.seed(7000 + z$i + round(100 * z$h2))
    out <- fix_config(sp, z$h2 / (1 - z$h2), tag)
  }
  saveRDS(out, ck)
  cat(sprintf("[%7.1f s] %s done\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), id))
  ck
}
st <- mclapply(cfgs, function(z) tryCatch(run_one(z), error = function(e)
  paste("ERROR", conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st)[grepl("ERROR", unlist(st))])
fx <- do.call(rbind, lapply(list.files(OUT, "rds$", full.names = TRUE), readRDS))
old <- read.csv("calibration/results/treatments_pilot.csv")
old <- old[!(old$method %in% unique(fx$method)), ]
old$nonmonotone_share <- NA
final <- rbind(old, fx[, names(old)])
write.csv(final, "calibration/results/treatments_pilot_final.csv", row.names = FALSE)
cat("configs fixed:", length(unique(fx$config)), "of", length(cfgs), "\n")
cat("done in", format(Sys.time() - t0), "\n")
