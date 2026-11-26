# ==============================================================================
# 04_variance_check.R
# Exact variance of the plug-in Wald statistic t(delta_hat) under H0, from the
# conditional representation Var(t) = E[m^2 + v] (m, v as in 01_exact_size_pilot),
# relative to the Student reference (nu0-1)/(nu0-3), compared with the
# first-order prediction 1 + kappa + 1/nu_delta (delta method of the initial
# sketches, done with REML and the exact Kackar-Harville decomposition).
# Usage (from the repository root):
#   nice -n 19 Rscript calibration/04_variance_check.R
# Output: calibration/results/variance_check.csv
# ==============================================================================
B_DRAW <- 10000L
source("calibration/00_functions.R")

mv_draws <- function(sp, delta0, B) {
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; nu0 <- sp$nu0; r <- sp$r
  Wg  <- 1 / (1 + outer(ell, GRID)); ldg <- colSums(log(1 + outer(ell, GRID)))
  S2g <- colSums(xi^2 * Wg) + xp^2
  d0  <- 1 + delta0 * ell
  res <- list()
  for (ch in seq_len(ceiling(B / CHUNK))) {
    b <- CHUNK
    s <- matrix(rnorm(r * b), r, b) * sqrt(d0)
    sp_ <- if (sp$nperp) rnorm(b) else numeric(b)
    ssr <- if (sp$nrest > 0) rchisq(b, sp$nrest) else numeric(b)
    s2 <- s * s
    Ag <- crossprod(s2, Wg) + (sp_^2 + ssr)
    Bg <- crossprod(xi * s, Wg) + xp * sp_
    O  <- -0.5 * (sweep((nu0 - 1) * log(Ag - sweep(Bg^2, 2, S2g, "/")), 2,
                        ldg + log(S2g), "+"))
    g <- max.col(O, ties.method = "first"); dh <- GRID[g]
    int <- g > 2 & g < length(GRID)
    if (any(int)) {
      ii <- which(int); lg <- log(GRID)
      x1 <- lg[g[ii] - 1]; x2 <- lg[g[ii]]; x3 <- lg[g[ii] + 1]
      f1 <- O[cbind(ii, g[ii] - 1)]; f2 <- O[cbind(ii, g[ii])]; f3 <- O[cbind(ii, g[ii] + 1)]
      den <- (x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1)
      num <- (x2 - x1)^2 * (f2 - f3) - (x2 - x3)^2 * (f2 - f1)
      xs <- ifelse(abs(den) > 1e-14, x2 - 0.5 * num / den, x2)
      dh[ii] <- exp(pmin(pmax(xs, x1), x3))
    }
    at <- function(dl) {
      iD <- 1 / (1 + outer(ell, dl))
      A  <- colSums(s2 * iD) + sp_^2 + ssr
      Bx <- colSums(xi * s * iD) + xp * sp_
      S2 <- colSums(xi^2 * iD) + xp^2
      list(Bx = Bx, S2 = S2, RSS = A - Bx^2 / S2)
    }
    Q <- at(dh); Q0 <- at(rep(delta0, b))
    sig2 <- Q$RSS / (nu0 - 1)
    m <- (Q$Bx / Q$S2 - Q0$Bx / Q0$S2) * sqrt(Q$S2 / sig2)
    res[[ch]] <- data.frame(m = m, v = Q$S2 / (Q0$S2 * sig2))
  }
  do.call(rbind, res)
}

vcheck <- function(sp, delta0, tag) {
  D  <- mv_draws(sp, delta0, B_DRAW)
  ix <- indices(sp, delta0)
  tv <- (sp$nu0 - 1) / (sp$nu0 - 3)
  c(tag, list(nu_delta = ix[["nu_delta"]], kappa = ix[["kappa"]],
    var_exact_rel = mean(D$m^2 + D$v) / tv,
    var_first_order_rel = 1 + ix[["kappa"]] + 1 / ix[["nu_delta"]]))
}

res <- list()
for (n in c(200L, 1000L)) for (nb in c(2L, 5L, 20L, 100L)) {
  set.seed(1000 * n + nb)
  Z  <- gen_blocks(n, nb, 50L, 0.95)
  sp <- spectral_setup(Z, j = 25L)
  for (h2 in c(0.05, 0.20)) res[[length(res) + 1]] <- vcheck(sp, h2 / (1 - h2),
    list(design = "AR1_rho0.95", n = n, blocks = nb, h2 = h2, candidate = "in_LD"))
}
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
for (i in seq_len(nrow(GENES))) {
  g <- GENES[i, ]
  L <- strsplit(readLines(file.path("data/geuvadis/loci", paste0(g$tag, ".gt"))), "\t",
                fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  sp <- spectral_setup(Z, j = match(g$lead_pos, pos))
  set.seed(91 + i)
  for (h2 in c(0.05, 0.20)) res[[length(res) + 1]] <- vcheck(sp, h2 / (1 - h2),
    list(design = "GEUVADIS_cis", n = nrow(Z), blocks = NA, h2 = h2,
         candidate = paste0("lead_", g$symbol)))
}
out <- do.call(rbind, lapply(res, as.data.frame))
write.csv(out, "calibration/results/variance_check.csv", row.names = FALSE)
print(out[, c("n", "blocks", "h2", "candidate", "nu_delta", "kappa", "var_exact_rel",
              "var_first_order_rel")], digits = 3)
