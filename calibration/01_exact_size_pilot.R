# Driver of the pilot; functions in 00_functions.R
args   <- commandArgs(TRUE)
B_DRAW <- if (length(args) >= 1) as.integer(args[1]) else 10000L
source("calibration/00_functions.R")

H2   <- c(0.05, 0.20)
res  <- list()
t0   <- Sys.time()
for (n in c(200L, 1000L)) for (nb in c(2L, 5L, 20L, 100L)) {
  set.seed(1000 * n + nb)
  Z  <- gen_blocks(n, nb, 50L, 0.95)
  xin <- as.numeric(scale(rnorm(n)))                     # control: unlinked candidate
  sp_ld  <- spectral_setup(Z, j = 25L)                   # candidate inside block 1
  sp_ind <- spectral_setup(Z, x = xin)                   # candidate outside the region
  for (h2 in H2) {
    d0 <- h2 / (1 - h2)
    for (ct in c("in_LD", "unlinked")) {
      sp <- if (ct == "in_LD") sp_ld else sp_ind
      res[[length(res) + 1]] <- run_config(sp, d0,
        list(design = "AR1_rho0.95", n = n, blocks = nb, m = 50L * nb,
             h2 = h2, candidate = ct))
      cat(sprintf("[%6.1f s] n=%4d blocks=%3d h2=%.2f %-8s done\n",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")), n, nb, h2, ct))
    }
  }
}
syn <- do.call(rbind, lapply(res, as.data.frame))
write.csv(syn, file.path(OUT, "exact_size_synthetic.csv"), row.names = FALSE)

## ----------------------------------------------------------------- GEUVADIS windows
D <- "data/geuvadis"; LOC <- file.path(D, "loci")
GENES <- data.frame(
  tag = c("ENSG00000164308", "ENSG00000197728", "ENSG00000166750", "ENSG00000203875",
          "ENSG00000198468", "ENSG00000124587", "ENSG00000230658", "ENSG00000174652"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region",
             "TRA2A-AS", "ZNF266"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L), stringsAsFactors = FALSE)
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
res2 <- list(); perSNP <- list()
for (i in seq_len(nrow(GENES))) {
  g  <- GENES[i, ]
  sp_ <- strsplit(readLines(file.path(LOC, paste0(g$tag, ".gt"))), "\t", fixed = TRUE)
  pos <- as.integer(vapply(sp_, `[`, "", 1L))
  G  <- vapply(sp_, function(r) gt2d(r[-1L]), integer(length(sp_[[1]]) - 1L))
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf)
  keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  G <- G[, keep, drop = FALSE]; pos <- pos[keep]
  Z <- scale(G); Z[!is.finite(Z)] <- 0
  jl <- match(g$lead_pos, pos)
  if (is.na(jl)) { cat("lead not found for", g$symbol, "\n"); next }
  sp <- spectral_setup(Z, j = jl)
  for (h2 in H2) {
    d0 <- h2 / (1 - h2)
    res2[[length(res2) + 1]] <- run_config(sp, d0,
      list(design = "GEUVADIS_cis", gene = g$symbol, n = nrow(Z), m = ncol(Z),
           h2 = h2, candidate = "lead"))
  }
  ## effective df for every SNP of the window (full-kernel approximation)
  spf <- spectral_setup(Z, x = rep(0, nrow(Z)) + 0)      # placeholder to get U, ell
  Uf <- spf$U; ellf <- spf$ell
  XI <- crossprod(Uf, sweep(Z, 2, colMeans(Z)))           # r x m loadings
  tot <- colSums(Z^2)
  for (h2 in H2) {
    d0 <- h2 / (1 - h2); w <- 1 / (1 + d0 * ellf); a <- ellf * w
    aU <- sum(a) / spf$nu0; vU <- sum(a^2) / spf$nu0 - aU^2
    S2 <- colSums(XI^2 * w) + pmax(tot - colSums(XI^2), 0)
    ap <- colSums(XI^2 * w * a) / S2
    nud <- spf$nu0 * vU / (ap - aU)^2
    perSNP[[length(perSNP) + 1]] <- data.frame(gene = g$symbol, h2 = h2, pos = pos,
      nu_delta = nud, is_lead = seq_along(pos) == jl)
  }
  cat(sprintf("[%6.1f s] GEUVADIS %s (m=%d) done\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), g$symbol, ncol(Z)))
}
geu <- do.call(rbind, lapply(res2, as.data.frame))
write.csv(geu, file.path(OUT, "exact_size_geuvadis.csv"), row.names = FALSE)
ps  <- do.call(rbind, perSNP)
write.csv(ps, file.path(OUT, "nu_delta_per_snp_geuvadis.csv"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "\n")
