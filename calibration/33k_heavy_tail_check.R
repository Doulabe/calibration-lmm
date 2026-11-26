# ==============================================================================
# 33k_heavy_tail_check.R
# For the Wald and Kenward-Roger tests, whose conditional probabilities given the error
# contrasts are heavy-tailed at intermediate and large delta0 (a single draw can dominate the
# average, relative standard errors near one), the refined evaluation of 33g is repeated with
# 25 times more draws: stage 1 with 100,000 draws (in chunks of 10,000) at the same grid of
# delta0, stage 2 with 200,000 draws near the maximum. Cheap, because these tests are
# conditionally Gaussian given the contrasts.
# Cases: the decisive Wald and Kenward-Roger cases of 33g.
# Usage (from the repository root): nice -n 10 Rscript calibration/33k_heavy_tail_check.R [cores=4]
# Output: results/protected_check_bigB/<case>.rds (same format as 33g)
# ==============================================================================
Sys.setenv(PC_SOURCE_ONLY = "1")
args0 <- commandArgs(TRUE)
try(source("calibration/33g_protected_check.R"), silent = TRUE)
stopifnot(exists("setup_case"))
NC <- if (length(args0) && grepl("^[0-9]+$", args0[1])) as.integer(args0[1]) else 4L
OUTB <- "calibration/results/protected_check_bigB"; dir.create(OUTB, showWarnings = FALSE, recursive = TRUE)
CH <- 10000L
lp_lin <- function(sp, delta0, obs, stat, B) {               # log conditional probabilities (Wald or KR), as in exact_one
  ell <- sp$ell; xi <- sp$xi; xp <- sp$xp; r <- sp$r
  d0 <- 1 + delta0 * ell; S0 <- sum(xi^2 / d0) + xp^2
  s <- matrix(rnorm(r * B), r, B) * sqrt(d0)
  spp <- if (sp$nperp) rnorm(B) else numeric(B); ssr <- if (sp$nrest > 0) rchisq(B, sp$nrest) else numeric(B)
  b <- (colSums((xi / d0) * s) + xp * spp) / S0
  Rk <- s - outer(xi, b); Rp <- spp - xp * b; rm(s)
  ev <- function(u) { bv <- u / sqrt(S0); gemma_stats(sp, Rk + outer(xi, bv), Rp + xp * bv, ssr) }
  z0 <- ev(rep(0, B)); zz <- ev(rep(1, B)); wa <- zz$t - z0$t; wg <- z0$t
  c0 <- if (stat == "wald") rep(abs(obs$t), B) else {
    kf <- kr_factors(ell, xi^2, xp^2, sp$nperp + sp$nrest, z0$dh1)
    qt(obs$p_kr / 2, df = kf$df, lower.tail = FALSE) * sqrt(kf$ratio) }
  a1 <- pnorm(-(c0 - wg) / wa, log.p = TRUE); a2 <- pnorm((-c0 - wg) / wa, log.p = TRUE)
  pmax(a1, a2) + log1p(exp(-abs(a1 - a2)))
}
big_one <- function(sp, d0, obs, stat, B) mc_summary(unlist(lapply(seq_len(ceiling(B / CH)), function(i) lp_lin(sp, d0, obs, stat, CH))))
ids <- c("PEX6_YRI_wald", "ZNF266_YRI_wald", "PEX6_EUR_wald", "ERAP2_EUR_kr")
run <- function(id) {
  ck <- file.path(OUTB, paste0(id, ".rds")); if (file.exists(ck)) return(paste("cached", id))
  k <- match(id, CASES$id); cs <- CASES[k, ]; set.seed(5100 + k); t0 <- Sys.time()
  S <- setup_case(cs)
  prof1 <- t(vapply(D1, function(d0) big_one(S$sp, d0, S$obs, cs$stat, 100000L), numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, function(d0) big_one(S$sp, d0, S$obs, cs$stat, 200000L), numeric(3)))
  prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2)); prof$id <- id
  i <- which.max(prof$log10p)
  res <- list(info = cbind(S$info, new_max = prof$log10p[i], new_argmax = prof$delta0[i], controlled = max(prof$log10up),
                           at_1e4 = prof1[D1 == 1e4, "log10p"], at_1e10 = prof1[D1 == 1e10, "log10p"],
                           elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))), profile = prof)
  saveRDS(res, ck); sprintf("done %s (%.1f min)", id, res$info$elapsed_min)
}
print(unlist(parallel::mclapply(ids, function(id) tryCatch(run(id), error = function(e) paste("FAIL", id, conditionMessage(e))),
                                mc.cores = NC, mc.preschedule = FALSE)))
for (id in ids) { r <- readRDS(file.path(OUTB, paste0(id, ".rds")))
  print(r$info[, c("id", "level", "nominal", "stored_protected", "new_max", "new_argmax", "controlled")], digits = 3) }
