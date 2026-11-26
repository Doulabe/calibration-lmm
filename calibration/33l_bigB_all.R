# ==============================================================================
# 33l_bigB_all.R
# Refined evaluation with 100,000 draws (200,000 near the maximum) of the protected p-values of
# the Wald and Kenward-Roger tests for all variants entering GEMMA's paths (38 variants), as in
# 33k_heavy_tail_check.R, which treated the four decisive cases. The conditional probabilities of
# these tests are heavy-tailed at intermediate and large delta0, so that 4,000 draws (33i) can
# produce isolated peaks from single draws. The four cases of 33k are reused.
# Same grid of delta0 as 33g/33i/33k, delta_hat on GEMMA's range.
# Usage (from the repository root):
#   nohup caffeinate -i nice -n 10 Rscript calibration/33l_bigB_all.R [cores=3] > log 2>&1 &
# Output: results/protected_check_bigB/<pop>_<gene>_<pos>_<step>_<stat>.rds (format of 33k)
# ==============================================================================
Sys.setenv(PC_SOURCE_ONLY = "1")
args0 <- commandArgs(TRUE)
try(source("calibration/33g_protected_check.R"), silent = TRUE)
stopifnot(exists("setup_case"), exists("gemma_stats"))
NC <- if (length(args0) && grepl("^[0-9]+$", args0[1])) as.integer(args0[1]) else 3L
OUTB <- "calibration/results/protected_check_bigB"; dir.create(OUTB, showWarnings = FALSE, recursive = TRUE)
CH <- 10000L

## step-aware setup (as setup_case2 of 33i), on the variants table of 33e
COMPV <- do.call(rbind, lapply(list.files(COMP, "^(EUR|YRI)_ENSG[0-9]+[.]rds$", full.names = TRUE), function(f) {
  r <- readRDS(f); if (is.null(r$variants)) return(NULL); cbind(pop = r$pop, gene = r$gene, r$variants) }))
setup_step <- function(pop, symbol, pos, step) {
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
  list(sp = sp, obs = list(t = o$t, Fs = o$Fs, lrt = o$lrt, p_kr = 10^v$log10p_kr_full), v = v,
       level = log10(2 * pt(-sqrt((nu_lev - 1) * expm1(pen / n)), nu_lev - 1)))
}
lp_lin <- function(sp, delta0, obs, stat, B) {               # as in 33k
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

## cases: every (variant, step) x {wald, kr}, minus the four cases of 33k
CA <- unique(COMPV[, c("pop", "gene", "pos", "step")])
CA <- do.call(rbind, lapply(seq_len(nrow(CA)), function(i) data.frame(CA[i, ], stat = c("wald", "kr"), row.names = NULL)))
CA$id <- sprintf("%s_%s_%d_%d_%s", CA$pop, gsub("[^A-Za-z0-9]", "", CA$gene), CA$pos, CA$step, CA$stat)
DONE33k <- c("YRI_PEX6region_42915920_1_wald", "YRI_ZNF266_9477868_1_wald", "EUR_PEX6region_42944850_1_wald", "EUR_ERAP2_96252432_1_kr")
CA <- CA[!CA$id %in% DONE33k, ]
CA <- CA[order(CA$pop == "EUR", decreasing = TRUE), ]           # the longer EUR cases first
cat("cases to compute:", nrow(CA), "\n")
if ("--test" %in% args0) {                                   # setup of a few cases, including repeated positions
  for (k in c(1, which(duplicated(CA$pos))[1], nrow(CA))) { cs <- CA[k, ]; S <- setup_step(cs$pop, cs$gene, cs$pos, cs$step)
    cat(cs$id, " nu =", S$sp$nu0, " level =", round(S$level, 2), " t check =",
        round((log(2) + pt(-abs(S$obs$t), S$sp$nu0 - 1, log.p = TRUE)) / log(10), 3), " gemma wald =", round(S$v$log10p_gemma_wald, 3), "\n") }
  quit(save = "no")
}
run <- function(k) {
  cs <- CA[k, ]; ck <- file.path(OUTB, paste0(cs$id, ".rds")); if (file.exists(ck)) return(paste("cached", cs$id))
  set.seed(6000 + k); t0 <- Sys.time()
  S <- setup_step(cs$pop, cs$gene, cs$pos, cs$step)
  prof1 <- t(vapply(D1, function(d0) big_one(S$sp, d0, S$obs, cs$stat, 100000L), numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, function(d0) big_one(S$sp, d0, S$obs, cs$stat, 200000L), numeric(3)))
  prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2)); prof$id <- cs$id
  i <- which.max(prof$log10p); v <- S$v
  info <- data.frame(id = cs$id, pop = cs$pop, gene = cs$gene, pos = cs$pos, step = cs$step, stat = cs$stat,
                     level = S$level, nominal = if (cs$stat == "wald") v$log10p_gemma_wald else v$log10p_kr_full,
                     stored_protected = NA, new_max = prof$log10p[i], new_argmax = prof$delta0[i],
                     controlled = max(prof$log10up), elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(info = info, profile = prof), ck); sprintf("done %s (%.1f min)", cs$id, info$elapsed_min)
}
print(unlist(parallel::mclapply(seq_len(nrow(CA)), function(k) tryCatch(run(k), error = function(e) paste("FAIL", CA$id[k], conditionMessage(e))),
                                mc.cores = NC, mc.preschedule = FALSE)))
cat("bigB files:", length(list.files(OUTB, "[.]rds$")), "\n")
