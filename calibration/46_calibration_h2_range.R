# ==============================================================================
# 46_calibration_h2_range.R
# The exact calibration study of 36_calibration_enlarged.R repeated at a true
# background heritability h2 = 0 (no regional background, delta0 = 0, the boundary of the
# parameter space) and h2 = 0.5 (delta0 = 1). The candidates are those of 36, so that each is
# compared with itself: h2 = 0 uses the ten candidates per GEUVADIS window chosen at h2 = 0.05,
# h2 = 0.5 those chosen at h2 = 0.20, and the seven AR(1) designs are the same. Each configuration
# runs calib_config() of 36 (Wald, Kenward-Roger, Satterthwaite, Edgeworth, likelihood ratio,
# score Q, test with delta known) and calib_score() of 36b (t(delta_hat_0) of EMMAX) on the same
# draws of the error contrasts (same seed before each), with B = 1000 draws.
# Usage (from the repository root): nice -n 10 Rscript calibration/46_calibration_h2_range.R [cores=3]
# Output: results/calib_h2range/<config>.rds, results/calib_h2range.csv, results/calib_h2range_extra.csv
# ==============================================================================
src <- readLines("calibration/36_calibration_enlarged.R")
eval(parse(text = src[seq_len(grep("^OUT <- ", src)[1] - 1L)]))   # functions, GENES, load_geu, cfgs of 36
s36b <- readLines("calibration/36b_score_reference.R")
eval(parse(text = s36b[grep("^calib_score <- ", s36b)[1]:(grep("^OUT2 <- ", s36b)[1] - 1L)]))
a0 <- commandArgs(TRUE)
NCORES <- if (length(a0) && grepl("^[0-9]+$", a0[1])) as.integer(a0[1]) else 3L
H2MAP <- c("0.05" = 0, "0.20" = 0.5)                        # design h2 of 36 -> new h2
cfgs <- lapply(cfgs, function(z) { z$h2 <- unname(H2MAP[sprintf("%.2f", z$h2)]); z })
stopifnot(all(vapply(cfgs, function(z) z$h2 %in% c(0, 0.5), TRUE)))
OUT_H2 <- "calibration/results/calib_h2range"
dir.create(OUT_H2, showWarnings = FALSE, recursive = TRUE)
run_h2 <- function(z) {                                      # set-up and seeds as run_one() of 36
  d0 <- z$h2 / (1 - z$h2)
  if (z$kind == "syn") {
    id <- sprintf("syn_n%d_b%d_h%.2f%s", z$n, z$nb, z$h2, if (z$unlinked) "_unlinked" else "")
    ck <- file.path(OUT_H2, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    set.seed(1000 * z$n + z$nb); Z <- gen_blocks(z$n, z$nb, 50L, 0.95)
    if (z$unlinked) { set.seed(77 + z$n); sp <- spectral_setup(Z, x = rnorm(z$n)) } else sp <- spectral_setup(Z, j = 25L)
    tag <- data.frame(config = id, design = "AR1_rho0.95", gene = NA, pos = NA, role = if (z$unlinked) "unlinked" else "block1",
                      n = z$n, blocks = z$nb, h2 = z$h2, nu_delta = indices(sp, d0)[["nu_delta"]])
    seed <- 4242 + 1000 * z$n + z$nb + round(100 * z$h2) + if (z$unlinked) 17L else 0L
  } else {
    g <- GENES[z$i, ]; id <- sprintf("geu_%s_%d_h%.2f", g$symbol, z$pos, z$h2)
    ck <- file.path(OUT_H2, paste0(id, ".rds")); if (file.exists(ck)) return(ck)
    W <- load_geu(z$i); j <- match(z$pos, W$pos); stopifnot(!is.na(j))
    sp <- spectral_setup(W$Z, j = j)
    tag <- data.frame(config = id, design = "GEUVADIS_cis", gene = g$symbol, pos = z$pos, role = z$role,
                      n = nrow(W$Z), blocks = NA, h2 = z$h2, nu_delta = indices(sp, d0)[["nu_delta"]])
    seed <- 9000 + 37 * z$i + round(100 * z$h2) + (z$pos %% 9973)
  }
  set.seed(seed); out <- calib_config(sp, d0, tag)
  set.seed(seed); sc <- calib_score(sp, d0, tag)            # same draws: emmax_t, and score_t as a check
  ex <- attr(out, "extra")
  out <- rbind(out, sc[sc$method == "emmax_t", names(out)])
  chk <- merge(sc[sc$method == "score_t", c("alpha", "size_ratio")], out[out$method == "score_t", c("alpha", "size_ratio")], by = "alpha")
  ex$score_check <- max(abs(chk$size_ratio.x / chk$size_ratio.y - 1), na.rm = TRUE)
  attr(out, "extra") <- ex
  saveRDS(out, ck); ck
}
cfgs <- c(cfgs[vapply(cfgs, function(z) z$kind == "syn", TRUE)], cfgs[vapply(cfgs, function(z) z$kind == "geu", TRUE)])
cat("configurations:", length(cfgs), "\n")
st <- mclapply(cfgs, function(z) tryCatch(run_h2(z), error = function(e) paste("ERROR", conditionMessage(e))),
               mc.cores = NCORES, mc.preschedule = FALSE)
print(unlist(st)[grepl("ERROR", unlist(st))])
fs <- list.files(OUT_H2, "rds$", full.names = TRUE); L <- lapply(fs, readRDS)
write.csv(do.call(rbind, L), "calibration/results/calib_h2range.csv", row.names = FALSE)
X <- do.call(rbind, lapply(L, attr, "extra"))
write.csv(X, "calibration/results/calib_h2range_extra.csv", row.names = FALSE)
cat("configs done:", length(fs), "of", length(cfgs), " max score check:", signif(max(X$score_check, na.rm = TRUE), 3), "\n")
## summary by class of nu_delta (as Table 2 of the article)
C <- do.call(rbind, L); C$cls <- cut(C$nu_delta, c(0, 20, 100, Inf), right = FALSE)
options(width = 200)
for (h in c(0, 0.5)) {
  cat("\n== h2 =", h, ": configurations per class\n"); print(table(unique(C[C$h2 == h, c("config", "cls")])$cls))
  print(aggregate(size_ratio ~ method + alpha + cls, C[C$h2 == h & C$method %in% c("oracle_t", "wald_t", "emmax_t", "lrt_chi2", "kr", "satterthwaite", "wald_edgeworth"), ],
                  function(v) signif(median(v), 3)))
}
