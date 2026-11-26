# ==============================================================================
# 33f_export_geuvadis_cov.R
# As 33a_export_geuvadis.R, with the expression and the standardized genotypes residualized
# on the intercept and the sample covariates (sex; FIN, TSI, CEU indicators for EUR), for
# SuSiE-inf, which takes no covariates: GEU_RESDIR=geuvadis_comp_cov python3 33b_susieinf_geuvadis.py
# Output: results/geuvadis_comp_cov/data/<pop>_<symbol>_{Z,y,pos}.txt
# Usage (from the repository root): Rscript calibration/33f_export_geuvadis_cov.R
# ==============================================================================
src <- readLines("calibration/33a_export_geuvadis.R")
stop_at <- grep("^OUT <- ", src)[1] - 1L
eval(parse(text = src[seq_len(stop_at)]))
COV <- read.csv(file.path(D, "geuvadis_covariates.csv"), stringsAsFactors = FALSE)
covariates <- function(pop, samples) {
  d <- COV[match(samples, COV$sample), ]; stopifnot(!anyNA(d$sex))
  C <- cbind(male = as.numeric(d$sex == "male"))
  if (pop == "EUR") C <- cbind(C, FIN = as.numeric(d$ancestry == "FIN"), TSI = as.numeric(d$ancestry == "TSI"),
                               CEU = as.numeric(d$ancestry == "CEU"))
  C
}
OUT <- "calibration/results/geuvadis_comp_cov/data"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
for (pop in names(POPS)) for (i in seq_len(nrow(GENES))) {
  f <- file.path(POPS[[pop]]$dir, paste0(GENES$tag[i], ".gt")); if (!file.exists(f)) next
  L <- strsplit(readLines(f), "\t", fixed = TRUE)
  pos <- as.integer(vapply(L, `[`, "", 1L))
  G <- vapply(L, function(r) gt2d(r[-1L]), integer(length(L[[1]]) - 1L))
  sf <- file.path(POPS[[pop]]$dir, paste0(GENES$tag[i], ".samples"))
  samples <- if (file.exists(sf)) readLines(sf) else POPS[[pop]]$samples
  maf <- colMeans(G) / 2; maf <- pmin(maf, 1 - maf); keep <- maf >= 0.05 & apply(G, 2, sd) > 0
  X <- scale(G[, keep, drop = FALSE]); X[!is.finite(X)] <- 0
  y <- as.numeric(scale(EXPR[[GENES$gene[i]]][samples]))
  qf <- qr(cbind(1, covariates(pop, samples)))
  Xr <- scale(qr.resid(qf, X)); Xr[!is.finite(Xr)] <- 0; yr <- as.numeric(scale(qr.resid(qf, y)))
  pre <- file.path(OUT, sprintf("%s_%s", pop, GENES$symbol[i]))
  write.table(Xr, paste0(pre, "_Z.txt"), row.names = FALSE, col.names = FALSE)
  writeLines(format(yr, digits = 17), paste0(pre, "_y.txt")); writeLines(as.character(pos[keep]), paste0(pre, "_pos.txt"))
  cat(pop, GENES$symbol[i], nrow(Xr), ncol(Xr), "\n")
}
