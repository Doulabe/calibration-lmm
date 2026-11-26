# ==============================================================================
# 58_geuvadis_application_int.R
# Sensitivity of the GEUVADIS application to the distribution of the phenotype. Reruns
# 27c_geuvadis_application_cov.R (same windows, covariates, kernels, paths, levels and seeds)
# after a rank-based inverse normal transformation of the expression of each gene within each
# sample, y_i = qnorm{(rank_i - 1/2) / n}, in place of its standardization.
# Usage (from the repository root): nice -n 10 Rscript calibration/58_geuvadis_application_int.R [cores=1] [--pop=EUR|YRI|ALL]
# Output: results/geuvadis_app_int/<pop>_<tag>.rds, results/geuvadis_app_int_{paths,variants}.csv
# ==============================================================================
src <- readLines("calibration/27c_geuvadis_application_cov.R")
old <- "y <- as.numeric(scale(EXPR[[g$gene]][W$samples])); stopifnot(!anyNA(y))"
k <- grep(old, src, fixed = TRUE); stopifnot(length(k) == 1L)
src[k] <- sub(old, "y0 <- EXPR[[g$gene]][W$samples]; stopifnot(!anyNA(y0)); y <- qnorm((rank(y0) - 0.5) / length(y0))", src[k], fixed = TRUE)
stopifnot(sum(grepl("geuvadis_app_cov", src[-seq_len(30)], fixed = TRUE)) >= 3L)
src <- gsub("geuvadis_app_cov", "geuvadis_app_int", src, fixed = TRUE)
eval(parse(text = src))
