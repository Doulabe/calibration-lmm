# ==============================================================================
# 33d_geuvadis_covariates.R
# Sample covariates of the GEUVADIS application, from the ArrayExpress sample table
# E-GEUV-1.sdrf.txt (path in GEUVADIS_SDRF): sex, ancestry category (British, Finnish,
# Tuscan, Utah, Yoruba) and sequencing laboratory. The analyses with covariates
# (27c_geuvadis_application_cov.R, 33e_geuvadis_comparators_cov.R) use sex for both
# samples and three ancestry indicators (FIN, TSI, CEU; GBR as reference) for EUR; the
# laboratory effects are absorbed by the PEER factors removed from the expression.
# Output: data/geuvadis/geuvadis_covariates.csv
# Usage (from the repository root): Rscript calibration/33d_geuvadis_covariates.R
# ==============================================================================
f <- Sys.getenv("GEUVADIS_SDRF", "data/geuvadis/E-GEUV-1.sdrf.txt")
s <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
u <- unique(s[, c("Source Name", "Characteristics[sex]", "Characteristics[ancestry category]", "Characteristics[laboratory]")])
names(u) <- c("sample", "sex", "ancestry", "lab")
stopifnot(!anyDuplicated(u$sample))
u$ancestry <- c(British = "GBR", Finnish = "FIN", Tuscan = "TSI", Utah = "CEU", Yoruba = "YRI")[u$ancestry]
stopifnot(!anyNA(u$ancestry))
write.csv(u, "data/geuvadis/geuvadis_covariates.csv", row.names = FALSE)
for (p in c("eur", "yri")) {
  smp <- readLines(sprintf("data/geuvadis/geuvadis_%s_overlap.txt", p)); d <- u[match(smp, u$sample), ]
  stopifnot(!anyNA(d$sex)); cat(toupper(p), length(smp), "samples\n"); print(table(d$ancestry, d$sex))
}
