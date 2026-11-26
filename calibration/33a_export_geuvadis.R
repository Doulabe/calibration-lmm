# ==============================================================================
# 33a_export_geuvadis.R
# Exports the standardized genotypes and expression of the GEUVADIS application
# (8 cis windows, EUR and YRI; as in 27_geuvadis_application.R) for the Python
# comparators of 33b_susieinf_geuvadis.py.
# Output: results/geuvadis_comp/data/<pop>_<symbol>_{Z,y}.txt and _pos.txt
# Usage (from the repository root): Rscript calibration/33a_export_geuvadis.R
# ==============================================================================
GENES <- data.frame(
  gene = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4", "ENSG00000203875.4",
           "ENSG00000198468.2", "ENSG00000124587.9", "ENSG00000230658.1", "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1", "PEX6-region", "TRA2A-AS", "ZNF266"),
  stringsAsFactors = FALSE)
GENES$tag <- sub("\\..*$", "", GENES$gene)
D <- "data/geuvadis"
EXPR <- { hdr <- strsplit(readLines(gzfile(file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")), n = 1L), "\t")[[1]]
  rows <- list(); con <- gzfile(file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")); open(con); readLines(con, n = 1L)
  while (length(l <- readLines(con, n = 4000L))) { for (g in GENES$gene) { h <- grep(paste0("^", g, "\t"), l)
    if (length(h)) rows[[g]] <- strsplit(l[h[1]], "\t")[[1]] }
    if (length(rows) == nrow(GENES)) break }
  close(con)
  lapply(rows, function(r) setNames(as.numeric(r[-(1:4)]), hdr[-(1:4)])) }
POPS <- list(EUR = list(dir = file.path(D, "loci"), samples = readLines(file.path(D, "geuvadis_eur_overlap.txt"))),
             YRI = list(dir = file.path(D, "loci_yri"), samples = readLines(file.path(D, "geuvadis_yri_overlap.txt"))))
gt2d <- function(v) { d <- integer(length(v))
  d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L; d[v %in% c("1|1", "1/1")] <- 2L; d }
OUT <- "calibration/results/geuvadis_comp/data"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
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
  pre <- file.path(OUT, sprintf("%s_%s", pop, GENES$symbol[i]))
  write.table(X, paste0(pre, "_Z.txt"), row.names = FALSE, col.names = FALSE)
  writeLines(format(y, digits = 17), paste0(pre, "_y.txt")); writeLines(as.character(pos[keep]), paste0(pre, "_pos.txt"))
  cat(pop, GENES$symbol[i], nrow(X), ncol(X), "\n")
}
