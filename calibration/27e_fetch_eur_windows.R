# ==============================================================================
# 27e_fetch_eur_windows.R
# European genotype windows of the eight GEUVADIS genes for 27_geuvadis_application.R and
# 27c_geuvadis_application_cov.R. For each gene the cis window (TSS +/- 500 kb, extended to
# cover the published lead eQTL by 100 kb) is streamed from the 1000 Genomes phase-3 release
# (EBI FTP) with bcftools, for the 358 samples present in both GEUVADIS and the phase-3
# European panel (data/geuvadis/geuvadis_eur_overlap.txt), biallelic SNPs only. The filter on
# the minor allele frequency is applied by the analysis scripts.
# Inputs: data/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz (TSS of each gene)
#         data/geuvadis/geuvadis_eur_overlap.txt
# Output: data/geuvadis/loci/<tag>.gt (POS + 358 GT columns), <tag> the Ensembl gene identifier
# Usage (from the repository root): Rscript calibration/27e_fetch_eur_windows.R
# ==============================================================================
D       <- "data/geuvadis"
LOC_DIR <- file.path(D, "loci"); dir.create(LOC_DIR, showWarnings = FALSE, recursive = TRUE)
SAMP    <- file.path(D, "geuvadis_eur_overlap.txt")
EXPR_GZ <- file.path(D, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz")
URLB    <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
W       <- 500000L

GENES <- data.frame(
  gene   = c("ENSG00000164308.12", "ENSG00000197728.5", "ENSG00000166750.4",
             "ENSG00000203875.4",  "ENSG00000198468.2", "ENSG00000124587.9",
             "ENSG00000230658.1",  "ENSG00000174652.12"),
  symbol = c("ERAP2", "RPS26", "SLFN5", "SNHG5", "FLVCR1-AS1(1q32)",
             "PEX6-region", "TRA2A-AS(7p15)", "ZNF266"),
  chr      = c("5", "12", "17", "6", "1", "6", "7", "19"),
  lead_pos = c(96252589L, 56401085L, 33571546L, 86387888L,
               213049214L, 42944850L, 23143113L, 9544276L),
  stringsAsFactors = FALSE)

## transcription start site of each gene, fourth column of the expression file
expr_rows <- list()
con <- gzfile(EXPR_GZ); open(con); readLines(con, n = 1L)
while (length(l <- readLines(con, n = 4000L))) {
  for (g in GENES$gene) {
    hit <- grep(paste0("^", g, "\t"), l)
    if (length(hit)) expr_rows[[g]] <- strsplit(l[hit[1]], "\t")[[1]]
  }
  if (length(expr_rows) == nrow(GENES)) break
}
close(con)
stopifnot(length(expr_rows) == nrow(GENES))
tss <- vapply(expr_rows, function(r) as.integer(r[4]), integer(1))[GENES$gene]

stream_locus <- function(chr, lo, hi, tag) {
  gt <- file.path(LOC_DIR, paste0(tag, ".gt"))
  url <- sprintf(URLB, chr)
  for (try in 1:3) {
    if (file.exists(gt) && file.size(gt) > 0) break
    cmd <- sprintf(paste0(
      "bcftools view -r %s:%d-%d -S %s --force-samples -m2 -M2 -v snps %s 2>/dev/null",
      " | bcftools query -f '%%POS[\\t%%GT]\\n' > %s 2>/dev/null"),
      chr, lo, hi, shQuote(SAMP), shQuote(url), shQuote(gt))
    system(cmd)
    if (!file.exists(gt) || file.size(gt) == 0) Sys.sleep(10)
  }
  if (!file.exists(gt) || file.size(gt) == 0) return(NULL)
  gt
}

for (i in seq_len(nrow(GENES))) {
  g <- GENES[i, ]; tag <- sub("\\..*$", "", g$gene)
  lo <- min(tss[g$gene] - W, g$lead_pos - 100000L)
  hi <- max(tss[g$gene] + W, g$lead_pos + 100000L)
  gt <- stream_locus(g$chr, lo, hi, tag)
  cat(sprintf("%-18s chr%s:%d-%d  %s\n", g$symbol, g$chr, lo, hi, if (is.null(gt)) "FAILED" else gt))
}
