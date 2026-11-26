# ==============================================================================
# 27a_prepare_yri_windows.R
# YRI genotype windows for 27_geuvadis_application.R, from the official GEUVADIS genotype
# release (local VCF subsets of chromosomes 1-7). bcftools cannot subset samples in these
# files (undeclared FORMAT tags), so all 465 samples are queried by the shell step below and
# the 89 YRI columns are kept here, in the order of data/geuvadis/geuvadis_yri_overlap.txt.
# SNPs with a missing YRI genotype are dropped. Same windows as the EUR files (27e).
#   shell step (per window): bcftools view -r chr:lo-hi -m2 -M2 -v snps VCF |
#                            bcftools query -f '%POS[\t%GT]\n' > all_<tag>.txt
# Usage: Rscript calibration/27a_prepare_yri_windows.R <dir with all_<tag>.txt> <vcf_samples.txt>
# Output: data/geuvadis/loci_yri/<tag>.gt (POS + 89 GT columns)
# ==============================================================================
a <- commandArgs(TRUE); IN <- a[1]; vs <- readLines(a[2])
yri <- readLines("data/geuvadis/geuvadis_yri_overlap.txt"); col <- match(yri, vs); stopifnot(!anyNA(col))
for (f in list.files(IN, "^all_ENSG.*\\.txt$", full.names = TRUE)) {
  tag <- sub("^all_(ENSG[0-9]+)\\.txt$", "\\1", basename(f))
  L <- strsplit(readLines(f), "\t", fixed = TRUE)
  pos <- vapply(L, `[`, "", 1L); G <- t(vapply(L, function(r) r[-1L][col], character(length(yri))))
  ok <- rowSums(G == "." | G == "./." | G == ".|.") == 0
  out <- paste(pos[ok], apply(G[ok, , drop = FALSE], 1, paste, collapse = "\t"), sep = "\t")
  writeLines(out, file.path("data/geuvadis/loci_yri", paste0(tag, ".gt")))
  cat(sprintf("%s: %d SNP rows, %d kept (no missing YRI genotype), %d dropped\n", tag, length(L), sum(ok), sum(!ok)))
}
