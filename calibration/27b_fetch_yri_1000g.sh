#!/bin/bash
# 27b_fetch_yri_1000g.sh
# YRI genotype windows for the three genes outside chromosomes 1-7 (RPS26, SLFN5, ZNF266), streamed from the
# 1000 Genomes phase-3 release (EBI FTP) as in 27e_fetch_eur_windows.R: biallelic SNPs of
# the same windows as the EUR files (data/geuvadis/loci/<tag>.gt), 87 GEUVADIS YRI samples present in phase 3
# (data/geuvadis/geuvadis_yri87_1000g.txt). The sample order of each output file is written to <tag>.samples.
# Run from the repository root (the remote .tbi indexes are cached there).
set -u
D=data/geuvadis; OUT=$D/loci_yri; S=$D/geuvadis_yri87_1000g.txt
URLB="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
for spec in "ENSG00000197728:12" "ENSG00000166750:17" "ENSG00000174652:19"; do
  tag=${spec%%:*}; chr=${spec##*:}; url=$(printf "$URLB" "$chr")
  lo=$(head -1 $D/loci/$tag.gt | cut -f1); hi=$(tail -1 $D/loci/$tag.gt | cut -f1)
  for try in 1 2 3; do
    bcftools view -r $chr:$lo-$hi -S $S --force-samples -m2 -M2 -v snps "$url" 2>/dev/null \
      | bcftools query -f '%POS[\t%GT]\n' > $OUT/$tag.gt 2>/dev/null
    [ -s $OUT/$tag.gt ] && break; sleep 10
  done
  bcftools query -l -S $S --force-samples "$url" 2>/dev/null > $OUT/$tag.samples
  echo "$tag chr$chr:$lo-$hi -> $(wc -l < $OUT/$tag.gt) SNPs, $(wc -l < $OUT/$tag.samples) samples"
done
