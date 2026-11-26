# Finite-sample calibration of regional mixed-model association tests

Code of the article *Finite-sample calibration of regional mixed-model association tests with an
estimated variance ratio*, by Kossi N. Doulabe and Lajmi Lakhal-Chaieb (2026).

The scripts compute the exact finite-sample distributions of the Wald, score, likelihood-ratio and
Bayes-factor statistics of a candidate variant when the variance ratio of a regional polygenic
background is estimated from the data under test, the candidate-specific index $\nu_j$, and the
plug-in and protected $p$-values. They run the published mixed-model and fine-mapping software on
the same simulated and real data, and they write every table and figure of the article.

All the code is in the folder `calibration/`. Scripts are numbered in the order of the analysis and are run from the repository root, for
instance `Rscript calibration/36_calibration_enlarged.R 3`. The header of each script gives its
purpose, inputs, outputs and usage.

## Contents

| Path in `calibration/` | Content |
|---|---|
| `NN_*.R`, `NN_*.py`, `NN_*.sh` | analysis scripts |
| `lib/` | functions of the Bayes-factor and conditional-regression stepwise paths, sourced by scripts 14, 27 and 27c |
| `results/` | per-configuration values read by the table and figure scripts (CSV, JSON lines, small RDS files) |
| `make_tables_and_figures.sh` | writes all tables and figures from `results/` |
| `output/` | tables (`tables/`) and figures (`figures/`) of the article, written by `make_tables_and_figures.sh` (not tracked) |
| `note/` | intermediate tables, written by the first table scripts and read by `35_make_article_tables.R` (not tracked) |
| `samples/` | GEUVADIS sample lists (358 European samples, 89 Yoruba samples, the 87 of them in 1000 Genomes phase 3) |
| `comparators/gemma_arm64_debug.patch` | patch used to build GEMMA 0.98.5 on Apple silicon, which disables the floating-point trapping of the debug mode with no effect on the computations (GPL-3, as GEMMA) |

The exact distributions are coded in `00_functions.R` (spectral reduction of the regional model,
conditional normal representation of the Wald statistic given the error contrasts,
Rao-Blackwellized tail probabilities, closed-form $\nu_\delta$), `06_treatments_pilot.R` (exact
critical values of the treatments of the ratio) and `36_calibration_enlarged.R` (calibration of the
Wald, score, likelihood-ratio, Kenward-Roger, Satterthwaite and Edgeworth tests).
`52_protected_design.R` computes the protected critical values, `14_selection_demo.R` and
`27c_geuvadis_application_cov.R` the stepwise paths and their exact $p$-values, and
`33g_protected_check.R` the protected $p$-values of the statistics of GEMMA.

## Requirements

- R 4.4.3 with `data.table`, `jsonlite`, `statmod` and `parallel`. The comparators also need `susieR` 0.14.2,
  `MM4LMM` 3.0.3, `rMVP` 1.4.6 (FarmCPU) with `bigmemory`, `GridLMM`
  ([deruncie/GridLMM](https://github.com/deruncie/GridLMM)), and MLMM 0.1.1
  ([MultLocMixMod](https://github.com/Gregor-Mendel-Institute/MultLocMixMod)) with EMMA 1.1.2.
- GEMMA 0.98.5 ([genetics-statistics/GEMMA](https://github.com/genetics-statistics/GEMMA)), also
  used for BSLMM, and QTLtools 1.3.1 ([qtltools](https://github.com/qtltools/qtltools)).
- Python 3 with `numpy`, `pandas` and `scipy`, and fine-mapping-inf
  ([FinucaneLab/fine-mapping-inf](https://github.com/FinucaneLab/fine-mapping-inf)) for
  SuSiE-inf and FINEMAP-inf.
- `bcftools` and `tabix` for the genotype windows, Ghostscript to embed the fonts of the figures.

The comparators are found through environment variables, with defaults relative to the
repository root (the folder `comparators/` at the root is not tracked).

| Variable | Default | Scripts |
|---|---|---|
| `GEMMA_BIN` | `comparators/GEMMA-0.98.5/bin/gemma` | 30, 33, 33e |
| `QTLTOOLS_BIN` | `comparators/qtltools/bin/QTLtools` | 30, 33, 33e |
| `COMPARATORS_RLIB` | `comparators/Rlib` (R library holding `mlmm` and `emma`) | 30, 33, 33e, 34b |
| `FINEMAPINF_DIR` | `comparators/fine-mapping-inf-master` | 31, 33b |
| `GEUVADIS_SDRF` | `data/geuvadis/E-GEUV-1.sdrf.txt` | 33d |
| `ART` | `calibration/output` | table and figure scripts |
| `ART_TABLES` | `calibration/output/tables` | 35 |

## Data

The real-data analyses use public data, placed in `data/geuvadis/` (not tracked).

- GEUVADIS expression (`GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`), the list of the best
  European cis-eQTLs (`EUR373.gene.cis.FDR5.best.rs137.txt.gz`) and the sample table
  (`E-GEUV-1.sdrf.txt`), from ArrayExpress E-GEUV-1.
- Genotypes of the 1000 Genomes phase 3 release. The European windows of the eight genes
  ($n=358$, `data/geuvadis/loci/`) are built by `27e_fetch_eur_windows.R`, the
  Yoruba windows (`data/geuvadis/loci_yri/`) by `27a_prepare_yri_windows.R` from the GEUVADIS
  genotype release (chromosomes 1 to 7) and by `27b_fetch_yri_1000g.sh` from 1000 Genomes
  (chromosomes 12, 17 and 19), and the sample covariates by `33d_geuvadis_covariates.R`. The
  sample lists of `calibration/samples/` go to `data/geuvadis/`.

## Tables and figures

`bash calibration/make_tables_and_figures.sh` writes every table and figure from
`calibration/results/`. The
tables of the selection study and of the GEUVADIS application also read the datasets exported by
`29b_export_purenull.R` and `33a_export_geuvadis.R`. These contain individual genotypes and are not
distributed, and the two scripts rebuild them from the data above.

| Article | File in `calibration/output/` | Written by | From the results of |
|---|---|---|---|
| Table 1 | (in the manuscript) | | |
| Table 2 | `tables/tab_corrections_main.tex` | `37` | `36`, `36b`, `49`, `50`, `52`, `56` |
| Table 3 | `tables/tab_selection_main.tex` | `35` | `14`, `14b`, `18`, `30`, `30b`, `31`, `34`, `34b`, `34c` |
| Tables 4, 5 | `tables/tab_geu_selection_main.tex`, `tables/tab_geu_examples.tex` | `35b` | `27c`, `27d`, `33a`, `33b`, `33e`, `33f`, `33g`, `33i`, `33k`, `33l` |
| Figure 1 | `figures/fig_calibration_nu.pdf` | `37` | `01`, `36`, `36b`, `52` |
| Table S1 | `tables/tabA_mono.tex` | `24b`, `35` | `24` |
| Tables S2, S3 | `tables/tabA_laplace.tex`, `tables/tabA_firstorder.tex` | `26b`, `35` | `26` |
| Table S4 | `tables/tabS_mc_error.tex` | `39b` | `39` |
| Table S5 | `tables/tabS_calib_enlarged.tex` | `37` | `36`, `36b`, `49`, `50`, `52`, `56`, `60` |
| Table S6 | `tables/tabA_calib_power.tex` | `35`, `35f` | `06`, `06b`, `08`, `22`, `32` |
| Table S7 | `tables/tabA_persnp.tex`, `tables/tabA_nu_n.tex` | `05`, `21`, `35` | `01`, `20` |
| Table S8 | `tables/tabA_variance.tex` | `05`, `35` | `04` |
| Table S9 | `tables/tabA_clr_sizepower.tex`, `tables/tabA_clr_nu.tex` | `11`, `35`, `35f` | `10` |
| Table S10 | `tables/tabA_ge.tex` | `13`, `35` | `10`, `12` |
| Table S11 | `tables/tabA_nonmono.tex` | `07`, `35` | `06`, `06b` |
| Table S12 | `tables/tabA_prior_compare.tex`, `tables/tabA_prior_u.tex` | `09`, `35`, `35f` | `06b`, `08` |
| Table S13 | `tables/tabA_misspecified.tex` | `42b` | `42` |
| Table S14 | `tables/tabS_calib_h2range.tex` | `46b` | `36`, `36b`, `46` |
| Table S15 | `tables/tabS_nonnormal.tex` | `57b` | `57` |
| Table S16 | `tables/tab_selection.tex` | `35` | `14`, `14b`, `18`, `30`, `30b`, `31`, `34`, `34b`, `34c` |
| Table S17 | `tables/tabS_purenull.tex` | `34d` | `14c`, `29b`, `30`, `31`, `34` |
| Table S18 | `tables/tabA_sel_delta.tex` | `15`, `35` | `14`, `16`, `17` |
| Tables S19, S20 | `tables/tabA_cond_frontier.tex`, `tables/tabA_cond_delta.tex` | `19`, `35` | `14`, `18` |
| Tables S21, S22 | `tables/tab_geu_selection.tex`, `tables/tab_geu_gemma_pvalues.tex` | `35b` | `27c`, `27d`, `33a`, `33b`, `33e`, `33f`, `33g`, `33i`, `33k`, `33l` |
| Table S23 | `tables/tabS_decisions.tex` | `35d` | `33e`, `33i`, `33l` |
| Table S24 | `tables/tabS_protected_check.tex` | `35c` | `27d`, `33g`, `33l`, `41` |
| Table S25 | (in the manuscript) | | |
| Figure S1 | `figures/fig_protected_profiles.pdf` | `35e` | `27d`, `33g`, `33l` |
| Illustration of the index and the tails (not in the article) | `figures/fig_nu_tails.pdf` | `53` | none (stylized spectrum and local approximation) |

On 29 October 2026, with the exported datasets in place, the scripts rewrote from
`calibration/results/` the 30 table files of the article byte for byte and its two figures pixel for pixel.

## Rerunning the analyses

The analyses run in the following order. Most simulation scripts take the number of cores as first
argument and save one file per configuration, so that an interrupted run resumes where it stopped.

1. **Exact distributions and their checks.** `01`-`04` (exact size of the plug-in Wald and score
   tests, variance of the Wald statistic, $\nu_\delta$ of every variant), `06`, `06b` and `08` (treatments of the
   ratio, marginal Bayes factor with a half-Cauchy prior), `10` and `12` (exact conditional tests),
   `20` ($\nu_j$ against $n$), `22` (Edgeworth corrections), `24`-`26` (numerical checks of the
   identities of the supplement).
2. **Calibration study.** `36` (174 configurations, size and power of every test), `36b` (the $t$
   test at the null estimate, as EMMAX), `39` (Monte Carlo error), `42` (misspecified covariance), `43`
   (index at the estimates, boundary mass), `46` (heritability 0 and 0.5), `32` (exact
   Kenward-Roger and Satterthwaite tests, checked against MM4LMM), `22` (Edgeworth thresholds),
   `49` and `50` (plug-in tests, thresholds at the estimated ratio evaluated at the true ratio with
   independent draws), `52` (protected critical values, computed once per design and
   evaluated with independent draws) and `52b` (position of the protected critical values of the
   likelihood ratio relative to the nominal quantiles) `55` (direct simulation of the Wald test
   at level $10^{-4}$ on two configurations), `56` (restricted likelihood ratio on a common residual
   space and size-adjusted powers), `57` (direct simulation with Student and gamma errors) `60`
   (modified signed likelihood root of Skovgaard) `61` (importance sampling of the error contrasts on
   four configurations) and `62` (protected Wald sizes of two designs with the mixture sampler). `52` ran for 10 hours on 30 processes of an Intel Xeon
   E5-2699 v3 server. The same seeds give other draws on arm64 macOS, since LAPACK returns
   eigenvectors with other signs, and the results then differ within Monte Carlo error. `47`, `48`
   and `51` are earlier versions of these computations, whose protected rows `52` replaced.
3. **Selection study.** `14`, `14b`, `14c`, `16`, `17`, `18` (stepwise paths on four GEUVADIS
   windows), `29` and `29b` (export of the datasets), `30`, `30b` and `31` (published software),
   `34`, `34b` and `34c` (common metrics), `44` (bootstrap of the recalls).
4. **GEUVADIS application.** `27a`, `27b`, `27e` and `33d` (data), `27` and `27c` (paths and exact
   $p$-values), `33`, `33a`, `33b`, `33c`, `33e` and `33f` (published software), `27d`, `33g`-`33l`,
   `40`, `41` and `45` (numerical control of the protected $p$-values), `54` (candidate tests of the
   likelihood-ratio paths with a nominal $p$-value between one tenth of the level and the level), `58` (paths
   after a rank-based inverse normal transformation of the expression), `59` (dependence of the two
   likelihood ratios on the upper end of the domain of the ratio).
5. **Tables and figures.** `make_tables_and_figures.sh`.

## Citation

Please cite the article (see [`CITATION.cff`](CITATION.cff)):

```
@article{Doulabe2026Calibration,
  author  = {Doulabe, Kossi N. and Lakhal-Chaieb, Lajmi},
  title   = {Finite-sample calibration of regional mixed-model association
             tests with an estimated variance ratio},
  note    = {Preprint},
  year    = {2026}
}
```

## License

MIT, see [`LICENSE`](LICENSE). The patch `calibration/comparators/gemma_arm64_debug.patch`
modifies GEMMA and follows its GPL-3 licence.
