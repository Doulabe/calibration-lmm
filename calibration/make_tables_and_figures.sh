#!/bin/bash
# Writes the tables and figures of the article from the per-configuration values in
# calibration/results: intermediate tables in calibration/note, then the tables and figures of
# the article in calibration/output (another folder with ART=<folder>). 38 and 38b print the
# numbers quoted in the calibration section. The tables of the selection study and of the GEUVADIS
# application also read the datasets exported by 29b and 33a (calibration/results/*/data), which
# are not distributed and are rebuilt from the GEUVADIS data (see README.md).
# Usage (from the repository root): bash calibration/make_tables_and_figures.sh
set -e
cd "$(dirname "$0")/.."
export ART="${ART:-calibration/output}"
export ART_TABLES="${ART_TABLES:-$ART/tables}"
mkdir -p calibration/note "$ART/tables" "$ART/figures"
for s in 05_make_note_tables 07_make_extension_tables 09_make_prior_tables 11_make_clr_tables \
         13_make_ge_tables 15_make_selection_tables 19_make_conditional_tables 21_make_nu_table \
         23_make_edgeworth_tables 24b_make_mono_table 26b_make_symmetric_tables 28_make_geuvadis_tables \
         34d_purenull_summary 35_make_article_tables 35b_geuvadis_cov_tables 35c_protected_check_table \
         35d_gemma_decisions_controlled 35e_protected_profiles_figure 35f_merge_supp_tables \
         37_make_calibration_figure 38_calibration_numbers 38b_within_window_spread 39b_mc_error_table \
         42b_misspecified_table 46b_h2range_table 57b_nonnormal_table; do
  echo "== $s"
  Rscript "calibration/$s.R"
done
