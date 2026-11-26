# ==============================================================================
# 05_make_note_tables.R
# LaTeX tables of the proposal note, generated from the pilot CSV files.
# Usage (from the repository root): Rscript calibration/05_make_note_tables.R
# Output: calibration/note/tab_*.tex
# ==============================================================================
R <- "calibration/results"; N <- "calibration/note"
dir.create(N, showWarnings = FALSE)

fm <- function(x) {                       # ratio formatting, math mode
  vapply(x, function(v) {
    if (is.na(v)) return("--")
    if (v >= 1000 || v < 0.01) {
      e <- floor(log10(v)); m <- v / 10^e
      return(sprintf("$%.1f\\times10^{%d}$", m, e))
    }
    if (v >= 100) return(sprintf("$%.0f$", v))
    if (v >= 10)  return(sprintf("$%.0f$", v))
    sprintf("$%.2f$", v)
  }, "")
}
fn <- function(x, d = 0) ifelse(is.na(x), "--", sprintf(paste0("$%.", d, "f$"), x))
wr <- function(file, head, body, cols) {
  writeLines(c(sprintf("\\begin{tabular}{%s}", cols), "\\toprule", head, "\\midrule",
               body, "\\bottomrule", "\\end{tabular}"), file.path(N, file))
}

syn <- read.csv(file.path(R, "exact_size_synthetic.csv"))
geu <- read.csv(file.path(R, "exact_size_geuvadis.csv"))
sco <- read.csv(file.path(R, "score_null_exact.csv"))
plg <- read.csv(file.path(R, "plugin_exact_calibration.csv"))
per <- read.csv(file.path(R, "nu_delta_per_snp_geuvadis.csv"))
vch <- if (file.exists(file.path(R, "variance_check.csv")))
  read.csv(file.path(R, "variance_check.csv")) else NULL

## ---- Table: synthetic designs --------------------------------------------------
ld  <- syn[syn$candidate == "in_LD", ]
unl <- syn[syn$candidate == "unlinked", ]
key <- function(d) paste(d$n, d$blocks, d$h2)
sc  <- sco[sco$design == "AR1_rho0.95", ]
ms  <- match(key(ld), key(sc)); mu <- match(key(ld), key(unl))
body <- sprintf("%d & %d & %.2f & %s & %s & %s & %s & %s & %s & %s \\\\",
  ld$n, ld$blocks, ld$h2, fn(ld$nu_delta), fm(ld$ratio_1e.03), fm(ld$ratio_5e.08),
  fm(ld$satt_ratio_5e.08), fm(sc$score_ratio_1e.03[ms]), fm(sc$score_ratio_5e.08[ms]),
  fm(unl$ratio_5e.08[mu]))
wr("tab_synthetic.tex",
   paste("$n$ & blocs & $h^2$ & $\\nu_j$ & Wald $10^{-3}$ & Wald $5\\cdot10^{-8}$ &",
         "Satt. $5\\cdot10^{-8}$ & Score $10^{-3}$ & Score $5\\cdot10^{-8}$ &",
         "Non li\\'e $5\\cdot10^{-8}$ \\\\"), body, "rrrrrrrrrr")

## ---- Table: GEUVADIS leads ---------------------------------------------------------
sg <- sco[sco$design == "GEUVADIS_cis", ]
sg$gene <- sub("^lead_", "", sg$candidate)
mg <- match(paste(geu$gene, geu$h2), paste(sg$gene, sg$h2))
body <- sprintf("%s & %.2f & %s & %s & %s & %s & %s & %s & %s \\\\",
  geu$gene, geu$h2, fn(geu$nu_delta, 1), fm(geu$ratio_5e.02), fm(geu$ratio_1e.03),
  fm(geu$ratio_5e.08), fm(geu$satt_ratio_5e.08), fm(sg$score_ratio_1e.03[mg]),
  fm(sg$score_ratio_5e.08[mg]))
wr("tab_geuvadis.tex",
   paste("G\\`ene & $h^2$ & $\\nu_j$ & Wald $0.05$ & Wald $10^{-3}$ & Wald $5\\cdot10^{-8}$ &",
         "Satt. $5\\cdot10^{-8}$ & Score $10^{-3}$ & Score $5\\cdot10^{-8}$ \\\\"),
   body, "lrrrrrrrr")

## ---- Table: per-SNP effective df (wide: h2 = 0.05 | 0.20) ---------------------------------
agg <- do.call(rbind, lapply(split(per, list(per$gene, per$h2)), function(z) data.frame(
  gene = z$gene[1], h2 = z$h2[1], m = nrow(z), med = median(z$nu_delta),
  f50 = 100 * mean(z$nu_delta < 50), f100 = 100 * mean(z$nu_delta < 100),
  lead = z$nu_delta[z$is_lead][1])))
a1 <- agg[agg$h2 == 0.05, ]; a2 <- agg[agg$h2 == 0.20, ]
a1 <- a1[order(a1$gene), ]; a2 <- a2[match(a1$gene, a2$gene), ]
body <- sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s & %s \\\\", a1$gene, a1$m,
  fn(a1$med, 1), fn(a1$f50), fn(a1$f100), fn(a1$lead, 1),
  fn(a2$med, 1), fn(a2$f50), fn(a2$f100), fn(a2$lead, 1))
wr("tab_persnp.tex",
   paste("& & \\multicolumn{4}{c}{$h^2=0.05$} & \\multicolumn{4}{c}{$h^2=0.20$} \\\\",
         "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
         "G\\`ene & SNP & m\\'ediane & \\% $<50$ & \\% $<100$ & lead &",
         "m\\'ediane & \\% $<50$ & \\% $<100$ & lead \\\\"), body, "lrrrrrrrrr")

## ---- Table: plug-in exact law ----------------------------------------------------------
lab <- ifelse(plg$design == "AR1_rho0.95", "AR(1), $n=1000$, 5 blocs",
              paste0("GEUVADIS ", sub("^lead_", "", plg$candidate)))
body <- sprintf("%s & %.2f & %s & %s & %s & %s & %s & %s & %s & %s \\\\", lab, plg$h2,
  fm(plg$naive_ratio_1e.03), fm(plg$plugin_exact_ratio_1e.03),
  fm(plg$naive_ratio_1e.05), fm(plg$plugin_exact_ratio_1e.05),
  fm(plg$naive_ratio_5e.08), fm(plg$plugin_exact_ratio_5e.08),
  fn(plg$c_naive_5e.08, 2), fn(plg$c_exact_at_delta0_5e.08, 2))
wr("tab_plugin.tex",
   paste("Configuration & $h^2$ & \\multicolumn{2}{c}{$10^{-3}$} &",
         "\\multicolumn{2}{c}{$10^{-5}$} & \\multicolumn{2}{c}{$5\\cdot10^{-8}$} &",
         "$c$ na\\\"if & $c$ exact \\\\",
         "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
         "& & na\\\"if & plug-in & na\\\"if & plug-in & na\\\"if & plug-in & & \\\\"),
   body, "lrrrrrrrrr")

## ---- Table: variance check -----------------------------------------------------------------
if (!is.null(vch)) {
  vch$lab <- ifelse(vch$design == "AR1_rho0.95",
                    sprintf("AR(1), $n=%d$, %d blocs", vch$n, vch$blocks),
                    paste0("GEUVADIS ", sub("^lead_", "", vch$candidate)))
  src <- rbind(transform(syn[syn$candidate == "in_LD", c("n", "blocks", "h2", "ratio_5e.08")],
                         k = paste("AR", n, blocks, h2)),
               transform(geu[, c("n", "h2", "ratio_5e.08")], blocks = NA,
                         k = paste("GE", geu$gene, h2))[, c("n", "blocks", "h2", "ratio_5e.08", "k")])
  vk <- ifelse(vch$design == "AR1_rho0.95", paste("AR", vch$n, vch$blocks, vch$h2),
               paste("GE", sub("^lead_", "", vch$candidate), vch$h2))
  rt <- src$ratio_5e.08[match(vk, src$k)]
  keep <- vch$h2 == 0.05
  v <- vch[keep, ]
  body <- sprintf("%s & %s & %s & %s & %s \\\\", v$lab, fn(v$nu_delta, 1),
    fn(v$var_first_order_rel, 3), fn(v$var_exact_rel, 3), fm(rt[keep]))
  wr("tab_variance.tex",
     paste("Configuration ($h^2=0.05$) & $\\nu_j$ & variance, 1\\textsuperscript{er} ordre &",
           "variance exacte & queue $5\\cdot10^{-8}$ \\\\"), body, "lrrrr")
}
cat("tables written to", N, "\n")
