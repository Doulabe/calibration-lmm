# ==============================================================================
# 30_comparators_selection.R
# Published software as comparators on the 1,800 datasets of the selection-level
# study (exported by 29_export_selection_datasets.R; four GEUVADIS cis windows,
# n = 358, regional background h2_bg in {0.10, 0.25}, K in {0, 1, 3} causals).
#   gemma_wald / gemma_lrt / gemma_score
#       GEMMA 0.98.5 (-lmm 4), regional kernel K = ZZ'/m held fixed, forward
#       conditional paths: the top variant of each test enters as a covariate.
#       Wald: delta re-estimated by REML under H1 for each variant; LRT: ML under
#       both models; score: delta estimated once under the null of the step.
#   mlmm     MLMM 0.1.1 (Segura et al. 2012, EMMA 1.1.2), mlmm(y, X, K, nbchunks = 2,
#            maxsteps = 11); native extBIC and mBonf models plus the forward path.
#   bslmm    GEMMA -bslmm 1 (Zhou et al. 2013) on the window genotypes, 10,000 burn-in
#            and 40,000 sampling iterations (chain-length control in 32_*); inclusion probabilities.
#   farmcpu  rMVP 1.4.6 MVP.FarmCPU (Liu et al. 2016), static bins of 10, 50 and 100 kb.
#   qtltools QTLtools 1.3.1 cis pass: 1,000 permutations, beta-approximated per-gene
#            threshold (as qtltools_runFDR_cis.R) at gene-level alpha 0.01, 0.05, 0.20,
#            then the conditional forward-backward --mapping pass.
# Every path is recorded down to a liberal level (p <= 1e-3), so the selection at any
# stricter per-test threshold is a prefix of it.
# Output: results/comparators/r_arms/<id>.rds
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/30_comparators_selection.R [cores=3] [--test]
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
args <- commandArgs(TRUE)
NCORES <- if (length(args) && grepl("^[0-9]+$", args[1])) as.integer(args[1]) else 3L
TEST <- "--test" %in% args
ARMS <- strsplit(Sys.getenv("ARMS", "gemma,mlmm,farmcpu,qtltools"), ",")[[1]]   # bslmm runs as a separate pass
ONLY_MAIN <- Sys.getenv("ONLY_MAIN", "0") == "1"                               # skip the 800 extra global nulls
suppressPackageStartupMessages({ library(parallel); library(bigmemory); library(rMVP) })
.libPaths(c(Sys.getenv("COMPARATORS_RLIB", "comparators/Rlib"), .libPaths()))
suppressPackageStartupMessages(library(mlmm))
GEMMA <- normalizePath(Sys.getenv("GEMMA_BIN", "comparators/GEMMA-0.98.5/bin/gemma"))
QTLT <- normalizePath(Sys.getenv("QTLTOOLS_BIN", "comparators/qtltools/bin/QTLtools"))
GENES <- data.frame(symbol = c("ERAP2", "RPS26", "SLFN5", "PEX6-region"), chr = c(5L, 12L, 17L, 6L),
                    tss = c(96211643L, 56435637L, 33570086L, 42946958L),
                    window = c(499538L, 499777L, 499908L, 499951L), stringsAsFactors = FALSE)
P_LIB <- 1e-3; K_MAX <- 10L
BS_W <- Sys.getenv("BS_W", "10000"); BS_S <- Sys.getenv("BS_S", "40000")   # BSLMM burn-in, sampling
ROOT <- file.path("calibration/results", Sys.getenv("COMP_ROOT", "comparators"))   # comparators_purenull: y = noise (29b)
OUT <- file.path(ROOT, if (identical(ARMS, "bslmm")) "r_arms_bslmm" else "r_arms"); WORK <- file.path(ROOT, "work"); TMP <- file.path(ROOT, "tmp")
for (d in c(OUT, WORK, TMP)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ---- per-window inputs, written once -------------------------------------------------------
D <- setNames(lapply(GENES$symbol, function(s) readRDS(file.path(ROOT, "data", paste0(s, ".rds")))), GENES$symbol)
smp <- function(n) sprintf("i%03d", seq_len(n))
for (s in GENES$symbol) {
  w <- file.path(WORK, s); dir.create(w, showWarnings = FALSE)
  d <- D[[s]]; n <- nrow(d$Z); m <- ncol(d$Z); ids <- paste0("s", seq_len(m))
  if (!file.exists(file.path(w, "geno.txt")))
    data.table::fwrite(data.frame(ids, "A", "B", t(d$G)), file.path(w, "geno.txt"), col.names = FALSE, sep = ",")
  if (!file.exists(file.path(w, "K.txt")))
    data.table::fwrite(as.data.frame(tcrossprod(d$Z) / m), file.path(w, "K.txt"), col.names = FALSE, sep = "\t")
  if (!file.exists(file.path(w, "genotypes.bed.gz.tbi"))) {
    g <- GENES[GENES$symbol == s, ]
    bed <- data.frame(g$chr, d$pos - 1L, d$pos, ids, ids, "+", t(d$G), check.names = FALSE)
    colnames(bed) <- c("#chr", "start", "end", "id", "gid", "strand", smp(n))
    bed <- bed[order(d$pos), ]
    data.table::fwrite(bed, file.path(w, "genotypes.bed"), sep = "\t", quote = FALSE)
    stopifnot(system2("bgzip", c("-f", shQuote(file.path(w, "genotypes.bed")))) == 0,
              system2("tabix", c("-f", "-p", "bed", shQuote(file.path(w, "genotypes.bed.gz")))) == 0)
  }
}

## ---- GEMMA ---------------------------------------------------------------------------------------
gemma_lmm4 <- function(gene, y, G_sel, tmp, tag) {
  w <- normalizePath(file.path(WORK, gene))
  write.table(y, file.path(tmp, "y.txt"), row.names = FALSE, col.names = FALSE)
  write.table(cbind(1, G_sel), file.path(tmp, "cov.txt"), row.names = FALSE, col.names = FALSE)
  st <- system2(GEMMA, c("-g", file.path(w, "geno.txt"), "-p", file.path(tmp, "y.txt"), "-k", file.path(w, "K.txt"),
                         "-c", file.path(tmp, "cov.txt"), "-lmm", "4", "-notsnp", "-outdir", tmp, "-o", tag),
                stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("GEMMA failed (", tag, ")")
  A <- data.table::fread(file.path(tmp, paste0(tag, ".assoc.txt")), data.table = FALSE)
  data.frame(idx = as.integer(sub("^s", "", A$rs)), p_wald = A$p_wald, p_lrt = A$p_lrt, p_score = A$p_score,
             l_remle = A$l_remle, l_mle = A$l_mle)
}
gemma_paths <- function(gene, y, tmp) {
  G <- D[[gene]]$G; cache <- list()
  run <- function(sel) {
    key <- paste0("k:", paste(sort(sel), collapse = ","))
    if (is.null(cache[[key]])) cache[[key]] <<- gemma_lmm4(gene, y, G[, sel, drop = FALSE], tmp, paste0("r", length(cache)))
    cache[[key]]
  }
  out <- list()
  for (test in c("wald", "lrt", "score")) {
    sel <- integer(0); pv <- numeric(0)
    for (k in seq_len(K_MAX)) {
      A <- run(sel); p <- A[[paste0("p_", test)]]; b <- which.min(p)
      if (!length(b) || p[b] > P_LIB) break
      sel <- c(sel, A$idx[b]); pv <- c(pv, p[b])
    }
    out[[test]] <- data.frame(idx = sel, p = pv)
  }
  out$n_runs <- length(cache)
  out
}

## ---- BSLMM ---------------------------------------------------------------------------------------
bslmm_pip <- function(gene, y, tmp, seed) {
  w <- normalizePath(file.path(WORK, gene))
  write.table(y, file.path(tmp, "yb.txt"), row.names = FALSE, col.names = FALSE)
  st <- system2(GEMMA, c("-g", file.path(w, "geno.txt"), "-p", file.path(tmp, "yb.txt"), "-bslmm", "1", "-notsnp",
                         "-w", BS_W, "-s", BS_S, "-seed", seed, "-outdir", tmp, "-o", "bs"),
                stdout = FALSE, stderr = FALSE)
  if (st != 0) stop("BSLMM failed")
  P <- data.table::fread(file.path(tmp, "bs.param.txt"), data.table = FALSE)
  hl <- strsplit(trimws(readLines(file.path(tmp, "bs.hyp.txt"))), "[[:space:]]+")   # trailing tab: parse by hand
  H <- as.data.frame(do.call(rbind, lapply(hl[-1], as.numeric))); names(H) <- hl[[1]]
  keep <- P$gamma > 0.01
  list(idx = as.integer(sub("^s", "", P$rs[keep])), pip = P$gamma[keep],
       pve = mean(H$pve), pge = mean(H$pge), n_gamma = mean(H$n_gamma))
}

## ---- MLMM ----------------------------------------------------------------------------------------
mlmm_run <- function(gene, y) {
  d <- D[[gene]]; n <- nrow(d$Z); m <- ncol(d$Z)
  X <- d$Z; colnames(X) <- paste0("s", seq_len(m)); rownames(X) <- smp(n)
  K <- tcrossprod(d$Z) / m; dimnames(K) <- list(smp(n), smp(n)); yy <- setNames(y, smp(n))
  sink(tempfile()); on.exit(sink())
  fit <- mlmm(yy, X, K, nbchunks = 2, maxsteps = 11)
  cof <- function(v) { v <- v[!is.na(v)]; as.integer(sub("^s", "", v)) }
  fw <- fit$step_table[grepl("^fwd", fit$step_table$step_) & fit$step_table$ncof > 0, ]
  entry_p <- vapply(seq_len(nrow(fw)), function(k) {                 # min p at the step before entry
    ps <- fit$pval_step[[k]]$out; min(ps$pval, na.rm = TRUE) }, 0)
  list(extBIC = cof(fit$opt_extBIC$cof), mbonf = cof(fit$opt_mbonf$cof),
       fwd = data.frame(idx = as.integer(sub("^\\+s", "", fw$cof)), p = entry_p, h2 = fw$h2),
       h2_null = fit$step_table$h2[1])
}

## ---- FarmCPU -------------------------------------------------------------------------------------
farmcpu_run <- function(gene, y) {
  d <- D[[gene]]; m <- ncol(d$G); n <- nrow(d$G)
  geno <- bigmemory::as.big.matrix(t(d$G))
  map <- data.frame(SNP = paste0("s", seq_len(m)), CHROM = GENES$chr[GENES$symbol == gene], POS = d$pos)
  fc <- MVP.FarmCPU(phe = data.frame(id = smp(n), y = y), geno = geno, map = map, CV = NULL, maxLoop = 10,
                    method.bin = "static", bin.size = c(1e4, 5e4, 1e5), bin.selection = seq(5, 50, 5),
                    p.threshold = 0.05 / m, ncpus = 1, verbose = FALSE)
  p <- fc[, ncol(fc)]; keep <- which(is.finite(p) & p <= P_LIB)
  data.frame(idx = keep, p = p[keep])
}

## ---- QTLtools -----------------------------------------------------------------------------------
qtltools_run <- function(gene, y, tmp, seed, alphas = c(0.01, 0.05, 0.20)) {
  g <- GENES[GENES$symbol == gene, ]; w <- normalizePath(file.path(WORK, gene)); n <- length(y)
  ph <- data.frame(g$chr, g$tss - 1L, g$tss, gene, gene, "+", t(y), check.names = FALSE)
  colnames(ph) <- c("#chr", "start", "end", "id", "gid", "strand", smp(n))
  fb <- file.path(tmp, "expression.bed")
  data.table::fwrite(ph, fb, sep = "\t", quote = FALSE)
  stopifnot(system2("bgzip", c("-f", shQuote(fb))) == 0, system2("tabix", c("-f", "-p", "bed", shQuote(paste0(fb, ".gz")))) == 0)
  common <- c("cis", "--vcf", file.path(w, "genotypes.bed.gz"), "--bed", paste0(fb, ".gz"), "--window", g$window,
              "--seed", seed)
  per <- file.path(tmp, "perm.txt")
  stopifnot(system2(QTLT, c(common, "--permute", "1000", "--out", per), stdout = FALSE, stderr = FALSE) == 0)
  p <- read.table(per)
  out <- list(perm = p)
  for (a in alphas) {
    pb <- qbeta(a, p$V14, p$V15); tt <- qf(pb, 1, p$V13, lower.tail = FALSE)
    corr2 <- tt / (p$V13 + tt); cutoff <- pf(p$V12 * corr2 / (1 - corr2), 1, p$V12, lower.tail = FALSE)
    thr <- file.path(tmp, "thr.txt"); write.table(data.frame(gene, cutoff), thr, row.names = FALSE, col.names = FALSE, quote = FALSE)
    mp <- file.path(tmp, "map.txt")
    stopifnot(system2(QTLT, c(common, "--mapping", thr, "--out", mp), stdout = FALSE, stderr = FALSE) == 0)
    tab <- if (file.exists(mp) && file.info(mp)$size > 0) read.table(mp) else data.frame()
    hits <- if (nrow(tab)) tab[tab$V21 == 1, , drop = FALSE] else tab
    out[[sprintf("a%.2f", a)]] <- if (nrow(hits)) as.integer(sub("^s", "", hits$V8)) else integer(0)
  }
  out
}

## ---- one dataset ---------------------------------------------------------------------------------
JOBS <- do.call(rbind, lapply(GENES$symbol, function(s) data.frame(gene = s, k = seq_len(ncol(D[[s]]$Y)),
                                                                     id = D[[s]]$meta$id, stringsAsFactors = FALSE)))
one_rep <- function(j) {
  gene <- JOBS$gene[j]; k <- JOBS$k[j]; id <- JOBS$id[j]
  ck <- file.path(OUT, paste0(id, ".rds"))
  if (!TEST && file.exists(ck) && file.size(ck) > 0) return("cached")
  y <- D[[gene]]$Y[, k]; tmp <- file.path(TMP, id); dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
  seed <- 4000L + j; tm <- c()
  arm <- function(nm, expr) { if (!nm %in% ARMS) return(NULL)
    t0 <- proc.time()[[3]]; v <- tryCatch(expr, error = function(e) conditionMessage(e))
    tm[nm] <<- proc.time()[[3]] - t0; v }
  gp <- arm("gemma", gemma_paths(gene, y, tmp)); ml <- arm("mlmm", mlmm_run(gene, y))
  bs <- arm("bslmm", bslmm_pip(gene, y, tmp, seed)); fc <- arm("farmcpu", farmcpu_run(gene, y))
  qt <- arm("qtltools", qtltools_run(gene, y, tmp, seed))
  res <- list(id = id, gene = gene, k = k, gemma = gp, mlmm = ml, bslmm = bs, farmcpu = fc, qtltools = qt, secs = tm)
  unlink(tmp, recursive = TRUE)
  if (!TEST) saveRDS(res, ck)
  if (TEST) res else "done"
}

if (TEST) {
  for (j in c(which(JOBS$id == "ERAP2_bg0.25_K3_c0.08_b001"), which(JOBS$id == "ERAP2_bg0.25_K0_c0.00_b001"))) {
    r <- one_rep(j); cat("\n", r$id, " causal:", D[[r$gene]]$meta$causal[r$k], "\n"); print(round(r$secs, 1))
    str(r[c("gemma", "mlmm", "farmcpu", "qtltools")], max.level = 2, vec.len = 6)
    cat("BSLMM top PIP:\n"); if (is.list(r$bslmm)) print(head(data.frame(idx = r$bslmm$idx, pip = r$bslmm$pip)[order(-r$bslmm$pip), ], 5))
  }
} else {
  ord <- order(JOBS$k, JOBS$gene)                       # interleave windows
  if (ONLY_MAIN) ord <- ord[!grepl("_b(0[2-9][6-9]|0[3-9][0-9]|1[0-9][0-9])$", JOBS$id[ord]) | grepl("_K[13]_", JOBS$id[ord])]
  st <- unlist(mclapply(ord, function(j) tryCatch(one_rep(j), error = function(e)
    paste("FAIL", JOBS$id[j], conditionMessage(e))), mc.cores = NCORES, mc.preschedule = FALSE))
  cat("datasets:", sum(st == "done"), "new,", sum(st == "cached"), "cached,", sum(grepl("FAIL", st)), "failed\n")
  print(head(st[grepl("FAIL", st)]))
  writeLines(format(Sys.time()), file.path(ROOT, "R_ARMS_DONE"))
}
