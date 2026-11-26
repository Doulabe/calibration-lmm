# ==============================================================================
# 41_replicate_decisive.R
# Two numerical checks of the refined protected p-values.
#  rep    independent replicate, with fresh random seeds and otherwise the same procedure, of
#         the decisive cases of Supplementary Table S25 and of the near-threshold decisions of
#         GEMMA's paths: the Wald and Kenward-Roger cases with 100,000 and 200,000 draws as in
#         33k/33l, the score and likelihood-ratio cases of GEMMA's paths with 4,000 and 10,000
#         draws as in 33g/33i (GEMMA's range for delta_hat), and the two decisive selections of
#         the likelihood-ratio path of this study on the domain of its decisions, [0, 1e3]
#         (27d --domain=1e3), run by the mode "paths" in a separate process because 27d and 33g
#         define functions with the same names;
#  caps   protected p-value of the decisive maximum-likelihood cases (the three score stops of
#         the Yoruba sample attributed to the reference distribution and the undetermined
#         second selection of GEMMA's likelihood-ratio path in the PEX6 region) when the upper
#         end of the domain of delta_hat is 1e4, 1e6 or 1e8 instead of GEMMA's 1e5, the
#         observed statistic being recomputed on the same domain (grid of 0.1 decade from 1e-5).
# In every case the supremum over the true ratio delta0 is searched on [0, 1e10] as before.
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/41_replicate_decisive.R gemma [cores=2]
#   nice -n 10 Rscript calibration/41_replicate_decisive.R paths
# Output: results/replicate_check/<id>.rds and results/ml_cap_protected/<id>.rds (format of 33g)
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
a0 <- commandArgs(TRUE); MODE <- a0[1]; stopifnot(MODE %in% c("gemma", "paths"))
NC <- if (length(a0) >= 2 && grepl("^[0-9]+$", a0[2])) as.integer(a0[2]) else 1L
OUT_REP <- "calibration/results/replicate_check"; OUT_CAP <- "calibration/results/ml_cap_protected"
dir.create(OUT_REP, showWarnings = FALSE, recursive = TRUE); dir.create(OUT_CAP, showWarnings = FALSE, recursive = TRUE)
grab <- function(file, from, to) {                       # lines of a script between two anchors (to excluded)
  s <- readLines(file); i <- grep(from, s)[1]; j <- if (is.null(to)) length(s) + 1L else grep(to, s)[1]
  stopifnot(!is.na(i), !is.na(j), j > i); s[i:(j - 1L)] }

if (MODE == "paths") {
  ## ---- likelihood-ratio path of this study, decision domain [0, 1e3] (27d --domain=1e3) ------------
  s27 <- readLines("calibration/27d_protected_check_paths.R")
  eval(parse(text = s27[seq_len(grep("^run_case <- ", s27)[1] - 1L)]))
  GRID2 <- c(0, 10^seq(-4, 3, length.out = 41))              # as 27d with --domain=1e3
  run_path <- function(k) {
    cs <- CASES[k, ]; id <- paste0(sub("_D1e3$", "", cs$id), "_D1e3_rep"); ck <- file.path(OUT_REP, paste0(id, ".rds"))
    if (file.exists(ck)) return(paste("cached", id))
    set.seed(77100 + k); t0 <- Sys.time()                    # 27d used 3390 + k
    W <- load_window(cs$tag); X <- scale(W$G); X[!is.finite(X)] <- 0
    y <- as.numeric(scale(EXPR[[cs$gene]][W$samples])); C <- covariates(W$samples)
    jb <- match(cs$pos_before, W$pos); j <- match(cs$pos, W$pos); stopifnot(!is.na(jb), !is.na(j))
    sp <- spectral_step(X, y, C, jb, j); o <- stats_from_coords(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
    n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m)
    prof1 <- t(vapply(D1, function(d0) exact_lrt(sp, d0, o$lrt, B1), numeric(3)))
    a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
    D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
    prof2 <- t(vapply(D2, function(d0) exact_lrt(sp, d0, o$lrt, B2), numeric(3)))
    prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2)); prof$id <- id
    i <- which.max(prof$log10p)
    info <- data.frame(id = id, pop = "YRI", gene = cs$symbol, pos = cs$pos, step = 2L, stat = "lrt (path of this study)",
                       cap = 1e3, level = pchisq(pen, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                       nominal = pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                       new_max = prof$log10p[i], new_argmax = prof$delta0[i], controlled = max(prof$log10up),
                       elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    saveRDS(list(info = info, profile = prof), ck); sprintf("done %s (%.1f min)", id, info$elapsed_min)
  }
  print(unlist(lapply(seq_len(nrow(CASES)), function(k) tryCatch(run_path(k), error = function(e) paste("FAIL", k, conditionMessage(e))))))
  quit(save = "no")
}

## ---- GEMMA's statistics: functions of 33g/33i (score, likelihood ratio) and 33l (Wald, KR) ------------
Sys.setenv(PC_SOURCE_ONLY = "1")
try(source("calibration/33g_protected_check.R"), silent = TRUE)
eval(parse(text = grab("calibration/33i_protected_check_all.R", "^## ---- lean statistics", "^## ---- cases")))
eval(parse(text = grab("calibration/33i_protected_check_all.R", "^setup_case2 <- ", "^run_one <- ")))
CH <- 10000L
eval(parse(text = grab("calibration/33l_bigB_all.R", "^COMPV <- ", "^## cases:")))
stopifnot(exists("exact_any"), exists("setup_case2"), exists("big_one"), exists("setup_step"))

REP <- data.frame(                                          # (pop, gene, pos, step, stat) of the replicated cases
  pop  = c("EUR", "EUR", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI"),
  gene = c("PEX6-region", "ERAP2", "ZNF266", "PEX6-region", "TRA2A-AS", "ZNF266", "PEX6-region", "PEX6-region"),
  pos  = c(42944850L, 96252432L, 9477868L, 42915920L, 23134790L, 9477868L, 43186801L, 43011766L),
  step = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L),
  stat = c("wald", "kr", "wald", "wald", "score", "score", "score", "lrt"), stringsAsFactors = FALSE)
REP$cap <- 1e5; REP$kind <- "rep"
CAPC <- REP[REP$stat %in% c("score", "lrt"), ]
CAPC <- do.call(rbind, lapply(c(1e4, 1e6, 1e8), function(cp) transform(CAPC, cap = cp, kind = "cap")))
JOBS <- rbind(REP, CAPC)
JOBS$id <- sprintf("%s_%s_%d_%d_%s_%s", JOBS$pop, gsub("[^A-Za-z0-9]", "", JOBS$gene), JOBS$pos, JOBS$step, JOBS$stat,
                   ifelse(JOBS$kind == "rep", "rep", sprintf("cap%g", JOBS$cap)))
JOBS <- JOBS[order(JOBS$pop == "EUR", decreasing = TRUE), ]  # the long EUR cases first

run_job <- function(k) {
  cs <- JOBS[k, ]; ck <- file.path(if (cs$kind == "rep") OUT_REP else OUT_CAP, paste0(cs$id, ".rds"))
  if (file.exists(ck)) return(paste("cached", cs$id))
  set.seed(77000 + k); t0 <- Sys.time()                      # 33g, 33i, 33k and 33l used other seeds
  GRID2 <<- c(0, 10^seq(-5, log10(cs$cap), by = 0.1))         # domain of delta_hat (GEMMA's range at 1e5)
  if (cs$stat %in% c("wald", "kr")) {
    S <- setup_step(cs$pop, cs$gene, cs$pos, cs$step)
    f1 <- function(d0) big_one(S$sp, d0, S$obs, cs$stat, 100000L); f2 <- function(d0) big_one(S$sp, d0, S$obs, cs$stat, 200000L)
    level <- S$level; nominal <- NA
  } else {
    S <- setup_case2(cs$pop, cs$gene, cs$pos, cs$step, cs$stat)
    f1 <- function(d0) exact_any(S, d0, cs$stat, B1); f2 <- function(d0) exact_any(S, d0, cs$stat, B2)
    level <- S$info$level; nominal <- S$info$ours_nominal
  }
  prof1 <- t(vapply(D1, f1, numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, f2, numeric(3)))
  prof <- rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2)); prof$id <- cs$id
  i <- which.max(prof$log10p)
  info <- data.frame(id = cs$id, pop = cs$pop, gene = cs$gene, pos = cs$pos, step = cs$step, stat = cs$stat, cap = cs$cap,
                     kind = cs$kind, level = level, nominal = nominal, new_max = prof$log10p[i], new_argmax = prof$delta0[i],
                     controlled = max(prof$log10up), elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  saveRDS(list(info = info, profile = prof), ck); sprintf("done %s (%.1f min)", cs$id, info$elapsed_min)
}
cat("jobs:", nrow(JOBS), "\n")
print(unlist(parallel::mclapply(seq_len(nrow(JOBS)), function(k) tryCatch(run_job(k), error = function(e) paste("FAIL", JOBS$id[k], conditionMessage(e))),
                                mc.cores = NC, mc.preschedule = FALSE)))
