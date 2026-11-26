# ==============================================================================
# 45_crn_decisive.R
# Two checks of the refined protected p-values.
#  gemma, paths  the decisive cases of 41_replicate_decisive.R (eight cases of GEMMA's paths, two
#                selections of the likelihood-ratio path of this study) recomputed with common random
#                numbers: the same draws of the error contrasts at every delta0 of a stage (the seed is
#                reset before each evaluation, the contrasts being drawn as standard normals scaled by
#                sqrt(1 + delta0 ell_i)), so that the profile p(delta0) is smooth and its maximum less
#                biased upward than with independent draws at each delta0 (33g, 33i, 33l, 27d, 41);
#  summary       for these cases and for every entry of GEMMA's paths (35d), the classification with
#                the pointwise margin of two standard errors and with a margin made simultaneous over
#                the G points of the grid by Bonferroni, z = qnorm(1 - pnorm(-2) / G).
# Usage (from the repository root):
#   nice -n 10 Rscript calibration/45_crn_decisive.R gemma [cores=2]
#   nice -n 10 Rscript calibration/45_crn_decisive.R paths
#   Rscript calibration/45_crn_decisive.R summary
# Output: results/crn_check/<id>.rds (format of 33g), results/crn_check_summary.csv,
#         results/simultaneous_decisions.csv
# ==============================================================================
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
a0 <- commandArgs(TRUE); MODE <- a0[1]; stopifnot(MODE %in% c("gemma", "paths", "summary"))
NC <- if (length(a0) >= 2 && grepl("^[0-9]+$", a0[2])) as.integer(a0[2]) else 1L
RES <- "calibration/results"; OUT_CRN <- file.path(RES, "crn_check")   # 33g and 27d define their own OUT
dir.create(OUT_CRN, showWarnings = FALSE, recursive = TRUE)
grab <- function(file, from, to) {                       # lines of a script between two anchors (to excluded)
  s <- readLines(file); i <- grep(from, s)[1]; j <- if (is.null(to)) length(s) + 1L else grep(to, s)[1]
  stopifnot(!is.na(i), !is.na(j), j > i); s[i:(j - 1L)] }
crn <- function(f, seed) function(d0) { set.seed(seed); f(d0) }   # same draws at every delta0 of a stage
profile2 <- function(f1, f2, seed) {                      # the two stages of 33g, with common random numbers
  prof1 <- t(vapply(D1, crn(f1, seed), numeric(3)))
  a <- log10(pmax(D1[which.max(prof1[, "log10p"])], 1e-3))
  D2 <- setdiff(10^seq(a - 0.3, a + 0.3, by = 0.05), D1)
  prof2 <- t(vapply(D2, crn(f2, seed + 1L), numeric(3)))
  rbind(data.frame(stage = 1L, delta0 = D1, prof1), data.frame(stage = 2L, delta0 = D2, prof2))
}

if (MODE == "paths") {
  ## ---- likelihood-ratio path of this study, decision domain [0, 1e3] (as 41, mode paths) ---------------
  s27 <- readLines("calibration/27d_protected_check_paths.R")
  eval(parse(text = s27[seq_len(grep("^run_case <- ", s27)[1] - 1L)]))
  GRID2 <- c(0, 10^seq(-4, 3, length.out = 41))
  for (k in seq_len(nrow(CASES))) {
    cs <- CASES[k, ]; id <- paste0(sub("_D1e3$", "", cs$id), "_D1e3_crn"); ck <- file.path(OUT_CRN, paste0(id, ".rds"))
    if (file.exists(ck)) { cat("cached", id, "\n"); next }
    t0 <- Sys.time()
    W <- load_window(cs$tag); X <- scale(W$G); X[!is.finite(X)] <- 0
    y <- as.numeric(scale(EXPR[[cs$gene]][W$samples])); C <- covariates(W$samples)
    sp <- spectral_step(X, y, C, match(cs$pos_before, W$pos), match(cs$pos, W$pos))
    o <- stats_from_coords(sp, matrix(sp$obs$sk), sp$obs$sp, sp$obs$ssr)
    n <- nrow(X); m <- ncol(X); pen <- log(n) + 2 * log(m)
    prof <- profile2(function(d0) exact_lrt(sp, d0, o$lrt, B1), function(d0) exact_lrt(sp, d0, o$lrt, B2), 45100L + 2L * k)
    prof$id <- id; i <- which.max(prof$log10p)
    info <- data.frame(id = id, pop = "YRI", gene = cs$symbol, pos = cs$pos, step = 2L, stat = "lrt (path of this study)",
                       cap = 1e3, level = pchisq(pen, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                       nominal = pchisq(o$lrt, 1, lower.tail = FALSE, log.p = TRUE) / log(10),
                       new_max = prof$log10p[i], new_argmax = prof$delta0[i], controlled = max(prof$log10up),
                       elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    saveRDS(list(info = info, profile = prof), ck); cat(sprintf("done %s (%.1f min)\n", id, info$elapsed_min))
  }
  quit(save = "no")
}

if (MODE == "gemma") {
  ## ---- GEMMA's statistics: functions of 33g/33i (score, likelihood ratio) and 33l (Wald, KR), as 41 ----
  Sys.setenv(PC_SOURCE_ONLY = "1")
  try(source("calibration/33g_protected_check.R"), silent = TRUE)
  eval(parse(text = grab("calibration/33i_protected_check_all.R", "^## ---- lean statistics", "^## ---- cases")))
  eval(parse(text = grab("calibration/33i_protected_check_all.R", "^setup_case2 <- ", "^run_one <- ")))
  CH <- 10000L
  eval(parse(text = grab("calibration/33l_bigB_all.R", "^COMPV <- ", "^## cases:")))
  stopifnot(exists("exact_any"), exists("setup_case2"), exists("big_one"), exists("setup_step"))
  JOBS <- data.frame(                                       # the replicated cases of 41 (GEMMA's range for delta_hat)
    pop  = c("EUR", "EUR", "YRI", "YRI", "YRI", "YRI", "YRI", "YRI"),
    gene = c("PEX6-region", "ERAP2", "ZNF266", "PEX6-region", "TRA2A-AS", "ZNF266", "PEX6-region", "PEX6-region"),
    pos  = c(42944850L, 96252432L, 9477868L, 42915920L, 23134790L, 9477868L, 43186801L, 43011766L),
    step = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L),
    stat = c("wald", "kr", "wald", "wald", "score", "score", "score", "lrt"), stringsAsFactors = FALSE)
  JOBS$id <- sprintf("%s_%s_%d_%d_%s_crn", JOBS$pop, gsub("[^A-Za-z0-9]", "", JOBS$gene), JOBS$pos, JOBS$step, JOBS$stat)
  GRID2 <- c(0, 10^seq(-5, 5, by = 0.1))
  run_job <- function(k) {
    cs <- JOBS[k, ]; ck <- file.path(OUT_CRN, paste0(cs$id, ".rds"))
    if (file.exists(ck)) return(paste("cached", cs$id))
    t0 <- Sys.time()
    if (cs$stat %in% c("wald", "kr")) {
      S <- setup_step(cs$pop, cs$gene, cs$pos, cs$step)
      f1 <- function(d0) big_one(S$sp, d0, S$obs, cs$stat, 100000L); f2 <- function(d0) big_one(S$sp, d0, S$obs, cs$stat, 200000L)
      level <- S$level
    } else {
      S <- setup_case2(cs$pop, cs$gene, cs$pos, cs$step, cs$stat)
      f1 <- function(d0) exact_any(S, d0, cs$stat, B1); f2 <- function(d0) exact_any(S, d0, cs$stat, B2)
      level <- S$info$level
    }
    prof <- profile2(f1, f2, 45000L + 2L * k); prof$id <- cs$id; i <- which.max(prof$log10p)
    info <- data.frame(id = cs$id, pop = cs$pop, gene = cs$gene, pos = cs$pos, step = cs$step, stat = cs$stat, cap = 1e5,
                       kind = "crn", level = level, new_max = prof$log10p[i], new_argmax = prof$delta0[i],
                       controlled = max(prof$log10up), elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
    saveRDS(list(info = info, profile = prof), ck); sprintf("done %s (%.1f min)", cs$id, info$elapsed_min)
  }
  JOBS <- JOBS[order(JOBS$pop == "EUR", decreasing = TRUE), ]
  print(unlist(parallel::mclapply(seq_len(nrow(JOBS)), function(k) tryCatch(run_job(k), error = function(e) paste("FAIL", JOBS$id[k], conditionMessage(e))),
                                  mc.cores = NC, mc.preschedule = FALSE)))
  quit(save = "no")
}

## ---- summary --------------------------------------------------------------------------------------------
## upper estimate: largest p_hat + z se over the grid; lower value: largest p_hat - z se over the points with
## z * rel_se < 1 (35d uses z = 2 and rel_se < 0.5); z = 2 pointwise, or Bonferroni over the G points visited
bounds <- function(p, simultaneous = FALSE) {
  z <- if (simultaneous) qnorm(pnorm(-2) / nrow(p), lower.tail = FALSE) else 2
  ok <- is.finite(p$log10p) & is.finite(p$rel_se)            # a zero estimate has no upper or lower value
  up <- ifelse(ok, p$log10p + log10(1 + z * p$rel_se), -Inf)
  lo <- ifelse(ok & z * p$rel_se < 1, p$log10p + log10(pmax(1 - z * p$rel_se, 1e-12)), -Inf)
  c(upper = max(up), lower = max(lo), z = z)
}
classify <- function(nom, lv, up, lo) if (nom <= lv) { if (up <= lv) "supported" else if (lo > lv) "not supported" else "undetermined" } else
  { if (up <= lv) "reference" else if (lo > lv) "confirmed" else "undetermined stop" }
key <- function(d) paste(d$pop, gsub("[^A-Za-z0-9]", "", d$gene), d$pos, d$step, d$stat)
rd_dir <- function(dir, pat = "[.]rds$") lapply(list.files(file.path(RES, dir), pat, full.names = TRUE), readRDS)
## refined profiles used in the article: 100,000 draws (33l) for the Wald and KR cases, else 33i
ORIG <- c(rd_dir("protected_check_all"), rd_dir("protected_check_bigB"))
okey <- vapply(ORIG, function(r) key(r$info), "")
DEC <- read.csv(file.path(RES, "gemma_decisions_controlled.csv"), stringsAsFactors = FALSE)   # nominal and level (35d)
DEC$stat <- DEC$path
## (1) the decisive cases: original refined evaluation, replicate (41) and common random numbers (45)
REPL <- rd_dir("replicate_check", "_rep[.]rds$"); CRN <- rd_dir("crn_check")
PATHS_O <- rd_dir("protected_check", "_D1e3[.]rds$")
cases <- unique(do.call(rbind, lapply(CRN, function(r) r$info[, c("pop", "gene", "pos", "step", "stat")])))
rows <- list()
for (q in seq_len(nrow(cases))) {
  cs <- cases[q, ]; isp <- grepl("path", cs$stat); k <- key(cs)
  if (isp) { kk <- paste("YRI", gsub("[^A-Za-z0-9]", "", cs$gene), cs$pos)
    fo <- Filter(function(r) paste(r$info$pop, gsub("[^A-Za-z0-9]", "", r$info$gene), r$info$pos) == kk, PATHS_O)
    fr <- Filter(function(r) grepl("lrpath", r$info$id) && r$info$pos == cs$pos, REPL)
    fc <- Filter(function(r) grepl("lrpath", r$info$id) && r$info$pos == cs$pos, CRN)
    o <- if (length(fo)) fo[[1]] else NULL; rp <- if (length(fr)) fr[[1]] else NULL; cr <- fc[[1]]
    nom <- cr$info$nominal; lv <- cr$info$level
  } else {
    hh <- which(okey == k); o <- if (length(hh)) ORIG[[tail(hh, 1)]] else NULL
    rp <- Filter(function(r) key(r$info) == k, REPL); rp <- if (length(rp)) rp[[1]] else NULL
    cr <- Filter(function(r) key(r$info) == k, CRN)[[1]]
    d <- DEC[DEC$pop == cs$pop & DEC$gene == cs$gene & DEC$pos == cs$pos & DEC$step == cs$step & DEC$stat == cs$stat, ]
    nom <- if (nrow(d)) d$nominal[1] else NA; lv <- if (nrow(d)) d$level[1] else cr$info$level
    if (cs$stat == "kr") { g <- rd_dir("protected_check", "ERAP2_EUR_kr[.]rds$"); if (length(g)) nom <- g[[1]]$info$nominal }
  }
  for (src in c("original", "replicate", "crn")) {
    r <- switch(src, original = o, replicate = rp, crn = cr); if (is.null(r)) next
    for (sim in c(FALSE, TRUE)) { b <- bounds(r$profile, sim)
      rows[[length(rows) + 1]] <- data.frame(cs, level = lv, nominal = nom, source = src, margin = if (sim) "simultaneous" else "pointwise",
        z = b[["z"]], points = nrow(r$profile), max_log10p = max(r$profile$log10p), upper = b[["upper"]], lower = b[["lower"]],
        class = if (is.finite(nom)) classify(nom, lv, b[["upper"]], b[["lower"]]) else NA) } }
}
S <- do.call(rbind, rows)
write.csv(S, file.path(RES, "crn_check_summary.csv"), row.names = FALSE)
options(width = 220)
print(S[, c("pop", "gene", "step", "stat", "level", "nominal", "source", "margin", "points", "max_log10p", "upper", "lower", "class")], digits = 3, row.names = FALSE)
## (2) every entry of GEMMA's paths (35d): classification with pointwise and simultaneous margins
PC <- do.call(rbind, lapply(seq_along(ORIG), function(i) { r <- ORIG[[i]]; b2 <- bounds(r$profile); bs <- bounds(r$profile, TRUE)
  data.frame(key = okey[i], up2 = b2[["upper"]], lo2 = b2[["lower"]], ups = bs[["upper"]], los = bs[["lower"]], nbig = nrow(r$profile)) }))
PC <- PC[!duplicated(PC$key, fromLast = TRUE), ]            # 33l (100,000 draws) listed after 33i replaces it
DEC$key <- key(DEC); m <- match(DEC$key, PC$key)
cat("\nentries of GEMMA's paths without a refined profile:", sum(is.na(m)), "\n")
DEC <- DEC[!is.na(m), ]; m <- m[!is.na(m)]
DEC$class_pointwise <- mapply(classify, DEC$nominal, DEC$level, PC$up2[m], PC$lo2[m])
DEC$class_simultaneous <- mapply(classify, DEC$nominal, DEC$level, PC$ups[m], PC$los[m])
write.csv(DEC, file.path(RES, "simultaneous_decisions.csv"), row.names = FALSE)
cat("\nentries of GEMMA's paths:", nrow(DEC), " matched:", sum(!is.na(m)), " class of 35d reproduced:", sum(DEC$class_pointwise == DEC$class, na.rm = TRUE), "\n")
print(table(pointwise = DEC$class_pointwise, simultaneous = DEC$class_simultaneous))
print(DEC[which(DEC$class_pointwise != DEC$class_simultaneous), c("pop", "gene", "path", "step", "level", "nominal", "class_pointwise", "class_simultaneous")], digits = 3)
