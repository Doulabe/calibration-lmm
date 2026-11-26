# 53_fig_nu_tails.R
#
# Purpose. Illustration of the spectral index nu_j and of the tail distortion of the Wald
# statistic on a stylized spectrum (Figure 2 of the short version of the article, and
# slide 8 of the presentation). No simulation output is read.
#
# Left panel. Eigenvalues l_i of a stylized regional kernel (nu = 40 directions, mean one),
# the sensitivities a_i = l_i / (1 + delta l_i) at delta = 0.25, the loadings p_i of two
# candidates (A concentrated on the directions of large background variance, B spread over
# the spectrum) and the resulting nu_j = nu V_U / Delta^2, Delta = abar_p - abar_U.
#
# Right panel. Two-sided null tail probability of the Wald statistic against the threshold c
# under the local approximation t(delta_hat) ~ N(0, 1 + rho(c)), rho(c) = kappa_j + c^2 / (2 nu_j)
# with kappa_j = 0.05, for nu_j in {10, 30, 100, 1000}, against the Student reference (nu large),
# with the variance factor and the tail ratio at alpha = 5e-8.
#
# Usage (from the repository root):
#   Rscript calibration/53_fig_nu_tails.R en calibration/output/figures/fig_nu_tails.pdf
#   Rscript calibration/53_fig_nu_tails.R fr presentation/fig_nu_tails.pdf
# Fonts are embedded with Ghostscript when it is available.

args <- commandArgs(trailingOnly = TRUE)
lang <- if (length(args) >= 1) args[1] else "en"
out  <- if (length(args) >= 2) args[2] else "fig_nu_tails.pdf"
stopifnot(lang %in% c("en", "fr"))
L <- if (lang == "fr") list(
  xlab1 = "directions propres i, variance de fond decroissante",
  leg_a = "(echelle de gauche)", leg_p = "charges", of = "du candidat",
  xlab2 = expression("seuil" ~ c ~ "sur" ~ "|" * t * "|"),
  ylab2 = expression(P * "(|" * t(hat(delta)) * "| >" ~ c * ")" ~ "sous" ~ H[0]),
  ref = "reference", variance = "variance x", tail = "queue x", big = " "
) else list(
  xlab1 = "eigen-directions i, decreasing background variance",
  leg_a = "(left scale)", leg_p = "loadings", of = "of candidate",
  xlab2 = expression("threshold" ~ c ~ "on" ~ "|" * t * "|"),
  ylab2 = expression(P * "(|" * t(hat(delta)) * "| >" ~ c * ")" ~ "under" ~ H[0]),
  ref = "reference", variance = "scale x", tail = "tail x", big = ","
)

nu    <- 40
delta <- 0.25
kappa <- 0.05
l <- 12 * exp(-0.45 * (0:(nu - 1))) + 0.05
l <- l / mean(l)
a <- l / (1 + delta * l)
aU <- mean(a); VU <- mean((a - aU)^2)
pA <- l^1.5; pA <- pA / sum(pA)
pB <- ifelse(seq_len(nu) <= 20, 1, 1.5); pB <- pB / sum(pB)
nuj <- function(p) { ap <- sum(p * a); c(ap = ap, nu = nu * VU / (ap - aU)^2) }
A <- nuj(pA); B <- nuj(pB)

col_gray <- "#B4B2A9"; col_A <- "#534AB7"; col_B <- "#0F6E56"
cols <- c(`10` = "#633806", `30` = "#BA7517", `100` = "#EF9F27", `1000` = "#F5C98A")

pdf(out, width = 8.5, height = 3.7, useDingbats = FALSE)   # scaled to the text width, labels stay above 7 pt
par(mfrow = c(1, 2), mar = c(4.0, 4.0, 1.4, 0.8), mgp = c(2.3, 0.6, 0), cex = 0.9)

## Left panel
bp <- barplot(a, col = col_gray, border = NA, ylim = c(0, 3.6), yaxs = "i",
              xlab = L$xlab1, ylab = expression(a[i] == l[i] / (1 + delta * l[i])),
              names.arg = rep("", nu))
abline(h = aU, lty = 2, col = "#5F5E5A")
abline(h = A["ap"], lty = 2, col = col_A)
abline(h = B["ap"], lty = 2, col = col_B)
sc <- 3.6 / 0.55
points(bp, pA * sc, pch = 19, col = col_A, cex = 0.9)
points(bp, pB * sc, pch = 15, col = col_B, cex = 0.9)
text(bp[nu], aU + 0.14, expression(bar(a)[U]), adj = 1, col = "#5F5E5A")
text(bp[nu], A["ap"] + 0.14, bquote(bar(a)[p] ~ .(L$of) ~ "A"), adj = 1, col = col_A)
text(bp[16], 0.62, bquote(bar(a)[p] ~ .(L$of) ~ "B"), adj = 0.5, col = col_B)
legend(bp[9], 2.55, bty = "n", cex = 0.9, xjust = 0, yjust = 1,
       pch = c(15, 19, 15), col = c(col_gray, col_A, col_B), pt.cex = c(1.4, 0.9, 0.9),
       legend = c(bquote(a[i] ~ .(L$leg_a)),
                  bquote(.(L$leg_p) ~ p[i] ~ .(L$of) ~ "A :" ~ nu[j] %~~% .(round(A["nu"], 1))),
                  bquote(.(L$leg_p) ~ p[i] ~ .(L$of) ~ "B :" ~ nu[j] %~~% .(round(B["nu"])))))

## Right panel
cs <- seq(2, 6.5, by = 0.02)
tail_ref <- function(c) 2 * pnorm(-c)
tail_loc <- function(c, n) { rho <- kappa + c^2 / (2 * n); 2 * pnorm(-c / sqrt(1 + rho)) }
plot(cs, log10(tail_ref(cs)), type = "n", ylim = c(-11, 0), xlim = c(2, 7.5), yaxt = "n",
     xlab = L$xlab2, ylab = L$ylab2)
axis(2, at = seq(-10, 0, 2), labels = c(expression(10^-10), expression(10^-8), expression(10^-6),
                                         expression(10^-4), expression(10^-2), "1"), las = 1)
abline(h = -10:0, col = "#E1E0D9")
levels <- c(1e-3, 1e-5, 5e-8)
lab <- c(expression(10^-3), expression(10^-5), expression(5 %.% 10^-8))
for (k in seq_along(levels)) {
  ck <- qnorm(1 - levels[k] / 2)
  abline(v = ck, lty = 3, col = "#5F5E5A")
  text(ck, 0, lab[[k]], pos = 3, cex = 0.85, xpd = TRUE, offset = 0.2)
}
lines(cs, log10(tail_ref(cs)), lty = 2, lwd = 2, col = "#5F5E5A")
for (n in c(1000, 100, 30, 10)) lines(cs, log10(tail_loc(cs, n)), lwd = 2.2, col = cols[as.character(n)])
for (n in c(10, 30, 100, 1000)) text(6.55, log10(tail_loc(6.5, n)), bquote(nu[j] == .(n)), adj = 0, cex = 0.9)
text(6.55, log10(tail_ref(6.5)), L$ref, adj = 0, cex = 0.9, col = "#5F5E5A")
c58 <- qnorm(1 - 5e-8 / 2)
txt <- lapply(c(10, 30, 100, 1000), function(n) {
  rho <- kappa + c58^2 / (2 * n); r <- tail_loc(c58, n) / 5e-8
  rtxt <- if (r >= 100) format(round(r), big.mark = L$big) else sprintf("%.1f", r)
  bquote(nu[j] == .(n) ~ ":" ~ .(L$variance) ~ .(sprintf("%.2f", 1 + rho)) * "," ~ .(L$tail) ~ .(rtxt))
})
legend("bottomleft", inset = c(0.02, 0.03), bty = "n", cex = 0.85,
       title = expression(alpha == 5 %.% 10^-8), legend = as.expression(txt))
invisible(dev.off())
## Embed Helvetica and Symbol (plotmath Greek letters): Ghostscript leaves the base-14 fonts
## unembedded unless its NeverEmbed list is cleared.
if (nzchar(Sys.which("gs"))) {
  tmp <- paste0(out, ".tmp")
  st <- system2("gs", c("-q", "-o", shQuote(tmp), "-sDEVICE=pdfwrite", "-dEmbedAllFonts=true",
                        "-dSubsetFonts=true", "-dCompatibilityLevel=1.5",
                        "-c", shQuote("<</NeverEmbed []>> setdistillerparams"), "-f", shQuote(out)))
  if (st == 0 && file.exists(tmp)) file.rename(tmp, out) else unlink(tmp)
}
cat(sprintf("%s written, nu_j(A) = %.2f, nu_j(B) = %.0f\n", out, A["nu"], B["nu"]))
