# ==============================================================================
# 25_verify_proofs.R
# Numerical verification of the exact identities used in the derivations
# (proofs), on random instances of the reduced model
# s ~ N(xi beta, sigma^2 diag(1 + delta l_i)):
#   A.2   Cov(b, r) = 0 and Var r = sigma^2 (D - xi xi'/S0)
#   A.3   Q(delta; r + xi kappa u) = sqrt(nu) z / sqrt(RSS/S + z^2)
#   A.4   l0 = ltilde + nu/2 log(1 - Q^2/nu);  l1 = l0 + nu/2 log(1 + t^2/(nu-1)) + 1/2 log(RSS/S)
#   4(a)  REML efficient information on delta: null model (nu/2) V_U, full model
#         (nu/2) V_U - V_p - nu Delta^2 / (2 (nu-1)), against direct matrix traces
#   4(b)  Giesbrecht-Burns df with the null-model information = (1/nu + 1/nu_j)^(-1),
#         and with the full-model information = closed form f1
#   4(c)  Var G = V_p exactly (G = -sum a w xi r / (sqrt(S0) sigma))
# Usage (from the repository root): Rscript calibration/25_verify_proofs.R
# ==============================================================================
set.seed(20260930)
tol <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))
ok <- function(name, err, eps = 1e-8) cat(sprintf("%-58s max rel. error %.2e  %s\n", name, err,
                                                   if (err < eps) "OK" else "FAIL"))
for (rep in 1:3) {
  nu <- sample(15:60, 1); ell <- c(rexp(nu - 5, 1 / 3), rep(0, 5)); xi <- rnorm(nu) * sqrt(ell + 0.2)
  delta0 <- runif(1, 0.05, 2); sigma2 <- runif(1, 0.5, 2); beta <- rnorm(1)
  cat(sprintf("\n-- instance %d: nu = %d, delta0 = %.3f, sigma2 = %.3f\n", rep, nu, delta0, sigma2))
  w <- 1 / (1 + delta0 * ell); D <- diag(1 / w); W <- diag(w); S0 <- sum(w * xi^2)
  a <- ell * w; p <- w * xi^2 / S0
  aU <- mean(a); VU <- mean(a^2) - aU^2; ap <- sum(p * a); Vp <- sum(p * a^2) - ap^2; Dl <- ap - aU
  nuj <- nu * VU / Dl^2; kap <- 2 * Vp / (nu * VU)
  ## A.2: b = xi' W s / S0, r = s - xi b = (I - xi xi' W / S0) s
  M <- diag(nu) - xi %*% t(xi) %*% W / S0; Vs <- sigma2 * D
  ok("A.2  Cov(b, r) = 0", max(abs(M %*% Vs %*% W %*% xi / S0)))
  ok("A.2  Var r = sigma2 (D - xi xi'/S0)", tol(M %*% Vs %*% t(M), sigma2 * (D - xi %*% t(xi) / S0)))
  ok("A.2  xi' W r = 0 for any s", max(abs(t(xi) %*% W %*% M)))
  ## draw data, build r
  s <- xi * beta + sqrt(sigma2 / w) * rnorm(nu); b <- sum(w * xi * s) / S0; r <- s - xi * b
  kap_u <- sqrt(sigma2 / S0)
  Afun <- function(dl, y) sum(y^2 / (1 + dl * ell)); Sfun <- function(dl) sum(xi^2 / (1 + dl * ell))
  Nfun <- function(dl, y) sum(xi * y / (1 + dl * ell))
  RSSf <- function(dl, y) Afun(dl, y) - Nfun(dl, y)^2 / Sfun(dl)
  Qf <- function(dl, y) Nfun(dl, y) / sqrt(Sfun(dl) * Afun(dl, y) / nu)
  tf <- function(dl, y) (Nfun(dl, y) / Sfun(dl)) * sqrt((nu - 1) * Sfun(dl) / RSSf(dl, y))
  L <- function(dl) sum(log1p(dl * ell))
  l0 <- function(dl, y) -0.5 * (L(dl) + nu * log(Afun(dl, y)))
  l1 <- function(dl, y) -0.5 * (L(dl) + log(Sfun(dl)) + (nu - 1) * log(RSSf(dl, y)))
  lt <- function(dl, y) -0.5 * (L(dl) + nu * log(RSSf(dl, y)))
  ## A.3 and A.4 on a grid of (delta, u)
  e3 <- e4a <- e4b <- 0
  for (dl in c(0, 0.01, 0.3, 1, 5)) for (u in c(-7, -1.5, 0, 0.8, 4)) {
    y <- r + xi * kap_u * u; z <- Nfun(dl, r) / Sfun(dl) + kap_u * u
    e3 <- max(e3, abs(Qf(dl, y) - sqrt(nu) * z / sqrt(RSSf(dl, r) / Sfun(dl) + z^2)))
    e4a <- max(e4a, abs(l0(dl, y) - (lt(dl, r) + nu / 2 * log(1 - Qf(dl, y)^2 / nu))))
    e4b <- max(e4b, abs(l1(dl, y) - (l0(dl, y) + nu / 2 * log1p(tf(dl, y)^2 / (nu - 1)) +
                                      0.5 * log(RSSf(dl, r) / Sfun(dl)))))
  }
  ok("A.3  Q(delta; r + xi kappa u) closed form", e3)
  ok("A.4  l0 = ltilde + nu/2 log(1 - Q^2/nu)", e4a)
  ok("A.4  l1 = l0 + nu/2 log(1 + t^2/(nu-1)) + 1/2 log(RSS/S)", e4b)
  ## 4(a): information by direct traces, parameters (delta, sigma2)
  V <- sigma2 * D; Vi <- solve(V); dVd <- sigma2 * diag(ell); dVs <- D
  info <- function(P) { m <- list(P %*% dVd, P %*% dVs)
    outer(1:2, 1:2, Vectorize(function(k, l) 0.5 * sum(diag(m[[k]] %*% m[[l]])))) }
  I0 <- info(Vi)
  Pi <- Vi - Vi %*% xi %*% t(xi) %*% Vi / c(t(xi) %*% Vi %*% xi); I1 <- info(Pi)
  eff <- function(I) I[1, 1] - I[1, 2]^2 / I[2, 2]
  ok("4(a) null-model efficient information = (nu/2) V_U", abs(eff(I0) - nu / 2 * VU) / (nu / 2 * VU))
  i1 <- nu / 2 * VU - Vp - nu * Dl^2 / (2 * (nu - 1))
  ok("4(a) full-model efficient information, closed form", abs(eff(I1) - i1) / i1)
  ok("4(a) full-model = (nu/2)V_U (1 - kappa - nu/((nu-1) nu_j))",
     abs(i1 - nu / 2 * VU * (1 - kap - nu / ((nu - 1) * nuj))) / i1)
  ## 4(b): Giesbrecht-Burns df, gradient of log V_b = log sigma2 - log S(delta)
  g <- c(ap, 1 / sigma2)
  fGB0 <- 2 / c(t(g) %*% solve(I0) %*% g); fGB1 <- 2 / c(t(g) %*% solve(I1) %*% g)
  ok("4(b) GB df, null-model information = (1/nu + 1/nu_j)^-1", abs(fGB0 - 1 / (1 / nu + 1 / nuj)) / fGB0)
  f1 <- ((nu - 1) * (nu * VU - 2 * Vp) - nu * Dl^2) / (nu * (VU + Dl^2) - 2 * Vp)
  ok("4(b) GB df, full-model information, closed form", abs(fGB1 - f1) / fGB1)
  ok("4(b) gradient: -S'(delta0)/S0 = abar_p",
     abs(-(Sfun(delta0 + 1e-7) - Sfun(delta0 - 1e-7)) / (2e-7) / S0 - ap) / abs(ap), 1e-6)
  ## 4(c): Var G = V_p exactly, G = -sum(a w xi r) / (sqrt(S0) sigma)
  cvec <- -a * w * xi / (sqrt(S0) * sqrt(sigma2))
  VarG <- c(t(cvec) %*% (sigma2 * (D - xi %*% t(xi) / S0)) %*% cvec)
  ok("4(c) Var G = V_p", abs(VarG - Vp) / Vp)
  dbeta <- (sum(xi * r / (1 + (delta0 + 1e-7) * ell)) / Sfun(delta0 + 1e-7) -
            sum(xi * r / (1 + (delta0 - 1e-7) * ell)) / Sfun(delta0 - 1e-7)) / 2e-7
  ok("4(c) d/d delta beta_hat(delta0; r) = -sum(a w xi r)/S0",
     abs(dbeta - (-sum(a * w * xi * r) / S0)) / max(1e-12, abs(dbeta)), 1e-5)
}
