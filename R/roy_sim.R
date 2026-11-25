# ---------------------------------------------------------------------------
# roy_sim: data-generating processes and Monte Carlo experiments
# ---------------------------------------------------------------------------

# invert u -> h(u | v; theta) in its first argument by vectorized bisection
.hinv_first <- function(w, v, fam, theta, iter = 60) {
  lo <- rep(0, length(w)); hi <- rep(1, length(w))
  for (i in seq_len(iter)) {
    mid <- (lo + hi) / 2
    up <- fam$h(pmin(pmax(mid, 1e-12), 1 - 1e-12), v, theta) < w
    lo[up] <- mid[up]; hi[!up] <- mid[!up]
  }
  (lo + hi) / 2
}

#' Simulate Model A (binary AI adoption)
#'
#' Covariates: `x1 ~ N(0,1)`, `x2 ~ Bernoulli(0.5)` and an excluded
#' instrument `z ~ N(0,1)` (e.g. exogenous variation in AI access cost).
#' \eqn{D = 1(\gamma_0 + \gamma_1 x_1 + \gamma_2 x_2 + \gamma_3 z + U_D > 0)},
#' \eqn{Y_s = \beta_{s0} + \beta_{s1} x_1 + \beta_{s2} x_2 + U_s}.
#'
#' With the Gaussian family \eqn{(U_0, U_1, U_D)} is trivariate normal with
#' correlations `rho0`, `rho1`, `rho01`. With another registered family the
#' draws follow a C-vine rooted at \eqn{U_D}: each \eqn{U_s} is drawn from
#' the family's conditional distribution given \eqn{U_D} with parameter
#' `theta0` / `theta1`, and \eqn{U_0 \perp U_1 \mid U_D}. The vine draw
#' assumes an exchangeable copula (true for Gaussian, Clayton, Frank, Gumbel
#' and Student-t); rotated families should supply their own simulator.
#'
#' @param n Sample size.
#' @param beta0,beta1 Outcome coefficients (intercept, x1, x2).
#' @param gamma Selection coefficients (intercept, x1, x2, z).
#' @param sigma0,sigma1 Outcome error standard deviations.
#' @param rho0,rho1,rho01 Gaussian correlations of (U0, UD), (U1, UD), (U0, U1).
#' @param dependence Family name or object.
#' @param theta0,theta1 Dependence parameters for non-Gaussian families.
#' @param seed Optional random seed.
#' @return A data frame with `y`, `d`, `x1`, `x2`, `z`, the latent `y0`,
#'   `y1`, and attribute `"truth"` holding the true parameter vector named as
#'   in `coef(roy(...))`.
#' @export
roy_sim_A <- function(n = 2000, beta0 = c(1, 0.5, 0.3), beta1 = c(1.4, 0.8, 0.3),
                      gamma = c(-0.2, 0.4, 0.2, 0.8), sigma0 = 1, sigma1 = 1.2,
                      rho0 = -0.3, rho1 = 0.4, rho01 = 0.2, dependence = "gaussian",
                      theta0 = NULL, theta1 = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fam <- roy_get_dependence(dependence)
  x1 <- stats::rnorm(n); x2 <- stats::rbinom(n, 1, 0.5); z <- stats::rnorm(n)
  X <- cbind(1, x1, x2); Z <- cbind(1, x1, x2, z)
  if (fam$name == "gaussian") {
    R <- matrix(c(1, rho01, rho0, rho01, 1, rho1, rho0, rho1, 1), 3)
    if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      stop("(rho0, rho1, rho01) do not form a valid correlation matrix")
    E <- matrix(stats::rnorm(3 * n), n) %*% chol(R)
    u0 <- sigma0 * E[, 1]; u1 <- sigma1 * E[, 2]; ud <- E[, 3]
    th0 <- rho0; th1 <- rho1
  } else {
    th0 <- if (is.null(theta0)) fam$start else theta0
    th1 <- if (is.null(theta1)) fam$start else theta1
    v <- stats::runif(n); ud <- stats::qnorm(v)
    u0 <- sigma0 * stats::qnorm(.hinv_first(stats::runif(n), v, fam, th0))
    u1 <- sigma1 * stats::qnorm(.hinv_first(stats::runif(n), v, fam, th1))
  }
  y0 <- drop(X %*% beta0) + u0; y1 <- drop(X %*% beta1) + u1
  d <- as.integer(drop(Z %*% gamma) + ud > 0)
  dat <- data.frame(y = ifelse(d == 1, y1, y0), d = d, x1 = x1, x2 = x2, z = z, y0 = y0, y1 = y1)
  xn <- c("(Intercept)", "x1", "x2")
  pn <- if (fam$name == "gaussian") "rho" else fam$par_names
  truth <- c(stats::setNames(gamma, paste0("sel:", c(xn, "z"))),
             stats::setNames(beta0, paste0("y0:", xn)), stats::setNames(beta1, paste0("y1:", xn)),
             sigma0 = sigma0, sigma1 = sigma1,
             if (fam$npar) stats::setNames(c(th0, th1), paste0(rep(pn, 2), rep(0:1, each = fam$npar))))
  attr(dat, "truth") <- truth
  dat
}

#' Simulate the classical (Roy 1951) and extended Roy models
#'
#' \eqn{Y_s = \beta_{s0} + \beta_{s1} x_1 + \beta_{s2} x_2 + U_s} with
#' \eqn{(U_0, U_1)} bivariate normal, and \eqn{D = 1(Y_1 - Y_0 - \pi w > 0)}.
#' With `pi = 0` this is Roy's original model: workers simply choose the
#' sector with the higher income. `w ~ N(0,1)` is an observed cost shifter
#' (e.g. AI training or subscription cost).
#'
#' @inheritParams roy_sim_A
#' @param rho01 Correlation of (U0, U1).
#' @param pi Cost coefficient on `w`.
#' @return A data frame with attribute `"truth"`.
#' @export
roy_sim_classical <- function(n = 2000, beta0 = c(1, 0.5, 0.3), beta1 = c(1.2, 0.9, 0.1),
                              sigma0 = 1, sigma1 = 1.3, rho01 = 0.3, pi = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rbinom(n, 1, 0.5); w <- stats::rnorm(n)
  X <- cbind(1, x1, x2)
  S <- matrix(c(sigma0^2, rho01 * sigma0 * sigma1, rho01 * sigma0 * sigma1, sigma1^2), 2)
  U <- matrix(stats::rnorm(2 * n), n) %*% chol(S)
  y0 <- drop(X %*% beta0) + U[, 1]; y1 <- drop(X %*% beta1) + U[, 2]
  d <- as.integer(y1 - y0 - pi * w > 0)
  dat <- data.frame(y = ifelse(d == 1, y1, y0), d = d, x1 = x1, x2 = x2, w = w, y0 = y0, y1 = y1)
  xn <- c("(Intercept)", "x1", "x2")
  m <- .classical_map(stats::setNames(beta0, xn), stats::setNames(beta1, xn), sigma0, sigma1, rho01,
                      if (pi != 0) c(w = pi) else numeric(0), xn, c(xn, if (pi != 0) "w"))
  attr(dat, "truth") <- c(stats::setNames(m$gamma, paste0("sel:", names(m$gamma))),
                          stats::setNames(beta0, paste0("y0:", xn)),
                          stats::setNames(beta1, paste0("y1:", xn)),
                          sigma0 = sigma0, sigma1 = sigma1, rho0 = m$rho0, rho1 = m$rho1,
                          rho01 = rho01, sigmaD = m$sigmaD, if (pi != 0) c(`cost:w` = pi))
  dat
}

#' Simulate Model B (human-only / AI-augmented / AI-dominant)
#'
#' Choice: \eqn{S = \arg\max_s \{Z'\gamma_s + \varepsilon_s\}}, i.i.d. Gumbel
#' \eqn{\varepsilon_s}, \eqn{\gamma_0 = 0}, with
#' \eqn{Z = (1, x_1, x_2, z_1, z_2)}; `z1` shifts the cost of AI augmentation
#' and `z2` the cost of AI-dominant production. Outcomes:
#' \eqn{Y_s = X'\beta_s + U_s},
#' \eqn{U_s = \sigma_s(\sum_j r_{sj}\varepsilon^*_j + \sqrt{1 - \sum_j r_{sj}^2}\,\eta_s)},
#' where \eqn{\varepsilon^*_j} are the standardized Gumbel errors. Under this
#' design the Dubin-McFadden correction is exactly correct, and the true
#' correction coefficients are \eqn{\lambda_{sj} = \sigma_s r_{sj}\sqrt{6}/\pi}.
#'
#' @param n Sample size.
#' @param beta List of three outcome coefficient vectors (intercept, x1, x2).
#' @param gamma 5 x 2 matrix of choice coefficients for regimes 1 and 2.
#' @param sigma Outcome error standard deviations (length 3).
#' @param r 3 x 3 loading matrix of outcome errors on choice errors.
#' @param seed Optional random seed.
#' @return A data frame with `y`, `s`, `x1`, `x2`, `z1`, `z2`, latent `y0`,
#'   `y1`, `y2`, and attribute `"truth"`.
#' @export
roy_sim_B <- function(n = 3000,
                      beta = list(c(1, 0.5, 0.3), c(1.3, 0.8, 0.3), c(1.1, 1.0, 0.1)),
                      gamma = cbind(c(0, 0.4, 0.2, 0.8, 0), c(-0.5, 0.6, 0, 0, 0.8)),
                      sigma = c(1, 1.1, 1.3),
                      r = rbind(c(0.5, -0.2, -0.2), c(-0.2, 0.5, -0.1), c(-0.2, -0.1, 0.6)),
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (any(rowSums(r^2) >= 1)) stop("rows of r must have sum of squares < 1")
  x1 <- stats::rnorm(n); x2 <- stats::rbinom(n, 1, 0.5)
  z1 <- stats::rnorm(n); z2 <- stats::rnorm(n)
  X <- cbind(1, x1, x2); Z <- cbind(1, x1, x2, z1, z2)
  eps <- -log(-log(matrix(stats::runif(3 * n), n)))
  gE <- -digamma(1)
  estar <- (eps - gE) / (pi / sqrt(6))
  V <- cbind(0, Z %*% gamma) + eps
  s <- max.col(V, ties.method = "first") - 1L
  Y <- sapply(1:3, function(k) {
    u <- sigma[k] * (drop(estar %*% r[k, ]) + sqrt(1 - sum(r[k, ]^2)) * stats::rnorm(n))
    drop(X %*% beta[[k]]) + u
  })
  dat <- data.frame(y = Y[cbind(seq_len(n), s + 1)], s = s, x1 = x1, x2 = x2, z1 = z1, z2 = z2,
                    y0 = Y[, 1], y1 = Y[, 2], y2 = Y[, 3])
  xn <- c("(Intercept)", "x1", "x2"); zn <- c(xn, "z1", "z2")
  lam <- sigma * r * sqrt(6) / pi
  truth <- c(stats::setNames(as.vector(gamma), paste0("sel", rep(1:2, each = 5), ":", rep(zn, 2))))
  for (k in 1:3) truth <- c(truth, stats::setNames(c(beta[[k]], lam[k, ]),
                                                   paste0("y", k - 1, ":", c(xn, paste0("m", 0:2)))))
  attr(dat, "truth") <- truth
  dat
}

#' Monte Carlo experiments
#'
#' Repeats `sim()` and `fit()` `R` times and summarizes bias, standard
#' deviation, RMSE, mean estimated standard error and coverage of nominal
#' `level` Wald intervals for every parameter present in `truth`.
#'
#' @param sim Function of no arguments returning a data set.
#' @param fit Function taking a data set and returning a `roy_fit`.
#' @param R Number of replications.
#' @param truth Named true parameter vector; defaults to `attr(sim(), "truth")`.
#' @param level Confidence level for coverage.
#' @param seed Optional seed.
#' @param progress Print a progress dot every 10 replications?
#' @return An object of class `roy_mc` with `summary` (data frame), `estimates`
#'   and `se` (R x p matrices) and `failures`.
#' @examples
#' \donttest{
#' mc <- roy_mc(function() roy_sim_A(1000),
#'              function(d) roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = d), R = 20)
#' mc
#' }
#' @export
roy_mc <- function(sim, fit, R = 100, truth = NULL, level = 0.95, seed = NULL, progress = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  est <- se <- NULL; fails <- 0L
  for (r in seq_len(R)) {
    dat <- sim()
    if (is.null(truth)) truth <- attr(dat, "truth")
    f <- tryCatch(suppressWarnings(suppressMessages(fit(dat))), error = function(e) NULL)
    if (is.null(f)) { fails <- fails + 1L; next }
    nm <- intersect(names(truth), names(f$par))
    if (is.null(est)) {
      est <- se <- matrix(NA_real_, R, length(nm), dimnames = list(NULL, nm))
    }
    est[r, ] <- f$par[colnames(est)]
    se[r, ] <- sqrt(pmax(diag(f$vcov)[colnames(est)], 0))
    if (progress && r %% 10 == 0) cat(".")
  }
  if (progress) cat("\n")
  if (is.null(est)) stop("every replication failed")
  ok <- stats::complete.cases(est)
  est <- est[ok, , drop = FALSE]; se <- se[ok, , drop = FALSE]
  tv <- truth[colnames(est)]
  q <- stats::qnorm(1 - (1 - level) / 2)
  cover <- colMeans(abs(sweep(est, 2, tv)) <= q * se, na.rm = TRUE)
  summ <- data.frame(parameter = colnames(est), truth = tv, mean = colMeans(est),
                     bias = colMeans(est) - tv, sd = apply(est, 2, stats::sd),
                     rmse = sqrt(colMeans(sweep(est, 2, tv)^2)),
                     mean_se = colMeans(se, na.rm = TRUE), coverage = cover,
                     row.names = NULL)
  structure(list(summary = summ, estimates = est, se = se, failures = fails + sum(!ok), R = R,
                 level = level), class = "roy_mc")
}

#' @export
print.roy_mc <- function(x, digits = 3, ...) {
  cat("RoyAI Monte Carlo:", x$R, "replications,", x$failures, "failed; coverage at",
      x$level, "\n")
  s <- x$summary; num <- vapply(s, is.numeric, TRUE)
  s[num] <- lapply(s[num], round, digits)
  print(s, row.names = FALSE)
  invisible(x)
}
