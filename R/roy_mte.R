# =============================================================================
# roy_mte.R  --  Treatment parameters for the copula generalized Roy model
# -----------------------------------------------------------------------------
# ATE, ATT, ATU, MTE and LATE for "roy_copula" fits (Model C).
#
# Notation.  W = F_V(V) is the selection uniform (higher W = more inclined
# to AI).  Heckman-Vytlacil "unobserved resistance" is U_S = 1 - W, so that
# D = 1{P(Z) >= U_S} in the binary case.  All treatment parameters are built
# from the copula conditional means
#     m_j(w) = E[U_j | W = w]
#     M_j(p) = int_0^p m_j(w) dw = E[U_j 1{W <= p}]
# which are computed by one-dimensional quadrature of the copula
# h-function, so no copula densities are needed:
#     E[U_j | W = w] = sigma_j [ int_0^inf (1 - G(z|w)) dz - int_-inf^0 G(z|w) dz ]
#     G(z|w)         = dC_j(Phi(z), w)/dw
# For the Gaussian copula m_j(w) = sigma_j rho_j Phi^{-1}(w) (checked in tests).
# =============================================================================

.zgrid <- function() {
  if (is.null(.roy_cache$zgrid)) {
    z <- seq(-9, 9, length.out = 401)       # 400 intervals; z = 0 is node 201
    i0 <- which(abs(z) < 1e-12); z[i0] <- 0
    simp <- function(n) { wt <- rep(c(2, 4), length.out = n); wt[c(1, n)] <- 1; wt }
    .roy_cache$zgrid <- list(z = z, i0 = i0, h = z[2] - z[1],
                             wneg = simp(i0), wpos = simp(length(z) - i0 + 1))
  }
  .roy_cache$zgrid
}

.sgrid <- function() {
  if (is.null(.roy_cache$sgrid)) .roy_cache$sgrid <- seq(-8.5, 8.5, length.out = 241)
  .roy_cache$sgrid
}

#' Conditional mean `E[U_j | W = w]` for a copula with `N(0, sigma^2)` margin
#'
#' Uses `E[Z | w] = int_0^inf (1 - G(z|w)) dz - int_-inf^0 G(z|w) dz` with
#' `G(z|w) = dC(Phi(z), w)/dw`, Simpson's rule on `z` in `[-9, 9]`.
#' @param w selection uniforms
#' @param cp copula parameter list from copula_par()
#' @param sigma outcome error sd
#' @return A numeric vector of conditional means, one per element of `w`.
#' @examples
#' cp <- copula_par("gaussian", tau = 0.4)
#' copula_cond_mean(c(0.1, 0.5, 0.9), cp, sigma = 1.3)
#' @export
copula_cond_mean <- function(w, cp, sigma = 1) {
  g <- .zgrid(); nz <- length(g$z)
  u <- rep(stats::pnorm(g$z), times = length(w))
  G <- matrix(copula_hfunc2(u, rep(w, each = nz), cp), nz, length(w))
  neg <- colSums(G[1:g$i0, , drop = FALSE] * g$wneg)
  pos <- colSums((1 - G[g$i0:nz, , drop = FALSE]) * g$wpos)
  sigma * (pos - neg) * g$h / 3
}

# Tabulate m_j and M_j on the s = Phi^{-1}(w) grid for every regime.
.roy_mgrid <- function(theta, family, J) {
  s <- .sgrid(); w <- stats::pnorm(s); ds <- s[2] - s[1]
  lapply(seq_len(J), function(j) {
    cp <- copula_par(family, theta$tau[j], theta$nu)
    m <- copula_cond_mean(w, cp, theta$sigma[j])
    f <- m * stats::dnorm(s)
    M <- c(0, cumsum((f[-1] + f[-length(f)]) / 2 * ds))
    M <- M - M[length(M)] * (seq_along(M) - 1) / (length(M) - 1) # force M(1) = 0
    list(s = s, m = m, M = M)
  })
}

.interp_M <- function(mg, p) {
  sp <- stats::qnorm(.clamp01(p, 1e-12))
  stats::approx(mg$s, mg$M, xout = sp, rule = 2)$y
}
.interp_m <- function(mg, p) {
  sp <- stats::qnorm(.clamp01(p, 1e-12))
  stats::approx(mg$s, mg$m, xout = sp, rule = 2)$y
}

# Regime intervals on the selection-uniform scale: S = s iff a_s < W <= b_s
.roy_intervals <- function(theta, Z, link) {
  zg <- drop(Z %*% theta$gamma)
  F <- if (link == "probit") stats::pnorm else stats::plogis
  J <- length(theta$kappa) - 1
  a <- sapply(seq_len(J), function(j) if (j == 1) rep(0, length(zg)) else F(theta$kappa[j] - zg))
  b <- sapply(seq_len(J), function(j) if (j == J) rep(1, length(zg)) else F(theta$kappa[j + 1] - zg))
  list(a = matrix(a, ncol = J), b = matrix(b, ncol = J))
}

.fit_link <- function(fit) fit$link

# Core computation of all mean treatment parameters for a parameter vector.
.roy_te_core <- function(par, fit, X = fit$X, Z = fit$Z, to = 1, from = 0) {
  th <- .roy_unpack(par, fit$sp)
  s1 <- to + 1; s0 <- from + 1
  if (identical(fit$selection_type, "mnl") && fit$J > 2) {
    # Lee (1983) single-index model: only (U_s, T_s) pairs are identified, so
    # E[U_r | S = s] for r != s (and hence ATT-type effects) is not.
    return(c(ATE = mean(drop(X %*% (th$beta[[s1]] - th$beta[[s0]])))))
  }
  mg <- .roy_mgrid(th, fit$family, fit$J)
  iv <- .roy_intervals(th, Z, .fit_link(fit))
  xdb <- drop(X %*% (th$beta[[s1]] - th$beta[[s0]]))
  out <- c(ATE = mean(xdb))
  for (k in seq_len(fit$J)) {
    a <- iv$a[, k]; b <- iv$b[, k]; p <- b - a
    dU <- (.interp_M(mg[[s1]], b) - .interp_M(mg[[s1]], a)) -
          (.interp_M(mg[[s0]], b) - .interp_M(mg[[s0]], a))
    out[paste0("TT_s", k - 1)] <- sum(p * xdb + dU) / sum(p)
  }
  if (fit$J == 2) {
    out["ATT"] <- out["TT_s1"]; out["ATU"] <- out["TT_s0"]
    out <- out[c("ATE", "ATT", "ATU")]
  }
  out
}

#' Average treatment parameters: ATE, ATT, ATU (binary) or effects of
#' regime `to` versus `from` conditional on each observed regime (J > 2)
#'
#' @param fit a "roy_copula" fit
#' @param to,from regimes compared (Y_to - Y_from)
#' @param se compute delta-method standard errors
#' @param newdata optional data to average over (default: estimation sample)
#' @return data.frame of estimates and standard errors
#' @examples
#' set.seed(2)
#' d <- roy_copula_sim(1500, list(c(0.5, 0.3), c(1, 0.5)), c(0.2, 0.3, 0.8, 0.6),
#'                     c(0.8, 1), c(-0.4, 0.5), "clayton")
#' fit <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "clayton")
#' roy_te(fit)
#' @export
roy_te <- function(fit, to = 1, from = 0, se = TRUE, newdata = NULL) {
  fit <- .as_roy_copula(fit)
  X <- fit$X; Z <- fit$Z
  if (!is.null(newdata)) {
    X <- stats::model.matrix(stats::delete.response(stats::terms(fit$outcome)), newdata)
    Z <- stats::model.matrix(stats::delete.response(stats::terms(fit$selection)), newdata)
  }
  f <- function(p) .roy_te_core(p, fit, X, Z, to, from)
  est <- f(fit$par)
  sev <- rep(NA_real_, length(est))
  if (se && !is.null(fit$vcov_raw)) {
    G <- numDeriv::jacobian(f, fit$par, method = "simple")
    sev <- sqrt(pmax(diag(G %*% fit$vcov_raw %*% t(G)), 0))
  }
  data.frame(parameter = names(est), estimate = unname(est), se = sev,
             row.names = NULL)
}

#' Marginal treatment effect curve MTE(x, u_S)
#'
#' `MTE(x, u) = x'(beta_to - beta_from) + E[U_to - U_from | U_S = u]`,
#' U_S = 1 - W the unobserved resistance to AI adoption.  For J > 2 (ordered)
#' this is the MTE along the single latent AI-propensity index.
#'
#' @param fit a "roy_copula" fit
#' @param u grid of resistance values in (0, 1)
#' @param x covariate vector (default: sample mean of X)
#' @param to,from regimes compared (MTE of Y_to - Y_from)
#' @param se delta-method pointwise standard errors
#' @return A data frame with columns `u`, `mte` and, if `se = TRUE`, `se`, `lower` and `upper` (pointwise 95% interval).
#' @examples
#' set.seed(2)
#' d <- roy_copula_sim(1500, list(c(0.5, 0.3), c(1, 0.5)), c(0.2, 0.3, 0.8, 0.6),
#'                     c(0.8, 1), c(-0.4, 0.5), "clayton")
#' fit <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "clayton")
#' roy_mte(fit, u = c(0.1, 0.5, 0.9))
#' @export
roy_mte <- function(fit, u = seq(0.01, 0.99, by = 0.01), x = colMeans(.as_roy_copula(fit)$X),
                    to = 1, from = 0, se = TRUE) {
  fit <- .as_roy_copula(fit)
  f <- function(p) {
    th <- .roy_unpack(p, fit$sp)
    mg <- .roy_mgrid(th, fit$family, fit$J)
    sum(x * (th$beta[[to + 1]] - th$beta[[from + 1]])) +
      .interp_m(mg[[to + 1]], 1 - u) - .interp_m(mg[[from + 1]], 1 - u)
  }
  est <- f(fit$par)
  sev <- rep(NA_real_, length(u))
  if (se && !is.null(fit$vcov_raw)) {
    G <- numDeriv::jacobian(f, fit$par, method = "simple")
    sev <- sqrt(pmax(diag(G %*% fit$vcov_raw %*% t(G)), 0))
  }
  data.frame(u = u, mte = est, se = sev, lower = est - 1.96 * sev,
             upper = est + 1.96 * sev)
}

#' Local average treatment effect for an instrument shift (binary model)
#'
#' Compares the propensity score under `data_from` and `data_to` (same
#' people, instrument changed, e.g. AI access z2 from 0 to 1).  The LATE is
#' the average gain among compliers,
#'   `LATE = sum_i [(P1_i - P0_i) x_i' dbeta + G(P1_i) - G(P0_i)] / sum_i (P1_i - P0_i)`
#' with `G(P) = E[(U_1 - U_0) 1{U_S <= P}]`.
#'
#' @param fit binary "roy_copula" fit
#' @param instrument name of the instrument column
#' @param from,to values of the instrument
#' @param data the estimation data (required)
#' @param se compute delta-method standard errors
#' @return A data frame with columns `parameter` (the LATE and the complier share), `estimate` and, if `se = TRUE`, `se`.
#' @examples
#' set.seed(2)
#' d <- roy_copula_sim(1500, list(c(0.5, 0.3), c(1, 0.5)), c(0.2, 0.3, 0.8, 0.6),
#'                     c(0.8, 1), c(-0.4, 0.5), "clayton")
#' fit <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "clayton")
#' roy_late(fit, "z2", data = d)
#' @export
roy_late <- function(fit, instrument, from = 0, to = 1, data = NULL, se = TRUE) {
  fit <- .as_roy_copula(fit)
  if (fit$J != 2) stop("roy_late() is defined for the binary (J = 2) model")
  if (is.null(data)) stop("supply the estimation data via `data`")
  d0 <- data; d1 <- data
  d0[[instrument]] <- from; d1[[instrument]] <- to
  tz <- stats::delete.response(stats::terms(fit$selection))
  Z0 <- stats::model.matrix(tz, d0); Z1 <- stats::model.matrix(tz, d1)
  X  <- stats::model.matrix(stats::delete.response(stats::terms(fit$outcome)), data)
  F <- if (.fit_link(fit) == "probit") stats::pnorm else stats::plogis
  f <- function(p) {
    th <- .roy_unpack(p, fit$sp)
    mg <- .roy_mgrid(th, fit$family, 2)
    P0 <- F(drop(Z0 %*% th$gamma)); P1 <- F(drop(Z1 %*% th$gamma))
    Gf <- function(P) -.interp_M(mg[[2]], 1 - P) + .interp_M(mg[[1]], 1 - P)
    xdb <- drop(X %*% (th$beta[[2]] - th$beta[[1]]))
    c(LATE = sum((P1 - P0) * xdb + Gf(P1) - Gf(P0)) / sum(P1 - P0),
      complier_share = mean(P1 - P0))
  }
  est <- f(fit$par)
  sev <- rep(NA_real_, 2)
  if (se && !is.null(fit$vcov_raw)) {
    G <- numDeriv::jacobian(f, fit$par, method = "simple")
    sev <- sqrt(pmax(diag(G %*% fit$vcov_raw %*% t(G)), 0))
  }
  data.frame(parameter = names(est), estimate = unname(est), se = sev)
}

#' Population ("true") treatment parameters at known parameter values
#'
#' Used by the Monte Carlo study: evaluates ATE/ATT/ATU/LATE at the true
#' structural parameters, averaging over a large covariate sample.
#' @param beta list of J outcome coefficient vectors
#' @param gamma selection coefficients (ordered/probit index)
#' @param sigma outcome error sds
#' @param tau Kendall's tau for each regime
#' @param family copula family
#' @param nu t-copula df
#' @param kappa interior cutpoints beyond 0 (J > 2)
#' @param covariates data.frame of covariates to average over (columns in the
#'   order used by [roy_copula_sim()])
#' @param instrument optional name of a binary instrument for the LATE (0 -> 1)
#' @param to,from regimes compared (Y_to - Y_from)
#' @param link "probit" or "logit"
#' @return named vector of true treatment parameters
#' @examples
#' set.seed(1)
#' cov <- data.frame(x1 = rnorm(5000), z1 = rnorm(5000), z2 = rbinom(5000, 1, 0.5))
#' roy_te_truth(list(c(0.5, 0.3), c(1, 0.5)), c(0.2, 0.3, 0.8, 0.6), c(0.8, 1),
#'              c(-0.3, 0.4), "clayton", covariates = cov, instrument = "z2")
#' @export
roy_te_truth <- function(beta, gamma, sigma, tau, family, nu = NULL,
                         kappa = NULL, covariates, instrument = NULL,
                         to = 1, from = 0, link = "probit") {
  J <- length(beta)
  th <- list(beta = beta, gamma = gamma, kappa = c(-Inf, 0, kappa, Inf),
             sigma = sigma, tau = tau, nu = nu)
  X <- cbind(1, as.matrix(covariates[, seq_len(length(beta[[1]]) - 1), drop = FALSE]))
  Z <- cbind(1, as.matrix(covariates[, seq_len(length(gamma) - 1), drop = FALSE]))
  mg <- .roy_mgrid(th, family, J)
  iv <- .roy_intervals(th, Z, link)
  s1 <- to + 1; s0 <- from + 1
  xdb <- drop(X %*% (beta[[s1]] - beta[[s0]]))
  out <- c(ATE = mean(xdb))
  for (k in seq_len(J)) {
    a <- iv$a[, k]; b <- iv$b[, k]; p <- b - a
    dU <- (.interp_M(mg[[s1]], b) - .interp_M(mg[[s1]], a)) -
          (.interp_M(mg[[s0]], b) - .interp_M(mg[[s0]], a))
    out[paste0("TT_s", k - 1)] <- sum(p * xdb + dU) / sum(p)
  }
  if (J == 2) { out["ATT"] <- out["TT_s1"]; out["ATU"] <- out["TT_s0"]; out <- out[c("ATE", "ATT", "ATU")] }
  if (!is.null(instrument) && J == 2) {
    F <- if (link == "probit") stats::pnorm else stats::plogis
    c0 <- covariates; c1 <- covariates; c0[[instrument]] <- 0; c1[[instrument]] <- 1
    Z0 <- cbind(1, as.matrix(c0[, seq_len(length(gamma) - 1), drop = FALSE]))
    Z1 <- cbind(1, as.matrix(c1[, seq_len(length(gamma) - 1), drop = FALSE]))
    P0 <- F(drop(Z0 %*% gamma)); P1 <- F(drop(Z1 %*% gamma))
    Gf <- function(P) -.interp_M(mg[[2]], 1 - P) + .interp_M(mg[[1]], 1 - P)
    out["LATE"] <- sum((P1 - P0) * xdb + Gf(P1) - Gf(P0)) / sum(P1 - P0)
  }
  out
}
