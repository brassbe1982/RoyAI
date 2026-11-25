# ---------------------------------------------------------------------------
# roy_select: binary and multinomial selection stages and correction terms
# ---------------------------------------------------------------------------

#' Estimate the selection equation
#'
#' Model A uses a probit, \eqn{D = 1(Z'\gamma + U_D > 0)} with
#' \eqn{U_D \sim N(0,1)}. Model B uses a multinomial logit,
#' \eqn{S = \arg\max_s \{Z'\gamma_s + \varepsilon_s\}} with i.i.d. Gumbel
#' \eqn{\varepsilon_s} and \eqn{\gamma_0 = 0} (human-only is the base).
#'
#' @param spec A [roy_spec()] object.
#' @param control List passed to [stats::optim()] for the multinomial logit.
#' @return An object of class `roy_selection` with `coef` (vector for Model A,
#'   `ncol(Z)` x `(J-1)` matrix for Model B), `vcov`, `prob` (n x J matrix of
#'   fitted choice probabilities), `loglik` and `type`.
#' @export
roy_select <- function(spec, control = list()) {
  stopifnot(inherits(spec, "roy_spec"))
  if (spec$J == 2) .fit_probit(spec$d, spec$Z) else .fit_mnl(spec$d, spec$Z, spec$J, control)
}

.fit_probit <- function(d, Z) {
  g <- suppressWarnings(stats::glm.fit(Z, d, family = stats::binomial(link = "probit")))
  coef <- g$coefficients; names(coef) <- colnames(Z)
  if (anyNA(coef)) stop("probit selection equation is rank deficient")
  eta <- drop(Z %*% coef)
  p1 <- stats::pnorm(eta)
  # expected information for the probit
  w <- stats::dnorm(eta)^2 / pmax(p1 * (1 - p1), 1e-300)
  V <- solve(crossprod(Z * sqrt(w)))
  dimnames(V) <- list(colnames(Z), colnames(Z))
  ll <- sum(stats::pnorm(ifelse(d == 1, eta, -eta), log.p = TRUE))
  structure(list(type = "probit", coef = coef, vcov = V, prob = cbind(1 - p1, p1),
                 index = eta, loglik = ll, converged = g$converged, n = length(d), J = 2L,
                 npar = length(coef)), class = "roy_selection")
}

# probabilities of an MNL with base category 0; G is kz x (J-1)
.mnl_prob <- function(G, Z) {
  V <- cbind(0, Z %*% G)
  V <- V - apply(V, 1, max)
  E <- exp(V)
  E / rowSums(E)
}

.fit_mnl <- function(d, Z, J, control = list()) {
  kz <- ncol(Z); Y <- outer(d, 0:(J - 1), "==") * 1
  nll <- function(th) {
    P <- .mnl_prob(matrix(th, kz, J - 1), Z)
    -sum(log(pmax(P[cbind(seq_along(d), d + 1L)], 1e-300)))
  }
  gr <- function(th) {
    P <- .mnl_prob(matrix(th, kz, J - 1), Z)
    -as.vector(crossprod(Z, Y[, -1, drop = FALSE] - P[, -1, drop = FALSE]))
  }
  ctrl <- utils::modifyList(list(maxit = 1000, reltol = 1e-12), control)
  o <- stats::optim(rep(0, kz * (J - 1)), nll, gr, method = "BFGS", control = ctrl)
  G <- matrix(o$par, kz, J - 1, dimnames = list(colnames(Z), paste0("alt", 1:(J - 1))))
  P <- .mnl_prob(G, Z)
  # analytic Hessian of the MNL log-likelihood
  H <- matrix(0, kz * (J - 1), kz * (J - 1))
  for (a in 1:(J - 1)) for (b in 1:(J - 1)) {
    w <- P[, a + 1] * ((a == b) - P[, b + 1])
    H[(a - 1) * kz + 1:kz, (b - 1) * kz + 1:kz] <- crossprod(Z * w, Z)
  }
  V <- solve(H)
  nm <- paste(rep(colnames(G), each = kz), rownames(G), sep = ":")
  dimnames(V) <- list(nm, nm)
  structure(list(type = "mnl", coef = G, vcov = V, prob = P, loglik = -o$value,
                 converged = o$convergence == 0, n = length(d), J = J,
                 npar = length(o$par)), class = "roy_selection")
}

#' @export
print.roy_selection <- function(x, ...) {
  cat("<roy_selection>", x$type, "| logLik", format(x$loglik, digits = 6),
      "| converged:", x$converged, "\n")
  print(x$coef, digits = 4)
  invisible(x)
}

#' Selection-correction terms
#'
#' `roy_imr()` returns the inverse Mills ratios for Model A:
#' \eqn{\lambda_1 = \phi(c)/\Phi(c)} for \eqn{D=1} and
#' \eqn{\lambda_0 = -\phi(c)/\Phi(-c)} for \eqn{D=0}, with \eqn{c = Z'\gamma}.
#'
#' `roy_dmf_terms()` returns, for regime `s`, the Bourguignon-Fournier-Gurgand
#' (2007) version of the Dubin-McFadden terms
#' \eqn{m_j = E[\varepsilon_j - \gamma_E \mid S = s]}: \eqn{-\log P_s} for
#' \eqn{j = s} and \eqn{P_j \log P_j / (1 - P_j)} for \eqn{j \neq s}.
#'
#' `roy_lee_terms()` returns Lee's (1983) term
#' \eqn{-\phi(\Phi^{-1}(P_s))/P_s}.
#'
#' @param index Selection index \eqn{Z'\gamma}.
#' @param d Selection indicator (0/1).
#' @param P n x J matrix of choice probabilities.
#' @param s Regime (0-based).
#' @return A numeric vector (`roy_imr`, `roy_lee_terms`) or an n x J matrix
#'   (`roy_dmf_terms`).
#' @export
roy_imr <- function(index, d) {
  ifelse(d == 1,
         exp(stats::dnorm(index, log = TRUE) - stats::pnorm(index, log.p = TRUE)),
         -exp(stats::dnorm(index, log = TRUE) - stats::pnorm(-index, log.p = TRUE)))
}

#' @rdname roy_imr
#' @export
roy_dmf_terms <- function(P, s) {
  P <- pmin(pmax(P, 1e-12), 1 - 1e-12)
  M <- P * log(P) / (1 - P)
  M[, s + 1] <- -log(P[, s + 1])
  colnames(M) <- paste0("m", seq_len(ncol(P)) - 1L)
  M
}

#' @rdname roy_imr
#' @export
roy_lee_terms <- function(P, s) {
  p <- pmin(pmax(P[, s + 1], 1e-12), 1 - 1e-12)
  -stats::dnorm(stats::qnorm(p)) / p
}

#' Fitted selection probabilities
#'
#' @param fit A `roy_fit` object.
#' @param newdata Optional data frame.
#' @return n x J matrix of probabilities of each regime.
#' @export
roy_propensity <- function(fit, newdata = NULL) {
  Z <- .design_new(fit$spec, newdata)$Z
  g <- fit$coef$selection
  if (fit$spec$J == 2) {
    p1 <- stats::pnorm(drop(Z %*% g)); P <- cbind(1 - p1, p1)
  } else P <- .mnl_prob(g, Z)
  colnames(P) <- fit$spec$regimes
  P
}
