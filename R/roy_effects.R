# ---------------------------------------------------------------------------
# Potential outcomes and average treatment parameters
# (MTE / LATE / PRTE live in roy_mte, which builds on roy_potential())
# ---------------------------------------------------------------------------

# Gauss-Hermite nodes and weights for the standard normal (Golub-Welsch)
.gh_normal <- function(n = 60) {
  Jm <- matrix(0, n, n)
  off <- sqrt(seq_len(n - 1))
  Jm[cbind(1:(n - 1), 2:n)] <- off; Jm[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(Jm, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

# E[e_s | D = k] for standardized outcome error e_s and selection index c,
# integrating the dependence family's h-function numerically.
.cond_mean_quad <- function(cidx, fam, theta, k, nodes = 60) {
  gh <- .gh_normal(nodes)
  eps <- 1e-12
  u <- pmin(pmax(stats::pnorm(gh$x), eps), 1 - eps)
  v <- pmin(pmax(stats::pnorm(-cidx), eps), 1 - eps)
  # P(D = 0 | e = x) = h(v | u)
  H <- vapply(u, function(ui) fam$h(v, rep(ui, length(v)), theta), numeric(length(v)))
  H <- matrix(H, nrow = length(v))
  Pk <- if (k == 1) 1 - H else H
  num <- drop(Pk %*% (gh$w * gh$x)); den <- drop(Pk %*% gh$w)
  num / den
}

.potential_A <- function(co, fam, X, Z, regime, given) {
  s <- as.character(regime)
  mu <- drop(X %*% co$outcome[[s]])
  if (is.null(given)) return(mu)
  cidx <- drop(Z %*% co$selection)
  if (!is.null(co$cov_sd)) {
    return(mu + co$cov_sd[[s]] * roy_imr(cidx, rep(given, length(cidx))))
  }
  if (!is.null(fam$cond_mean))
    return(mu + co$sigma[[s]] * fam$cond_mean(cidx, co$dependence[[s]], select = given == 1))
  mu + co$sigma[[s]] * .cond_mean_quad(cidx, fam, co$dependence[[s]], given)
}

.potential_B <- function(co, spec, correction, X, Z, regime, given) {
  s <- as.character(regime)
  mu <- drop(X %*% co$outcome[[s]])
  if (is.null(given)) return(mu)
  P <- .mnl_prob(co$selection, Z)
  if (correction == "dmf") {
    M <- roy_dmf_terms(P, given)
    return(mu + drop(M %*% co$correction[[s]]))
  }
  if (regime != given) return(rep(NA_real_, length(mu)))
  mu + co$correction[[s]] * roy_lee_terms(P, given)
}

#' Potential outcomes
#'
#' Returns \eqn{E[Y_s \mid X]} (`given = NULL`) or
#' \eqn{E[Y_s \mid X, Z, S = k]} (`given = k`) for every row. With
#' `regime != given` this is the counterfactual mean that the Roy model is
#' built to recover: what an adopter would have earned without AI, and so on.
#'
#' For Model A with a non-Gaussian dependence family the conditional mean is
#' integrated numerically from the family's h-function. For Model B the
#' counterfactual means require `correction = "dmf"`; with Lee's correction
#' only own-regime conditional means are identified and `NA` is returned.
#'
#' @param fit A `roy_fit` object.
#' @param regime Potential outcome to predict (0-based regime code).
#' @param given Conditioning regime, or `NULL` for the unconditional mean.
#' @param newdata Optional data frame (defaults to the estimation sample).
#' @param coef_list Internal: coefficient list to use instead of `fit$coef`.
#' @return Numeric vector.
#' @export
roy_potential <- function(fit, regime, given = NULL, newdata = NULL, coef_list = NULL) {
  dm <- .design_new(fit$spec, newdata)
  co <- if (is.null(coef_list)) fit$coef else coef_list
  if (fit$model == "A") .potential_A(co, fit$family, dm$X, dm$Z, regime, given)
  else .potential_B(co, fit$spec, fit$type, dm$X, dm$Z, regime, given)
}

#' @rdname roy_fit_methods
#' @param newdata Optional data frame.
#' @param type For `predict()`: `"observed"` (\eqn{E[Y \mid X, Z, S]} for
#'   the observed regime; needs the selection variable), `"potential"`
#'   (matrix of unconditional potential outcomes) or `"propensity"`.
#' @export
predict.roy_fit <- function(object, newdata = NULL, type = c("potential", "observed", "propensity"), ...) {
  type <- match.arg(type)
  J <- object$spec$J
  if (type == "propensity") return(roy_propensity(object, newdata))
  if (type == "potential") {
    M <- vapply(0:(J - 1), function(s) roy_potential(object, s, NULL, newdata), numeric(
      if (is.null(newdata)) object$spec$n else nrow(newdata)))
    M <- matrix(M, ncol = J); colnames(M) <- object$spec$regimes
    return(M)
  }
  d <- if (is.null(newdata)) object$spec$d else {
    dv <- newdata[[all.vars(object$spec$selection)[1]]]
    match(as.character(dv), object$spec$levels) - 1L
  }
  out <- numeric(length(d))
  for (s in 0:(J - 1)) {
    i <- d == s
    if (any(i)) out[i] <- roy_potential(object, s, s, newdata)[i]
  }
  out
}

#' Average treatment parameters
#'
#' For Model A: ATE, ATT (adopters) and ATU (non-adopters) of AI adoption,
#' \eqn{E[Y_1 - Y_0 \mid \cdot]}. For Model B: for every pair of regimes
#' \eqn{(k, j)} with \eqn{k > j}, the average effect of moving from regime
#' \eqn{j} to \eqn{k} in the population and within each observed regime.
#' Standard errors use the delta method on the fitted parameter covariance
#' (or the bootstrap draws when the fit used `se = "bootstrap"`).
#'
#' @param fit A `roy_fit` object.
#' @param se Compute standard errors?
#' @return A data frame with columns `contrast`, `population`, `estimate`,
#'   `se`.
#' @export
roy_effects <- function(fit, se = TRUE) {
  spec <- fit$spec; J <- spec$J; d <- spec$d
  pairs <- if (J == 2) list(c(1, 0)) else {
    cmb <- utils::combn(0:(J - 1), 2); lapply(seq_len(ncol(cmb)), function(i) rev(cmb[, i]))
  }
  rows <- list()
  for (pr in pairs) {
    k <- pr[1]; j <- pr[2]
    lab <- sprintf("%s vs %s", spec$regimes[k + 1], spec$regimes[j + 1])
    rows[[length(rows) + 1]] <- list(lab, "all", k, j, NULL)
    for (g in 0:(J - 1)) {
      pop <- if (J == 2) c("non-adopters (ATU)", "adopters (ATT)")[g + 1] else paste0("S = ", spec$regimes[g + 1])
      rows[[length(rows) + 1]] <- list(lab, pop, k, j, g)
    }
  }
  if (J == 2) rows[[1]][[2]] <- "all (ATE)"
  eff <- function(par) {
    co <- fit$relist(par)
    vapply(rows, function(r) {
      k <- r[[3]]; j <- r[[4]]; g <- r[[5]]
      if (is.null(g)) {
        mean(roy_potential(fit, k, NULL, coef_list = co) - roy_potential(fit, j, NULL, coef_list = co))
      } else {
        i <- d == g
        mean((roy_potential(fit, k, g, coef_list = co) - roy_potential(fit, j, g, coef_list = co))[i])
      }
    }, 1)
  }
  est <- eff(fit$par)
  sev <- rep(NA_real_, length(est))
  if (se) {
    if (!is.null(fit$boot)) {
      bd <- apply(fit$boot, 1, eff)
      sev <- apply(matrix(bd, nrow = length(est)), 1, stats::sd)
    } else {
      ok <- is.finite(diag(fit$vcov))
      f2 <- function(p) { pp <- fit$par; pp[ok] <- p; eff(pp) }
      Jm <- numDeriv::jacobian(f2, fit$par[ok])
      sev <- sqrt(pmax(diag(Jm %*% fit$vcov[ok, ok] %*% t(Jm)), 0))
    }
  }
  data.frame(contrast = vapply(rows, `[[`, "", 1), population = vapply(rows, `[[`, "", 2),
             estimate = est, se = sev, stringsAsFactors = FALSE)
}
