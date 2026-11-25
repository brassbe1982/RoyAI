# ---------------------------------------------------------------------------
# roy_diagnostics: identification, convergence and specification checks
# ---------------------------------------------------------------------------

#' Diagnostics for a fitted Roy model
#'
#' Collects the checks that decide whether a Roy model estimate can be
#' trusted:
#' \describe{
#'   \item{convergence}{optimizer status, largest absolute score, smallest
#'     Hessian eigenvalue and its condition number (ML fits).}
#'   \item{selection}{Wald test of no selection on unobservables
#'     (\eqn{\rho_0 = \rho_1 = 0}, or all correction coefficients zero).}
#'   \item{instruments}{likelihood-ratio test that the excluded selection
#'     variables (exclusion restrictions) have no effect on selection; a weak
#'     result means identification leans on functional form.}
#'   \item{overlap}{distribution of fitted regime probabilities and the share
#'     of observations with probabilities near 0 or 1 (limited support for
#'     counterfactuals and for MTE estimation).}
#'   \item{boundary}{dependence parameters close to the edge of the
#'     parameter space, and (Model B) near-collinear correction terms.}
#' }
#'
#' @param fit A `roy_fit` object.
#' @param overlap_eps Threshold for "near 0 or 1" probabilities.
#' @param rho_bound Absolute correlation above which a boundary flag is raised.
#' @param kappa_bound Condition number above which second-stage designs are
#'   flagged as near-collinear.
#' @return An object of class `roy_diagnostics`: a list of the components
#'   above plus `flags`, a character vector of warnings in plain language.
#' @export
roy_diagnostics <- function(fit, overlap_eps = 0.02, rho_bound = 0.95, kappa_bound = 1e3) {
  stopifnot(inherits(fit, "roy_fit"))
  spec <- fit$spec; flags <- character(0)

  # convergence --------------------------------------------------------------
  conv <- list(converged = fit$converged)
  if (!is.null(fit$optim)) {
    ev <- eigen(fit$optim$hessian, symmetric = TRUE, only.values = TRUE)$values
    conv$max_abs_score <- max(abs(fit$optim$gradient))
    conv$min_hessian_eigen <- min(ev)
    conv$hessian_condition <- max(ev) / min(ev)
    if (min(ev) <= 0) flags <- c(flags, "Hessian is not positive definite: the likelihood is flat or the optimum is a saddle point; standard errors are unreliable.")
    if (conv$max_abs_score > 1e-2) flags <- c(flags, sprintf("largest score is %.3g: the optimizer may have stopped early; refit with tighter control or new starting values.", conv$max_abs_score))
  }
  if (!isTRUE(fit$converged)) flags <- c(flags, "the optimizer did not report convergence.")

  # selection on unobservables ----------------------------------------------
  sel_par <- if (fit$model == "A") {
    if (fit$method == "twostep") c("cov_sd0", "cov_sd1")
    else if (!is.null(fit$family) && fit$family$name == "gaussian") c("rho0", "rho1")
    else if (!is.null(fit$family) && fit$family$npar > 0) paste0(rep(fit$family$par_names, 2), rep(0:1, each = fit$family$npar))
    else character(0)
  } else grep(":(m[0-9]+|lee)$", names(fit$par), value = TRUE)
  selection <- if (length(sel_par) && (fit$model == "B" || fit$method == "twostep" ||
                                         identical(fit$family$name, "gaussian"))) {
    .wald(fit$par[sel_par], fit$vcov[sel_par, sel_par, drop = FALSE])
  } else NULL
  if (!is.null(selection)) selection$parameters <- sel_par
  if (!is.null(fit$family) && !is.null(fit$family$tau) && fit$family$npar > 0) {
    selection$kendall_tau <- c(`0` = fit$family$tau(fit$coef$dependence[["0"]]),
                               `1` = fit$family$tau(fit$coef$dependence[["1"]]))
  }

  # instruments --------------------------------------------------------------
  instruments <- NULL
  if (!length(spec$excluded)) {
    if (!identical(fit$type, "classical"))
      flags <- c(flags, "no exclusion restriction: identification of the selection correction rests on functional form only.")
  } else {
    keep <- setdiff(colnames(spec$Z), spec$excluded)
    Zr <- spec$Z[, keep, drop = FALSE]
    full <- if (spec$J == 2) .fit_probit(spec$d, spec$Z) else .fit_mnl(spec$d, spec$Z, spec$J)
    rest <- if (spec$J == 2) .fit_probit(spec$d, Zr) else .fit_mnl(spec$d, Zr, spec$J)
    lr <- 2 * (full$loglik - rest$loglik); df <- full$npar - rest$npar
    instruments <- list(excluded = spec$excluded, lr = lr, df = df,
                        p_value = stats::pchisq(lr, df, lower.tail = FALSE))
    if (lr / df < 10) flags <- c(flags, sprintf("exclusion restrictions are weak (LR/df = %.1f < 10).", lr / df))
  }

  # overlap ------------------------------------------------------------------
  P <- roy_propensity(fit)
  by_regime <- lapply(seq_len(spec$J), function(s) {
    p <- P[, s]
    c(min = min(p), q05 = unname(stats::quantile(p, 0.05)), median = stats::median(p),
      q95 = unname(stats::quantile(p, 0.95)), max = max(p))
  })
  names(by_regime) <- spec$regimes
  extreme <- rowSums(P < overlap_eps | P > 1 - overlap_eps) > 0
  overlap <- list(summary = do.call(rbind, by_regime), share_extreme = mean(extreme), eps = overlap_eps)
  if (mean(extreme) > 0.05)
    flags <- c(flags, sprintf("%.1f%% of observations have a regime probability below %.2f or above %.2f: counterfactuals there are extrapolations; consider trimming.",
                              100 * mean(extreme), overlap_eps, 1 - overlap_eps))

  # boundary -----------------------------------------------------------------
  boundary <- list()
  if (fit$model == "A" && !is.null(fit$family) && fit$family$name == "gaussian") {
    r <- fit$par[intersect(c("rho0", "rho1"), names(fit$par))]
    boundary$rho <- r
    if (any(abs(r) > rho_bound))
      flags <- c(flags, sprintf("a correlation is near the boundary (|rho| > %.2f): the likelihood may be poorly identified.", rho_bound))
  }
  if (fit$model == "B") {
    boundary$condition_number <- stats::setNames(fit$cond_W, spec$regimes)
    if (any(fit$cond_W > kappa_bound))
      flags <- c(flags, sprintf("second-stage design is near-collinear in regime(s) %s: correction terms are weakly identified (common with Dubin-McFadden terms); report effects with bootstrap SEs.",
                                paste(spec$regimes[fit$cond_W > kappa_bound], collapse = ", ")))
  }

  structure(list(convergence = conv, selection = selection, instruments = instruments,
                 overlap = overlap, boundary = boundary, flags = flags),
            class = "roy_diagnostics")
}

.wald <- function(b, V) {
  ok <- is.finite(b) & is.finite(diag(V))
  b <- b[ok]; V <- V[ok, ok, drop = FALSE]
  if (!length(b)) return(NULL)
  stat <- drop(t(b) %*% solve(V, b))
  list(estimate = b, stat = stat, df = length(b), p_value = stats::pchisq(stat, length(b), lower.tail = FALSE))
}

#' @export
print.roy_diagnostics <- function(x, ...) {
  cat("RoyAI diagnostics\n=================\n")
  cv <- x$convergence
  cat("Convergence:", cv$converged)
  if (!is.null(cv$max_abs_score))
    cat(sprintf(" | max |score| = %.2e | min Hessian eigenvalue = %.3g | condition = %.3g",
                cv$max_abs_score, cv$min_hessian_eigen, cv$hessian_condition))
  cat("\n")
  if (!is.null(x$selection))
    cat(sprintf("No selection on unobservables (H0: %s = 0): Wald = %.2f, df = %d, p = %.3g\n",
                paste(x$selection$parameters, collapse = ", "), x$selection$stat, x$selection$df,
                x$selection$p_value))
  if (!is.null(x$selection$kendall_tau))
    cat("Kendall's tau (U0,UD), (U1,UD):", paste(round(x$selection$kendall_tau, 3), collapse = ", "), "\n")
  if (!is.null(x$instruments))
    cat(sprintf("Exclusion restrictions (%s): LR = %.1f, df = %d, p = %.3g\n",
                paste(x$instruments$excluded, collapse = ", "), x$instruments$lr, x$instruments$df,
                x$instruments$p_value))
  cat(sprintf("Overlap: %.1f%% of observations with a regime probability outside [%.2f, %.2f]\n",
              100 * x$overlap$share_extreme, x$overlap$eps, 1 - x$overlap$eps))
  print(round(x$overlap$summary, 3))
  if (length(x$flags)) {
    cat("\nFlags:\n"); for (f in x$flags) cat(" -", f, "\n")
  } else cat("\nNo flags raised.\n")
  invisible(x)
}

#' Likelihood-ratio test between nested Roy models
#'
#' Typical use: the classical or extended Roy model against the generalized
#' Roy model, or the independence family against the Gaussian family.
#'
#' @param restricted,full Fitted `roy_fit` objects estimated by ML on the
#'   same data, with `restricted` nested in `full`.
#' @return A list of class `htest`.
#' @export
roy_lrtest <- function(restricted, full) {
  if (is.na(restricted$loglik) || is.na(full$loglik)) stop("both models must be ML fits")
  if (restricted$spec$n != full$spec$n) stop("models were fitted on different samples")
  stat <- 2 * (full$loglik - restricted$loglik); df <- full$npar - restricted$npar
  if (df <= 0) stop("'full' must have more parameters than 'restricted'")
  structure(list(statistic = c(LR = stat), parameter = c(df = df),
                 p.value = stats::pchisq(stat, df, lower.tail = FALSE),
                 method = "Likelihood-ratio test for nested Roy models",
                 data.name = sprintf("%s (restricted) vs %s (full)",
                                     .model_label(restricted), .model_label(full))),
            class = "htest")
}

.model_label <- function(f) {
  if (f$model == "A") sprintf("%s/%s", f$type, if (is.null(f$family)) "-" else f$family$name)
  else sprintf("B/%s", f$type)
}

#' Information criteria for competing Roy models
#'
#' @param ... Fitted `roy_fit` objects (ML), optionally named.
#' @return A data frame with log-likelihood, number of parameters, AIC and
#'   BIC, sorted by AIC.
#' @export
roy_ic <- function(...) {
  fits <- list(...)
  nm <- names(fits)
  if (is.null(nm)) nm <- vapply(fits, .model_label, "")
  nm[nm == ""] <- vapply(fits[nm == ""], .model_label, "")
  out <- data.frame(model = nm,
                    loglik = vapply(fits, function(f) f$loglik, 1),
                    npar = vapply(fits, function(f) f$npar, 1),
                    n = vapply(fits, function(f) f$spec$n, 1))
  out$AIC <- -2 * out$loglik + 2 * out$npar
  out$BIC <- -2 * out$loglik + log(out$n) * out$npar
  out[order(out$AIC), ]
}
