# =============================================================================
# roy_continuous.R  --  Model D: copula Roy model for continuous AI-use intensity
# -----------------------------------------------------------------------------
# Monograph Section 5.4 ("A copula Roy model for intensity"), Proposition
# prop:D-lik (likelihood) and Theorem thm:D-id (identification).
#
#   latent propensity  A* = Z'gamma + V,                V ~ N(0, 1)
#   observed intensity A  = 0 if A* <= 0, Lambda(A*) if A* > 0, A in [0, 1)
#   potential outcomes Y(a) = X'beta_x + D(a)'beta_a + sigma(a) eps,
#                      log sigma(a) = W(a)'alpha,        eps ~ N(0, 1)
#   (Phi(eps), Phi(V)) ~ C(. , . ; tau)   one of roy_copula_families
#
# Lambda is known: "exp" Lambda(t) = 1 - exp(-t) or "tanh" Lambda(t) = tanh(t/2).
# tau > 0 means workers with a strong pull towards AI (high V) also have high
# outcome shocks, so with sigma increasing in a they gain more from intensity.
# Rank invariance: the same eps enters Y(a) for every a.
# =============================================================================

.rc_links <- list(
  exp  = list(inv = function(a) -log1p(-a), ljac = function(a) -log1p(-a),
              fwd = function(t) -expm1(-t)),
  tanh = list(inv = function(a) 2 * atanh(a), ljac = function(a) log(2) - log1p(-a^2),
              fwd = function(t) tanh(t / 2)))
# ljac(a) = log d Lambda^{-1}(a) / da: exp: 1/(1-a);  tanh: 2/(1-a^2)

# ---- copula log-density (needed for users, A > 0) -----------------------------
# Base (unrotated) families; rotation by 90 degrees: c_rot(u, v) = c(1 - u, v),
# consistent with copula_hfunc().
.copula_logdens_base <- function(u, v, cp) {
  u <- .clamp01(u, 1e-12); v <- .clamp01(v, 1e-12); th <- cp$theta
  switch(cp$family,
    gaussian = {
      x <- stats::qnorm(u); y <- stats::qnorm(v); r <- th
      -0.5 * log1p(-r^2) - (r^2 * (x^2 + y^2) - 2 * r * x * y) / (2 * (1 - r^2))
    },
    t = {
      nu <- cp$nu; r <- th; x <- stats::qt(u, nu); y <- stats::qt(v, nu)
      lgamma((nu + 2) / 2) - lgamma(nu / 2) - log(nu * pi) - 0.5 * log1p(-r^2) -
        (nu + 2) / 2 * log1p((x^2 - 2 * r * x * y + y^2) / (nu * (1 - r^2))) -
        stats::dt(x, nu, log = TRUE) - stats::dt(y, nu, log = TRUE)
    },
    clayton = {
      if (th < 1e-10) return(numeric(length(u)) * 0)
      log1p(th) - (th + 1) * (log(u) + log(v)) -
        (1 / th + 2) * log(u^(-th) + v^(-th) - 1)
    },
    gumbel = {
      if (th - 1 < 1e-10) return(numeric(length(u)) * 0)
      x <- -log(u); y <- -log(v); S <- x^th + y^th; A <- S^(1 / th)
      -A - log(u) - log(v) + (th - 1) * (log(x) + log(y)) +
        (2 / th - 2) * log(S) + log(1 + (th - 1) / A)
    },
    frank = {
      if (abs(th) < 1e-8) return(numeric(length(u)) * 0)
      den <- expm1(-th) + expm1(-th * u) * expm1(-th * v)
      log(-th * expm1(-th)) - th * (u + v) - 2 * log(abs(den))
    })
}

#' Copula log-density
#'
#' Log of the copula density c(u, v) for the families of
#' [roy_copula_families], parameterised by Kendall's tau through
#' [copula_par()]. Rotated (negative-tau) Clayton and Gumbel copulas use
#' c(1 - u, v), consistent with [copula_hfunc()].
#' @param u outcome uniform(s)
#' @param v selection uniform(s)
#' @param cp copula parameter list from [copula_par()]
#' @return numeric vector of log-densities
#' @export
copula_logdens <- function(u, v, cp) {
  if (cp$rot == 90) .copula_logdens_base(1 - u, v, cp) else .copula_logdens_base(u, v, cp)
}

# ---- design matrices ------------------------------------------------------------
.rc_design <- function(sp, data, a = NULL) {
  if (!is.null(a)) data[[sp$aname]] <- a
  M <- stats::model.matrix(stats::delete.response(stats::terms(sp$outcome)), data)
  Da <- stats::model.matrix(sp$dose, data)
  Da <- Da[, colnames(Da) != "(Intercept)", drop = FALSE]
  W <- stats::model.matrix(sp$scale, data)
  list(X = M, D = Da, W = W)
}

.rc_layout <- function(sp) {
  c(paste0("beta:", sp$xnames), paste0("beta:", sp$dnames),
    paste0("gamma:", sp$znames), paste0("alpha:", sp$wnames), "atanh_tau",
    if (sp$family == "t" && is.null(sp$fix_nu)) "log_nu_m2")
}

.rc_unpack <- function(par, sp) {
  i <- 0; take <- function(n) { r <- par[i + seq_len(n)]; i <<- i + n; unname(r) }
  list(bx = take(length(sp$xnames)), bd = take(length(sp$dnames)),
       gamma = take(length(sp$znames)), alpha = take(length(sp$wnames)),
       tau = tanh(take(1)),
       nu = if (sp$family == "t") { if (is.null(sp$fix_nu)) 2 + exp(take(1)) else sp$fix_nu } else NULL)
}

.rc_loglik_i <- function(par, dat, sp) {
  th <- .rc_unpack(par, sp)
  cp <- copula_par(sp$family, th$tau, th$nu)
  zg <- drop(dat$Z %*% th$gamma)
  mu <- drop(dat$X %*% th$bx) + (if (length(th$bd)) drop(dat$D %*% th$bd) else 0)
  lsig <- drop(dat$W %*% th$alpha)
  e <- (dat$y - mu) / exp(lsig)
  u <- stats::pnorm(e)
  ll <- stats::dnorm(e, log = TRUE) - lsig
  z0 <- dat$a <= 0
  if (any(z0)) {
    h <- copula_hfunc(stats::pnorm(-zg[z0]), u[z0], cp)
    ll[z0] <- ll[z0] + log(pmax(h, 1e-300))
  }
  if (any(!z0)) {
    v <- dat$tA[!z0] - zg[!z0]
    ll[!z0] <- ll[!z0] + stats::dnorm(v, log = TRUE) +
      copula_logdens(u[!z0], stats::pnorm(v), cp) + dat$ljac[!z0]
  }
  ll
}

# Two-step control-function estimator: censored first stage, then OLS of y on
# the outcome design and the generalized residual E[V | A, Z] (Gaussian
# copula, homoskedastic). Used for starting values and as a benchmark.
.rc_twostep <- function(dat, sp) {
  pos <- dat$a > 0
  g0 <- stats::coef(stats::glm.fit(dat$Z, as.numeric(pos), family = stats::binomial("probit")))
  nll <- function(g) {
    zg <- drop(dat$Z %*% g)
    -sum(stats::dnorm(dat$tA[pos] - zg[pos], log = TRUE)) -
      sum(stats::pnorm(-zg[!pos], log.p = TRUE))
  }
  gamma <- stats::optim(g0, nll, method = "BFGS")$par
  zg <- drop(dat$Z %*% gamma)
  vhat <- numeric(length(dat$y))
  vhat[pos] <- dat$tA[pos] - zg[pos]
  vhat[!pos] <- -exp(stats::dnorm(zg[!pos], log = TRUE) - stats::pnorm(-zg[!pos], log.p = TRUE))
  R <- cbind(dat$X, dat$D, vhat)
  ols <- stats::lm.fit(R, dat$y)
  b <- ols$coefficients; b[!is.finite(b)] <- 0
  s <- sqrt(sum(ols$residuals^2) / (length(dat$y) - ncol(R)) + b[length(b)]^2)
  rho <- max(min(b[length(b)] / s, 0.9), -0.9)
  list(beta = b[-length(b)], gamma = gamma, sigma = s, rho = rho,
       tau = 2 / pi * asin(rho), vhat = vhat)
}

#' Fit Model D: copula Roy model for continuous AI-use intensity
#'
#' Maximum-likelihood estimation of the censored-intensity copula Roy model
#' of the monograph (Section 5.4, Proposition "Likelihood of Model D").
#' Intensity `A` in \[0, 1) is zero for non-users (latent propensity
#' `Z'gamma + V <= 0`) and `Lambda(Z'gamma + V)` for users, where `Lambda`
#' is a known link. Outcomes are `Y(a) = X'beta_x + D(a)'beta_a +
#' sigma(a) eps` with `log sigma(a) = W(a)'alpha` and the copula of
#' `(Phi(eps), Phi(V))` given by `family`, parameterised by Kendall's tau.
#'
#' The same safeguards as [roy_copula_fit()] apply: for non-Gaussian
#' families the likelihood is maximised from several starts (two-step
#' control function, independence, flipped tau) and the best *regular*
#' maximum is kept (positive definite Hessian, max |score| <= `score_tol`,
#' |tau| <= 0.98). The runs are stored in `fit$multistart`.
#'
#' @param outcome formula `y ~ x1 + ...` for the part of `mu(x, a)` that does
#'   not involve intensity (the intercept lives here)
#' @param intensity formula `a ~ z1 + ...`: intensity on the left, selection
#'   regressors (including at least one excluded variable) on the right
#' @param data data.frame
#' @param dose one-sided formula in the intensity variable for the dose part of
#'   `mu`, e.g. `~ a + I(a^2)`; intercept is dropped. Default: linear and
#'   quadratic. Avoid data-dependent bases such as `poly()`.
#' @param scale one-sided formula for `log sigma(a)`, e.g. `~ a` (default) or
#'   `~ 1` (homoskedastic)
#' @param family copula family, one of [roy_copula_families]
#' @param link `"exp"` (`Lambda(t) = 1 - exp(-t)`) or `"tanh"`
#'   (`Lambda(t) = tanh(t/2)`)
#' @param fix_nu fix the t-copula degrees of freedom (NULL = estimate)
#' @param start optional start (unconstrained scale, see `fit$par` names)
#' @param multistart `NULL` (default) means `TRUE` for non-Gaussian families
#'   when `start` is not given
#' @param score_tol largest absolute score accepted at a regular maximum
#' @param control passed to [stats::optim()]
#' @return object of class `"roy_continuous"`
#' @references Niankara, I. et al. (2025) <doi:10.1109/GCET68529.2025.11450624>.
#' @examples
#' set.seed(1)
#' d <- roy_continuous_sim(800, tau = 0.4)
#' fit <- roy_continuous(y ~ x1, a ~ x1 + z1 + z2, d, family = "gaussian")
#' fit
#' roy_continuous_effects(fit, a = 0.5, se = FALSE)
#' @export
roy_continuous <- function(outcome, intensity, data, dose = NULL, scale = NULL,
                           family = c("gaussian", "t", "clayton", "gumbel", "frank"),
                           link = c("exp", "tanh"), fix_nu = NULL, start = NULL,
                           multistart = NULL, score_tol = 1,
                           control = list(maxit = 3000, reltol = 1e-10)) {
  family <- match.arg(family); link <- match.arg(link)
  cl <- match.call()
  aname <- all.vars(intensity[[2]])
  if (length(aname) != 1) stop("the left-hand side of `intensity` must be one variable")
  if (is.null(dose)) dose <- stats::as.formula(paste0("~ ", aname, " + I(", aname, "^2)"))
  if (is.null(scale)) scale <- stats::as.formula(paste0("~ ", aname))
  data <- as.data.frame(data)
  vars <- unique(c(all.vars(outcome), all.vars(intensity)))
  data <- data[stats::complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  a <- data[[aname]]
  if (any(a < 0 | a >= 1)) stop("intensity must lie in [0, 1)")
  if (mean(a == 0) == 0 || mean(a > 0) == 0) stop("need both non-users (a = 0) and users (a > 0)")
  sp <- list(outcome = outcome, intensity = intensity, dose = dose, scale = scale,
             aname = aname, family = family, link = link, fix_nu = fix_nu)
  des <- .rc_design(sp, data)
  Z <- stats::model.matrix(intensity, data)
  sp$xnames <- colnames(des$X); sp$dnames <- colnames(des$D)
  sp$znames <- colnames(Z); sp$wnames <- colnames(des$W)
  lk <- .rc_links[[link]]
  pos <- a > 0
  tA <- numeric(length(a)); tA[pos] <- lk$inv(a[pos])
  ljac <- numeric(length(a)); ljac[pos] <- lk$ljac(a[pos])
  des0 <- .rc_design(sp, data, a = 0)
  # Non-users are evaluated at a = 0: mu(x, 0), sigma(0) -- the design already is.
  dat <- list(y = stats::model.response(stats::model.frame(outcome, data)),
              a = a, tA = tA, ljac = ljac, X = des$X, D = des$D, W = des$W, Z = Z)
  nms <- .rc_layout(sp)

  ts <- .rc_twostep(dat, sp)
  user_start <- !is.null(start)
  if (is.null(multistart)) multistart <- !user_start && family != "gaussian"
  if (!user_start) {
    start <- c(ts$beta, ts$gamma, log(ts$sigma), rep(0, length(sp$wnames) - 1),
               atanh(max(min(ts$tau, 0.9), -0.9)),
               if (family == "t" && is.null(fix_nu)) log(6))
  }
  if (length(start) != length(nms)) stop("`start` must have length ", length(nms))
  names(start) <- nms
  itau <- which(nms == "atanh_tau")
  starts <- list(twostep = start)
  if (isTRUE(multistart)) {
    st <- start; st[itau] <- 0; starts$independence <- st
    st <- start; st[itau] <- -start[itau]; starts$flipped_tau <- st
  }
  nll <- function(p) {
    v <- -sum(.rc_loglik_i(p, dat, sp))
    if (!is.finite(v)) 1e100 else v
  }
  run1 <- function(st) {
    opt <- stats::optim(st, nll, method = "BFGS", control = control)
    if (opt$convergence != 0) {
      o2 <- stats::nlminb(opt$par, nll, control = list(eval.max = 3000, iter.max = 2000))
      if (o2$objective <= opt$value) {
        opt$par <- o2$par; opt$value <- o2$objective; opt$convergence <- o2$convergence
      }
    }
    names(opt$par) <- nms
    opt
  }
  runs <- lapply(starts, run1)
  vals <- vapply(runs, function(o) o$value, numeric(1))
  ord <- order(vals); reg <- rep(NA, length(runs)); Hs <- vector("list", length(runs))
  maxsc <- rep(NA_real_, length(runs)); pick <- NA_integer_
  for (i in ord) {
    Hs[[i]] <- stats::optimHess(runs[[i]]$par, nll)
    ev <- tryCatch(eigen(Hs[[i]], symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    g <- tryCatch(numDeriv::grad(nll, runs[[i]]$par), error = function(e) NA_real_)
    maxsc[i] <- max(abs(g))
    reg[i] <- all(is.finite(ev)) && min(ev) > 0 && all(is.finite(g)) &&
      maxsc[i] <= score_tol && abs(tanh(runs[[i]]$par[itau])) <= 0.98
    if (reg[i]) { pick <- i; break }
  }
  if (is.na(pick)) {
    pick <- ord[1]
    warning("no start reached a regular maximum (positive definite Hessian, small score, ",
            "|tau| <= 0.98); the best run is returned, treat its standard errors with caution",
            call. = FALSE)
  } else if (pick != ord[1]) {
    warning(sprintf("the highest log-likelihood (start '%s') is irregular; the best regular maximum (start '%s') is used",
                    names(starts)[ord[1]], names(starts)[pick]), call. = FALSE)
  }
  opt <- runs[[pick]]
  H <- Hs[[pick]]
  V <- tryCatch(solve(H), error = function(e) .ginv(H)); dimnames(V) <- list(nms, nms)
  ms <- data.frame(start = names(starts), loglik = -unname(vals), regular = reg,
                   max_score = maxsc, chosen = seq_along(runs) == pick)
  structure(list(call = cl, sp = sp, family = family, link = link, par = opt$par,
                 vcov_raw = V, theta = .rc_unpack(opt$par, sp), loglik = -opt$value,
                 n = length(a), npar = length(opt$par), convergence = opt$convergence,
                 regular = isTRUE(reg[pick]), multistart = ms,
                 twostep = ts[c("beta", "gamma", "sigma", "rho", "tau")],
                 data = data, dat = dat, share_users = mean(pos)),
            class = "roy_continuous")
}

# Natural-scale parameter vector (tau and nu instead of atanh/log scale).
.rc_natural <- function(par, sp) {
  th <- .rc_unpack(par, sp)
  out <- c(th$bx, th$bd, th$gamma, th$alpha, th$tau, th$nu)
  names(out) <- c(paste0("beta:", sp$xnames), paste0("beta:", sp$dnames),
                  paste0("gamma:", sp$znames), paste0("alpha:", sp$wnames), "tau",
                  if (sp$family == "t") "nu")
  out
}

#' Methods for Model D fits
#'
#' @param object,x a `"roy_continuous"` fit
#' @param ... unused
#' @name roy_continuous_methods
NULL

#' @rdname roy_continuous_methods
#' @export
coef.roy_continuous <- function(object, ...) .rc_natural(object$par, object$sp)

#' @rdname roy_continuous_methods
#' @export
vcov.roy_continuous <- function(object, ...) {
  J <- numDeriv::jacobian(function(p) .rc_natural(p, object$sp), object$par)
  V <- J %*% object$vcov_raw %*% t(J)
  nm <- names(coef(object)); dimnames(V) <- list(nm, nm); V
}

#' @rdname roy_continuous_methods
#' @export
logLik.roy_continuous <- function(object, ...)
  structure(object$loglik, df = object$npar, nobs = object$n, class = "logLik")

#' @rdname roy_continuous_methods
#' @export
nobs.roy_continuous <- function(object, ...) object$n

#' @rdname roy_continuous_methods
#' @export
print.roy_continuous <- function(x, ...) {
  cat("Model D: copula Roy model for continuous AI-use intensity\n")
  cat(sprintf("family %s, link %s, n = %d (users %.1f%%), logLik %.2f, %s maximum\n",
              x$family, x$link, x$n, 100 * x$share_users, x$loglik,
              if (x$regular) "regular" else "IRREGULAR"))
  print(round(coef(x), 4))
  invisible(x)
}

#' @rdname roy_continuous_methods
#' @export
summary.roy_continuous <- function(object, ...) {
  b <- coef(object); se <- sqrt(pmax(diag(vcov(object)), 0))
  tab <- cbind(Estimate = b, `Std. Error` = se, `z value` = b / se,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(b / se)))
  structure(list(fit = object, coefficients = tab), class = "summary.roy_continuous")
}

#' @rdname roy_continuous_methods
#' @export
print.summary.roy_continuous <- function(x, ...) {
  print(x$fit)
  cat("\n"); stats::printCoefmat(x$coefficients, digits = 4)
  ms <- x$fit$multistart
  if (nrow(ms) > 1) { cat("\nMultistart runs:\n"); print(ms, row.names = FALSE) }
  invisible(x)
}

# ---- effects ------------------------------------------------------------------
# E[eps | V = v] for the fitted copula (Gaussian: rho * v).
.rc_m <- function(v, cp) {
  if (cp$family == "gaussian") return(cp$theta * v)
  uv <- sort(unique(round(v, 3)))
  m <- copula_cond_mean(stats::pnorm(uv), cp, 1)
  stats::approx(uv, m, xout = v, rule = 2)$y
}

.rc_effects_core <- function(par, sp, data, Z, a, delta, grid) {
  th <- .rc_unpack(par, sp)
  cp <- copula_par(sp$family, th$tau, th$nu)
  lk <- .rc_links[[sp$link]]
  zg <- drop(Z %*% th$gamma)
  mu_at <- function(aa) {
    d <- .rc_design(sp, data, a = rep(aa, nrow(data)))
    drop(d$X %*% th$bx) + (if (length(th$bd)) drop(d$D %*% th$bd) else 0)
  }
  sig_at <- function(aa) exp(drop(.rc_design(sp, data, a = rep(aa, nrow(data)))$W %*% th$alpha))
  dr <- vapply(grid, function(g) mean(mu_at(g)), numeric(1))
  tt <- vapply(a, function(aa) {
    v <- lk$inv(aa) - zg
    w <- stats::dnorm(v)                              # f(A = aa | z), up to the Jacobian
    g <- mu_at(min(aa + delta, 0.999)) - mu_at(aa) +
      (sig_at(min(aa + delta, 0.999)) - sig_at(aa)) * .rc_m(v, cp)
    sum(w * g) / sum(w)
  }, numeric(1))
  ate <- vapply(a, function(aa) mean(mu_at(min(aa + delta, 0.999)) - mu_at(aa)), numeric(1))
  # Non-users moved to intensity aa: E[Y(aa) - Y(0) | A = 0]
  p0 <- stats::pnorm(-zg)
  m0 <- -exp(stats::dnorm(zg, log = TRUE) - stats::pnorm(-zg, log.p = TRUE))  # Gaussian only
  if (cp$family != "gaussian") {
    # E[eps 1{V <= c}] = int_{-inf}^{c} m(v) phi(v) dv by quadrature
    vg <- seq(-8, 8, length.out = 321); mg <- .rc_m(vg, cp) * stats::dnorm(vg)
    Mg <- c(0, cumsum((mg[-1] + mg[-length(mg)]) / 2 * diff(vg)))
    m0 <- stats::approx(vg, Mg, xout = -zg, rule = 2)$y / pmax(p0, 1e-12)
  } else m0 <- cp$theta * m0
  atu <- vapply(a, function(aa) {
    sum(p0 * (mu_at(aa) - mu_at(0) + (sig_at(aa) - sig_at(0)) * m0)) / sum(p0)
  }, numeric(1))
  c(dr, tt, ate, atu)
}

#' Treatment effects for Model D
#'
#' Computes, from a [roy_continuous()] fit,
#' * the average dose-response function `mu_bar(a) = E[Y(a)]` on `grid`;
#' * `TT(a)`: the effect of raising intensity by `delta` for workers currently
#'   at intensity `a`, `E[Y(a + delta) - Y(a) | A = a]` (monograph Section 5.4);
#' * `ATE(a)`: the same increase averaged over all workers;
#' * `TUT(a)`: the effect of moving current non-users from 0 to `a`,
#'   `E[Y(a) - Y(0) | A = 0]`.
#' Conditioning on `A = a` integrates over the sample distribution of
#' `(X, Z)` weighted by the density of intensity at `a`. Delta-method
#' standard errors use the observed-information covariance.
#' @param fit a `"roy_continuous"` fit
#' @param a intensities at which to evaluate TT, ATE and TUT
#' @param delta increase in intensity
#' @param grid intensity grid for the dose-response function
#' @param se compute delta-method standard errors
#' @return data.frame with columns `estimand`, `a`, `estimate`, `se`
#' @export
roy_continuous_effects <- function(fit, a = c(0.2, 0.5, 0.8), delta = 0.1,
                                   grid = seq(0, 0.95, by = 0.05), se = TRUE) {
  stopifnot(inherits(fit, "roy_continuous"))
  f <- function(p) .rc_effects_core(p, fit$sp, fit$data, fit$dat$Z, a, delta, grid)
  est <- f(fit$par)
  s <- if (se) {
    J <- numDeriv::jacobian(f, fit$par, method = "simple")
    sqrt(pmax(diag(J %*% fit$vcov_raw %*% t(J)), 0))
  } else NA_real_
  data.frame(estimand = c(rep("dose_response", length(grid)), rep("TT", length(a)),
                          rep("ATE", length(a)), rep("TUT", length(a))),
             a = c(grid, a, a, a), estimate = est, se = s, row.names = NULL)
}

#' Simulate data from Model D
#'
#' Draws covariates `x1 ~ N(0,1)`, `z1 ~ N(0,1)`, `z2 ~ Bernoulli(0.5)`,
#' a latent propensity `Z'gamma + V` on `Z = (1, x1, z1, z2)`, intensity
#' `A = Lambda(Z'gamma + V)` when positive and 0 otherwise, and outcomes
#' `Y = b0 + b1 x1 + beta_a[1] a + beta_a[2] a^2 + exp(alpha[1] + alpha[2] a) eps`
#' with copula `family` and Kendall's `tau` between `Phi(eps)` and `Phi(V)`.
#' The default `gamma` gives about 30 percent non-users.
#' @param n sample size
#' @param beta_x outcome coefficients on (1, x1)
#' @param beta_a dose coefficients on (a, a^2); `beta_a[2] < 0` gives
#'   diminishing returns
#' @param alpha coefficients of `log sigma(a)` on (1, a)
#' @param gamma selection coefficients on (1, x1, z1, z2)
#' @param tau Kendall's tau
#' @param family copula family
#' @param nu t-copula degrees of freedom
#' @param link `"exp"` or `"tanh"`
#' @param covariates optional data.frame with x1, z1, z2
#' @return data.frame with `y`, `a`, covariates, `eps`, `v` and the
#'   potential-outcome function is recoverable from the true parameters
#' @export
roy_continuous_sim <- function(n, beta_x = c(1, 0.5), beta_a = c(0.8, -0.4),
                               alpha = c(log(0.6), 0.5), gamma = c(0.6, 0.3, 0.8, 0.6),
                               tau = 0.4, family = "gaussian", nu = 5,
                               link = c("exp", "tanh"), covariates = NULL) {
  link <- match.arg(link)
  if (is.null(covariates))
    covariates <- data.frame(x1 = stats::rnorm(n), z1 = stats::rnorm(n),
                             z2 = stats::rbinom(n, 1, 0.5))
  n <- nrow(covariates)
  Z <- cbind(1, as.matrix(covariates[, c("x1", "z1", "z2")]))
  w <- stats::runif(n); v <- stats::qnorm(w)
  cp <- copula_par(family, tau, if (family == "t") nu)
  u <- .hinv2(stats::runif(n), w, cp)
  eps <- stats::qnorm(u)
  ast <- drop(Z %*% gamma) + v
  a <- ifelse(ast > 0, .rc_links[[link]]$fwd(pmax(ast, 0)), 0)
  a <- pmin(a, 1 - 1e-9)
  mu <- beta_x[1] + beta_x[2] * covariates$x1 + beta_a[1] * a + beta_a[2] * a^2
  y <- mu + exp(alpha[1] + alpha[2] * a) * eps
  data.frame(y = y, a = a, covariates, eps = eps, v = v)
}

#' True Model D effects for a given parameter vector
#'
#' Evaluates the estimands of [roy_continuous_effects()] at known parameters
#' of the [roy_continuous_sim()] design (dose `~ a + I(a^2)`, scale `~ a`),
#' for Monte Carlo studies.
#' @inheritParams roy_continuous_sim
#' @inheritParams roy_continuous_effects
#' @param covariates data.frame with x1, z1, z2 over which to average
#' @export
roy_continuous_truth <- function(beta_x = c(1, 0.5), beta_a = c(0.8, -0.4),
                                 alpha = c(log(0.6), 0.5), gamma = c(0.6, 0.3, 0.8, 0.6),
                                 tau = 0.4, family = "gaussian", nu = 5, link = "exp",
                                 covariates, a = c(0.2, 0.5, 0.8), delta = 0.1,
                                 grid = seq(0, 0.95, by = 0.05)) {
  sp <- list(outcome = y ~ x1, dose = ~ a + I(a^2), scale = ~ a, aname = "a",
             family = family, link = link, fix_nu = if (family == "t") nu,
             xnames = c("(Intercept)", "x1"), dnames = c("a", "I(a^2)"),
             znames = c("(Intercept)", "x1", "z1", "z2"), wnames = c("(Intercept)", "a"))
  dd <- covariates; dd$a <- 0; dd$y <- 0
  Z <- cbind(1, as.matrix(covariates[, c("x1", "z1", "z2")]))
  par <- c(beta_x, beta_a, gamma, alpha, atanh(tau))
  est <- .rc_effects_core(par, sp, dd, Z, a, delta, grid)
  data.frame(estimand = c(rep("dose_response", length(grid)), rep("TT", length(a)),
                          rep("ATE", length(a)), rep("TUT", length(a))),
             a = c(grid, a, a, a), truth = est, row.names = NULL)
}
