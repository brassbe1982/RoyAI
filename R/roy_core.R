# ---------------------------------------------------------------------------
# roy_core: classical and generalized Roy models (Model A and Model B)
# ---------------------------------------------------------------------------

#' Fit a generalized Roy model
#'
#' One entry point for the binary AI-adoption Roy model (Model A) and the
#' multinomial Roy model of human-only, AI-augmented and AI-dominant
#' production (Model B). The model is chosen from the number of values of the
#' selection variable.
#'
#' **Model A.** \eqn{Y_s = X'\beta_s + U_s} (s = 0, 1),
#' \eqn{D = 1(Z'\gamma + U_D > 0)}, \eqn{U_D \sim N(0,1)},
#' \eqn{U_s \sim N(0, \sigma_s^2)}, with each pair \eqn{(U_s, U_D)} joined by
#' a dependence family (Gaussian by default, which gives the normal
#' endogenous-switching / generalized Roy model). `type` restricts the
#' selection rule:
#' \describe{
#'   \item{`"generalized"`}{selection index unrestricted (Heckman-Vytlacil
#'     generalized Roy model).}
#'   \item{`"extended"`}{\eqn{D = 1(Y_1 - Y_0 - W'\pi > 0)}: selection on
#'     gains net of an observed cost \eqn{W'\pi}, where \eqn{W} are the
#'     selection variables excluded from the outcome.}
#'   \item{`"classical"`}{Roy (1951): \eqn{D = 1(Y_1 > Y_0)}.}
#' }
#' The extended and classical models are parameter restrictions on the
#' generalized model, so [roy_lrtest()] tests them directly.
#'
#' **Model B.** \eqn{Y_s = X'\beta_s + U_s}, s = 0, 1, 2;
#' \eqn{S = \arg\max_s\{Z'\gamma_s + \varepsilon_s\}} with Gumbel
#' \eqn{\varepsilon_s}. Estimated in two steps: multinomial logit, then
#' regime-specific regressions with Dubin-McFadden terms in the
#' Bourguignon-Fournier-Gurgand (2007) form (`correction = "dmf"`, which also
#' identifies counterfactual means) or Lee's (1983) term (`"lee"`).
#'
#' @param outcome,selection,data See [roy_spec()].
#' @param method `"ml"` (full information maximum likelihood, Model A) or
#'   `"twostep"` (Heckman / multinomial correction). Methods registered with
#'   [roy_register_method()] are also accepted.
#' @param dependence Name of a registered dependence family or a
#'   [roy_dependence()] object (Model A, `method = "ml"`).
#' @param type `"generalized"`, `"extended"` or `"classical"` (Model A).
#' @param correction `"dmf"` or `"lee"` (Model B).
#' @param se `"analytic"` (inverse Hessian for ML; Murphy-Topel corrected
#'   sandwich for two-step) or `"bootstrap"`.
#' @param B Bootstrap replications when `se = "bootstrap"`.
#' @param control List of options passed to [stats::optim()].
#' @param start Model A ML only. Either a full parameter vector on the working
#'   scale (as in `fit$optim$par`), or a list of two natural-scale dependence
#'   parameter vectors (regimes 0 and 1), for example `list(-0.5, -0.4)` for a
#'   Kendall-tau family, combined with the two-step coefficients.
#' @param multistart Model A ML with a non-Gaussian dependence family. If
#'   `TRUE` (the default for such families), the likelihood is maximised from
#'   the family's default start (independence), from the Kendall tau implied by
#'   the two-step correlations, \eqn{2\arcsin(\rho)/\pi}, and from its sign
#'   flip, and the best maximum is kept. Copula likelihoods are often bimodal
#'   in the sign of the dependence. Runs that end at an irregular point (Hessian
#'   not positive definite, large score, or |dependence| > 0.98) are rejected
#'   in favour of the best regular maximum, with a warning. The runs are stored
#'   in `fit$multistart` (columns `loglik`, `regular`, `chosen`).
#' @param ... Passed to registered plug-in fitters.
#' @return An object of class `c("roy_fit_A", "roy_fit")` or
#'   `c("roy_fit_B", "roy_fit")`.
#' @examples
#' dat <- roy_sim_A(n = 2000, seed = 1)
#' fit <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
#' summary(fit)
#' roy_effects(fit)
#' @export
roy <- function(outcome, selection, data,
                method = c("ml", "twostep"), dependence = "gaussian",
                type = c("generalized", "extended", "classical"),
                correction = c("dmf", "lee"), se = c("analytic", "bootstrap"),
                B = 199, control = list(), start = NULL, multistart = NULL, ...) {
  cl <- match.call()
  # the classical Roy model has no exclusion restriction by construction
  quiet <- identical(type[1], "classical")
  spec <- withCallingHandlers(roy_spec(outcome, selection, data), warning = function(w)
    if (quiet && grepl("exclusion restriction", conditionMessage(w))) invokeRestart("muffleWarning"))
  if (length(method) > 1) method <- if (spec$model == "A") "ml" else "twostep"
  se <- match.arg(se)
  key <- paste(spec$model, method, sep = ":")
  if (!is.null(.roy_registry$methods[[key]])) {
    fit <- .roy_registry$methods[[key]](spec, dependence = dependence, se = se, B = B,
                                        control = control, ...)
  } else if (spec$model == "A") {
    type <- match.arg(type)
    fam <- roy_get_dependence(dependence)
    fit <- switch(method,
      ml = .fit_A_ml(spec, fam, type, control, start = start, multistart = multistart),
      twostep = {
        if (type != "generalized") stop("type = '", type, "' is only available with method = 'ml'")
        .fit_A_twostep(spec)
      },
      stop("unknown method '", method, "' for Model A"))
  } else {
    correction <- match.arg(correction)
    if (method != "twostep")
      stop("Model B is estimated by method = 'twostep' in the core; full-likelihood ",
           "versions are provided by plug-ins (registered: ",
           paste(roy_methods(), collapse = ", "), ")")
    fit <- .fit_B_twostep(spec, correction)
  }
  fit$call <- cl
  if (se == "bootstrap") fit <- .bootstrap_se(fit, B, method, dependence,
                                              if (spec$model == "A") match.arg(type) else NULL,
                                              if (spec$model == "B") match.arg(correction) else NULL,
                                              control)
  fit
}

# ---------------------------------------------------------------------------
# Model A: full-information maximum likelihood
# ---------------------------------------------------------------------------

# Map from the classical / extended Roy parameters to the generalized model.
# theta = (beta0, beta1, log s0, log s1, atanh rho01, pi)
.classical_map <- function(b0, b1, s0, s1, r01, pi_w, xcols, zcols) {
  sD <- sqrt(s0^2 + s1^2 - 2 * r01 * s0 * s1)
  g <- stats::setNames(numeric(length(zcols)), zcols)
  g[xcols] <- (b1 - b0) / sD
  if (length(pi_w)) g[names(pi_w)] <- -pi_w / sD
  list(gamma = g, rho0 = (r01 * s0 * s1 - s0^2) / (s0 * sD),
       rho1 = (s1^2 - r01 * s0 * s1) / (s1 * sD), sigmaD = sD)
}

.unpack_A <- function(th, spec, fam, type) {
  kx <- ncol(spec$X); kz <- ncol(spec$Z); np <- fam$npar
  if (type == "generalized") {
    i <- 0
    g <- th[i + 1:kz]; i <- i + kz
    b0 <- th[i + 1:kx]; i <- i + kx
    b1 <- th[i + 1:kx]; i <- i + kx
    s0 <- exp(unname(th[i + 1])); s1 <- exp(unname(th[i + 2])); i <- i + 2
    d0 <- fam$linkinv(unname(th[i + seq_len(np)])); i <- i + np
    d1 <- fam$linkinv(unname(th[i + seq_len(np)]))
    names(g) <- colnames(spec$Z); names(b0) <- names(b1) <- colnames(spec$X)
    list(gamma = g, beta0 = b0, beta1 = b1, sigma0 = s0, sigma1 = s1, dep0 = d0, dep1 = d1)
  } else {
    if (!all(colnames(spec$X) %in% colnames(spec$Z)))
      stop("classical/extended Roy needs every outcome regressor in the selection equation")
    b0 <- th[1:kx]; b1 <- th[kx + 1:kx]
    names(b0) <- names(b1) <- colnames(spec$X)
    s0 <- exp(unname(th[2 * kx + 1])); s1 <- exp(unname(th[2 * kx + 2])); r01 <- tanh(unname(th[2 * kx + 3]))
    pi_w <- if (type == "extended" && length(spec$excluded))
      stats::setNames(th[2 * kx + 3 + seq_along(spec$excluded)], spec$excluded) else numeric(0)
    m <- .classical_map(b0, b1, s0, s1, r01, pi_w, colnames(spec$X), colnames(spec$Z))
    list(gamma = m$gamma, beta0 = b0, beta1 = b1, sigma0 = s0, sigma1 = s1,
         dep0 = m$rho0, dep1 = m$rho1, rho01 = r01, pi = pi_w, sigmaD = m$sigmaD)
  }
}

# log-likelihood contributions of Model A for unpacked parameters p
.ll_A_obs <- function(p, spec, fam) {
  y <- spec$y; d <- spec$d
  cidx <- drop(spec$Z %*% p$gamma)
  e0 <- (y - drop(spec$X %*% p$beta0)) / p$sigma0
  e1 <- (y - drop(spec$X %*% p$beta1)) / p$sigma1
  ll <- numeric(length(y)); i1 <- d == 1; i0 <- !i1
  ll[i1] <- stats::dnorm(e1[i1], log = TRUE) - log(p$sigma1)
  ll[i0] <- stats::dnorm(e0[i0], log = TRUE) - log(p$sigma0)
  if (fam$name == "gaussian") {
    r0 <- p$dep0; r1 <- p$dep1
    if (!is.finite(r0) || !is.finite(r1) || abs(r0) >= 1 || abs(r1) >= 1) return(rep(-Inf, length(y)))
    ll[i1] <- ll[i1] + stats::pnorm((cidx[i1] + r1 * e1[i1]) / sqrt(1 - r1^2), log.p = TRUE)
    ll[i0] <- ll[i0] + stats::pnorm(-(cidx[i0] + r0 * e0[i0]) / sqrt(1 - r0^2), log.p = TRUE)
  } else {
    eps <- 1e-12
    clamp <- function(u) pmin(pmax(u, eps), 1 - eps)
    v <- clamp(stats::pnorm(-cidx))
    ll[i1] <- ll[i1] + log(pmax(1 - fam$h(v[i1], clamp(stats::pnorm(e1[i1])), p$dep1), 1e-300))
    ll[i0] <- ll[i0] + log(pmax(fam$h(v[i0], clamp(stats::pnorm(e0[i0])), p$dep0), 1e-300))
  }
  ll
}

.fit_A_ml <- function(spec, fam, type, control, start = NULL, multistart = NULL) {
  if (type != "generalized" && fam$name != "gaussian")
    stop("the classical and extended Roy restrictions are implemented for the Gaussian family")
  ts <- .fit_A_twostep(spec)
  clip <- function(r) pmin(pmax(r, -0.9), 0.9)
  st0 <- ts$coef
  starts <- list()
  if (is.numeric(start)) {
    starts$user <- start
  } else if (type == "generalized") {
    base <- c(st0$selection, st0$outcome[["0"]], st0$outcome[["1"]], log(st0$sigma))
    rho_ts <- clip(c(st0$dependence[["0"]], st0$dependence[["1"]]))
    dep <- function(d0, d1) c(base, fam$link(d0), fam$link(d1))
    if (fam$npar == 0) {
      starts$twostep <- base
    } else if (fam$name == "gaussian") {
      starts$twostep <- c(base, atanh(rho_ts))
    } else {
      if (is.list(start)) starts$user <- dep(start[[1]], start[[2]])
      if (is.null(multistart)) multistart <- TRUE
      starts$default <- dep(fam$start, fam$start)
      if (isTRUE(multistart) && identical(fam$par_names[1], "tau")) {
        # Kendall's tau implied by the two-step correlations (Gaussian map),
        # and its mirror image: copula likelihoods often have one mode per sign
        tau_ts <- 2 / pi * asin(rho_ts)
        mk <- function(tau) { th <- fam$start; th[1] <- tau; th }
        starts$twostep_tau <- dep(mk(tau_ts[1]), mk(tau_ts[2]))
        starts$flipped_tau <- dep(mk(-tau_ts[1]), mk(-tau_ts[2]))
      }
    }
  } else {
    starts$twostep <- c(st0$outcome[["0"]], st0$outcome[["1"]], log(st0$sigma), 0,
                        if (type == "extended") rep(0, length(spec$excluded)))
  }
  nll <- function(th) {
    p <- .unpack_A(th, spec, fam, type)
    v <- -sum(.ll_A_obs(p, spec, fam))
    if (!is.finite(v)) 1e100 else v
  }
  ctrl <- utils::modifyList(list(maxit = 2000, reltol = 1e-12), control)
  run1 <- function(st) {
    o <- stats::optim(st, nll, method = "BFGS", control = ctrl)
    if (o$convergence != 0) o <- stats::optim(o$par, nll, method = "BFGS", control = ctrl)
    o
  }
  runs <- lapply(starts, run1)
  vals <- vapply(runs, function(o) o$value, numeric(1))

  # Newton polish with step halving: BFGS on numerical gradients can stop
  # with a sizeable score when regressors are on very different scales.
  polish <- function(o) {
    grad <- numDeriv::grad(nll, o$par)
    for (it in seq_len(50)) {
      if (max(abs(grad)) < 1e-4) break
      H <- numDeriv::hessian(nll, o$par)
      step <- tryCatch(solve(H, grad), error = function(e) MASS_ginv(H) %*% grad)
      if (any(!is.finite(step))) break
      t <- 1; f0 <- o$value; improved <- FALSE
      while (t > 1e-6) {
        cand <- o$par - t * drop(step); fc <- nll(cand)
        if (is.finite(fc) && fc < f0) { o$par <- cand; o$value <- fc; improved <- TRUE; break }
        t <- t / 2
      }
      if (!improved) break
      grad <- numDeriv::grad(nll, o$par)
    }
    o$gradient <- grad
    o$hessian <- numDeriv::hessian(nll, o$par)
    o
  }
  # A run is a regular interior maximum if the Hessian of the negative
  # log-likelihood is positive definite, the score is small and the
  # dependence parameters are away from their bounds. Copula likelihoods can
  # spike as |tau| -> 1 (the density degenerates), giving a higher but
  # spurious log-likelihood with a non-PD Hessian; such runs are rejected.
  regular <- function(o) {
    ev <- tryCatch(eigen(o$hessian, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    if (!all(is.finite(ev)) || min(ev) <= 0 || max(abs(o$gradient)) > 1) return(FALSE)
    if (fam$npar > 0 && type == "generalized") {
      dp <- .unpack_A(o$par, spec, fam, type)
      dd <- c(dp$dep0[1], dp$dep1[1])
      if (length(dd) && any(abs(dd) > 0.98)) return(FALSE)
    }
    TRUE
  }
  ord <- order(vals)
  reg <- rep(NA, length(runs)); polished <- vector("list", length(runs))
  pick <- NA_integer_
  for (i in ord) {
    polished[[i]] <- polish(runs[[i]])
    reg[i] <- regular(polished[[i]])
    if (reg[i]) { pick <- i; break }
    if (length(runs) == 1) break
  }
  if (is.na(pick)) {
    pick <- ord[1]
    if (is.null(polished[[pick]])) polished[[pick]] <- polish(runs[[pick]])
    if (length(runs) > 1)
      warning("no start reached a regular maximum (positive definite Hessian, dependence ",
              "away from its bounds); the best run is returned, treat its standard errors ",
              "with caution", call. = FALSE)
  } else if (pick != ord[1]) {
    warning(sprintf(paste("the highest log-likelihood (start '%s') is an irregular boundary point;",
                          "the best regular maximum (start '%s') is used"),
                    names(starts)[ord[1]], names(starts)[pick]), call. = FALSE)
  }
  o <- polished[[pick]]; grad <- o$gradient; H <- o$hessian
  Vth <- tryCatch(solve(H), error = function(e) MASS_ginv(H))

  # natural-scale parameter vector and delta-method vcov
  nat <- function(th) .natural_A(.unpack_A(th, spec, fam, type), fam, type)
  par <- nat(o$par)
  Jm <- numDeriv::jacobian(nat, o$par)
  V <- Jm %*% Vth %*% t(Jm); dimnames(V) <- list(names(par), names(par))

  fit <- .new_fit_A(spec, fam, type, method = "ml", par = par, vcov = V,
                    loglik = -o$value, converged = o$convergence == 0,
                    npar = length(o$par))
  fit$optim <- list(par = o$par, value = o$value, counts = o$counts,
                    convergence = o$convergence, gradient = grad, hessian = H)
  fit$first_stage <- ts$first_stage
  fit$multistart <- data.frame(start = names(starts), loglik = -unname(vals),
                               regular = reg, chosen = seq_along(vals) == pick)
  fit
}

# natural-scale named parameter vector for Model A
.natural_A <- function(p, fam, type) {
  dn <- function(x, s) if (length(x)) stats::setNames(x, paste0(fam$par_names, s)) else NULL
  v <- c(stats::setNames(p$gamma, paste0("sel:", names(p$gamma))),
         stats::setNames(p$beta0, paste0("y0:", names(p$beta0))),
         stats::setNames(p$beta1, paste0("y1:", names(p$beta1))),
         sigma0 = p$sigma0, sigma1 = p$sigma1)
  if (type == "generalized") {
    v <- c(v, dn(p$dep0, "0"), dn(p$dep1, "1"))
  } else {
    v <- c(v, rho0 = p$dep0, rho1 = p$dep1, rho01 = p$rho01, sigmaD = p$sigmaD)
    if (length(p$pi)) v <- c(v, stats::setNames(p$pi, paste0("cost:", names(p$pi))))
  }
  v
}

# rebuild a coefficient list from a natural parameter vector
.relist_A <- function(par, spec, fam) {
  g <- par[paste0("sel:", colnames(spec$Z))]; names(g) <- colnames(spec$Z)
  b0 <- par[paste0("y0:", colnames(spec$X))]; names(b0) <- colnames(spec$X)
  b1 <- par[paste0("y1:", colnames(spec$X))]; names(b1) <- colnames(spec$X)
  pn <- if (fam$name == "gaussian") "rho" else fam$par_names
  d0 <- if (fam$npar) unname(par[paste0(pn, "0")]) else numeric(0)
  d1 <- if (fam$npar) unname(par[paste0(pn, "1")]) else numeric(0)
  sig <- c(`0` = unname(par["sigma0"]), `1` = unname(par["sigma1"]))
  cs <- if ("cov_sd0" %in% names(par)) c(`0` = unname(par["cov_sd0"]), `1` = unname(par["cov_sd1"]))
        else if (fam$name == "gaussian") sig * c(d0, d1)
        else if (fam$name == "independence") c(`0` = 0, `1` = 0) else NULL
  list(selection = g, outcome = list(`0` = b0, `1` = b1), sigma = sig,
       dependence = list(`0` = d0, `1` = d1), cov_sd = cs)
}

.new_fit_A <- function(spec, fam, type, method, par, vcov, loglik, converged, npar) {
  structure(list(spec = spec, model = "A", method = method, type = type, family = fam,
                 coef = .relist_A(par, spec, fam), par = par, vcov = vcov,
                 loglik = loglik, converged = converged, npar = npar,
                 relist = function(p) .relist_A(p, spec, fam), se_type = "analytic"),
            class = c("roy_fit_A", "roy_fit"))
}

MASS_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X); pos <- s$d > max(tol * s$d[1L], 0)
  s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

# ---------------------------------------------------------------------------
# Generic two-step second stage with Murphy-Topel corrected sandwich
# ---------------------------------------------------------------------------
# For regime s: y = X b + C(gamma) delta + e on the subsample S = s.
# V(b) = (W'W)^-1 (sum e^2 w w') (W'W)^-1 + A Vg A',  A = -(W'W)^-1 W' G,
# G = d(C(gamma) delta)/d gamma'.  Cov(b, gamma) = A Vg.
.second_stage <- function(y, X, Cfun, gamma, idx) {
  C <- Cfun(gamma)
  W <- cbind(X, C)[idx, , drop = FALSE]
  yy <- y[idx]
  qrW <- qr(W)
  if (qrW$rank < ncol(W)) stop("second-stage design is rank deficient (weak or collinear correction terms)")
  b <- qr.coef(qrW, yy); names(b) <- colnames(W)
  e <- drop(yy - W %*% b)
  WWi <- chol2inv(qrW$qr[seq_len(ncol(W)), seq_len(ncol(W)), drop = FALSE])
  WWi <- WWi[order(qrW$pivot), order(qrW$pivot)]
  meat <- crossprod(W * e)
  kx <- ncol(X)
  delta <- b[-(1:kx)]
  G <- numDeriv::jacobian(function(g) drop(Cfun(g)[idx, , drop = FALSE] %*% delta), gamma)
  A <- -WWi %*% crossprod(W, G)
  list(b = b, e = e, W = W, WWi = WWi, Vhc = WWi %*% meat %*% WWi, A = A, n = length(idx))
}

.assemble_vcov <- function(Vg, stages) {
  # joint vcov of (gamma, b_0, b_1, ...)
  kg <- nrow(Vg); ks <- vapply(stages, function(s) length(s$b), 1L)
  K <- kg + sum(ks); V <- matrix(0, K, K)
  V[1:kg, 1:kg] <- Vg
  off <- cumsum(c(kg, ks))
  for (a in seq_along(stages)) {
    ia <- off[a] + seq_len(ks[a]); Aa <- stages[[a]]$A
    V[ia, 1:kg] <- Aa %*% Vg; V[1:kg, ia] <- t(V[ia, 1:kg])
    for (b in seq_along(stages)) {
      ib <- off[b] + seq_len(ks[b])
      V[ia, ib] <- Aa %*% Vg %*% t(stages[[b]]$A) + if (a == b) stages[[a]]$Vhc else 0
    }
  }
  V
}

# ---------------------------------------------------------------------------
# Model A: Heckman two-step
# ---------------------------------------------------------------------------
.fit_A_twostep <- function(spec) {
  fs <- .fit_probit(spec$d, spec$Z)
  st <- lapply(0:1, function(s) {
    Cfun <- function(g) { cc <- drop(spec$Z %*% g)
      matrix(roy_imr(cc, rep(s, length(cc))), ncol = 1, dimnames = list(NULL, "lambda")) }
    .second_stage(spec$y, spec$X, Cfun, fs$coef, which(spec$d == s))
  })
  cidx <- fs$index
  sig <- rho <- numeric(2)
  for (s in 0:1) {
    idx <- spec$d == s
    lam <- roy_imr(cidx[idx], rep(s, sum(idx)))
    bl <- st[[s + 1]]$b["lambda"]
    s2 <- mean(st[[s + 1]]$e^2) + bl^2 * mean(lam * (lam + cidx[idx]))
    sig[s + 1] <- sqrt(max(s2, 1e-8)); rho[s + 1] <- max(min(bl / sig[s + 1], 0.999), -0.999)
  }
  kx <- ncol(spec$X)
  V <- .assemble_vcov(fs$vcov, st)
  par <- c(stats::setNames(fs$coef, paste0("sel:", colnames(spec$Z))),
           stats::setNames(st[[1]]$b[1:kx], paste0("y0:", colnames(spec$X))),
           cov_sd0 = unname(st[[1]]$b["lambda"]),
           stats::setNames(st[[2]]$b[1:kx], paste0("y1:", colnames(spec$X))),
           cov_sd1 = unname(st[[2]]$b["lambda"]))
  dimnames(V) <- list(names(par), names(par))
  # sigma and rho are derived (no standard errors in the two-step)
  par_all <- c(par, sigma0 = sig[1], sigma1 = sig[2], rho0 = rho[1], rho1 = rho[2])
  Vall <- matrix(NA_real_, length(par_all), length(par_all), dimnames = list(names(par_all), names(par_all)))
  Vall[names(par), names(par)] <- V
  fam <- roy_dep_gaussian()
  fit <- .new_fit_A(spec, fam, "generalized", "twostep", par_all, Vall,
                    loglik = NA_real_, converged = fs$converged, npar = length(par))
  fit$first_stage <- fs
  fit
}

# ---------------------------------------------------------------------------
# Model B: multinomial selection correction
# ---------------------------------------------------------------------------
.fit_B_twostep <- function(spec, correction) {
  fs <- .fit_mnl(spec$d, spec$Z, spec$J)
  J <- spec$J; kz <- ncol(spec$Z); kx <- ncol(spec$X)
  st <- lapply(0:(J - 1), function(s) {
    Cfun <- function(g) {
      P <- .mnl_prob(matrix(g, kz, J - 1), spec$Z)
      if (correction == "dmf") roy_dmf_terms(P, s)
      else matrix(roy_lee_terms(P, s), ncol = 1, dimnames = list(NULL, "lee"))
    }
    .second_stage(spec$y, spec$X, Cfun, as.vector(fs$coef), which(spec$d == s))
  })
  V <- .assemble_vcov(fs$vcov, st)
  gnm <- paste0("sel", rep(1:(J - 1), each = kz), ":", rep(colnames(spec$Z), J - 1))
  par <- stats::setNames(as.vector(fs$coef), gnm)
  for (s in 0:(J - 1)) {
    b <- st[[s + 1]]$b
    par <- c(par, stats::setNames(b, paste0("y", s, ":", names(b))))
  }
  dimnames(V) <- list(names(par), names(par))
  relist <- function(p) .relist_B(p, spec, correction)
  fit <- structure(list(spec = spec, model = "B", method = "twostep", type = correction,
                        family = NULL, par = par, vcov = V, coef = relist(par),
                        relist = relist, loglik = NA_real_, converged = fs$converged,
                        npar = length(par), first_stage = fs, se_type = "analytic",
                        sigma_resid = vapply(st, function(s) sqrt(mean(s$e^2)), 1),
                        cond_W = vapply(st, function(s) kappa(s$W, exact = FALSE), 1)),
                   class = c("roy_fit_B", "roy_fit"))
  fit
}

.relist_B <- function(p, spec, correction) {
  J <- spec$J; kz <- ncol(spec$Z)
  G <- matrix(p[seq_len(kz * (J - 1))], kz, J - 1,
              dimnames = list(colnames(spec$Z), paste0("alt", 1:(J - 1))))
  out <- lam <- list()
  for (s in 0:(J - 1)) {
    pre <- paste0("y", s, ":")
    b <- p[startsWith(names(p), pre)]; names(b) <- sub(pre, "", names(b), fixed = TRUE)
    out[[as.character(s)]] <- b[colnames(spec$X)]
    lam[[as.character(s)]] <- b[setdiff(names(b), colnames(spec$X))]
  }
  list(selection = G, outcome = out, correction = lam)
}

# ---------------------------------------------------------------------------
# Bootstrap standard errors (nonparametric, resampling observations)
# ---------------------------------------------------------------------------
.bootstrap_se <- function(fit, B, method, dependence, type, correction, control) {
  dat <- fit$spec$data; n <- nrow(dat)
  draws <- matrix(NA_real_, B, length(fit$par), dimnames = list(NULL, names(fit$par)))
  for (b in seq_len(B)) {
    db <- dat[sample.int(n, n, replace = TRUE), , drop = FALSE]
    fb <- tryCatch(suppressWarnings(suppressMessages(
      roy(fit$spec$outcome, fit$spec$selection, db, method = method, dependence = dependence,
          type = if (is.null(type)) "generalized" else type,
          correction = if (is.null(correction)) "dmf" else correction, control = control))),
      error = function(e) NULL)
    if (!is.null(fb) && identical(names(fb$par), names(fit$par))) draws[b, ] <- fb$par
  }
  ok <- stats::complete.cases(draws)
  if (sum(ok) < B / 2) warning("more than half of the bootstrap replications failed")
  fit$vcov <- stats::cov(draws[ok, , drop = FALSE])
  fit$boot <- draws[ok, , drop = FALSE]
  fit$se_type <- paste0("bootstrap (B = ", sum(ok), ")")
  fit
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Methods for fitted Roy models
#'
#' @param object,x A `roy_fit` object.
#' @param part For `coef()`: `"all"` (natural-scale parameter vector) or one
#'   of `"selection"`, `"outcome"`, `"sigma"`, `"dependence"`, `"correction"`.
#' @param ... Unused.
#' @name roy_fit_methods
NULL

#' @rdname roy_fit_methods
#' @export
coef.roy_fit <- function(object, part = "all", ...) {
  if (part == "all") object$par else object$coef[[part]]
}

#' @rdname roy_fit_methods
#' @export
vcov.roy_fit <- function(object, ...) object$vcov

#' @rdname roy_fit_methods
#' @export
logLik.roy_fit <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$spec$n, class = "logLik")
}

#' @rdname roy_fit_methods
#' @export
nobs.roy_fit <- function(object, ...) object$spec$n

#' @rdname roy_fit_methods
#' @export
summary.roy_fit <- function(object, ...) {
  est <- object$par
  se <- sqrt(pmax(diag(object$vcov)[names(est)], 0))
  z <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(fit = object, table = tab), class = "summary.roy_fit")
}

#' @rdname roy_fit_methods
#' @export
print.summary.roy_fit <- function(x, ...) {
  f <- x$fit
  .print_header(f)
  tab <- x$table
  out <- data.frame(Estimate = formatC(tab[, 1], digits = 4, format = "f"),
                    `Std. Error` = ifelse(is.na(tab[, 2]), "", formatC(tab[, 2], digits = 4, format = "f")),
                    `z value` = ifelse(is.na(tab[, 3]), "", formatC(tab[, 3], digits = 2, format = "f")),
                    `Pr(>|z|)` = ifelse(is.na(tab[, 4]), "", format.pval(tab[, 4], digits = 3, eps = 1e-4)),
                    ` ` = ifelse(is.na(tab[, 4]), "", as.character(cut(tab[, 4], c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                                                                      c("***", "**", "*", ".", "")))),
                    check.names = FALSE, row.names = rownames(tab))
  print(out, right = TRUE)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")
  if (!is.na(f$loglik)) cat("\nlog-likelihood:", format(f$loglik, digits = 7),
                            " (df =", f$npar, ")\n")
  cat("Standard errors:", f$se_type, "\n")
  invisible(x)
}

#' @rdname roy_fit_methods
#' @export
print.roy_fit <- function(x, ...) {
  .print_header(x)
  print(round(x$par, 4))
  invisible(x)
}

.print_header <- function(f) {
  lab <- if (f$model == "A") {
    sprintf("Model A (binary AI adoption), %s Roy, %s, dependence: %s",
            f$type, f$method, if (is.null(f$family)) "-" else f$family$name)
  } else sprintf("Model B (multinomial: %s), %s, correction: %s",
                 paste(f$spec$regimes, collapse = " / "), f$method, f$type)
  cat("RoyAI fit:", lab, "\n")
  cat("n =", f$spec$n, "| regime sizes:", paste(tabulate(f$spec$d + 1L, f$spec$J), collapse = ", "),
      "| converged:", f$converged, "\n\n")
}
