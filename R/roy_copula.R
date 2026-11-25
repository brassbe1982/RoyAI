# =============================================================================
# roy_copula.R  --  Model C: generalized Roy model with copula dependence
# -----------------------------------------------------------------------------
# Part of the RoyAI package ("Generalized Roy Model of Human-AI
# Complementarity").  Implements the copula-based endogenous switching
# likelihood that replaces joint normality of (U_0, ..., U_{J-1}, V).
#
# Model (J ordered regimes, J = 2 is the binary Roy model):
#   Y_j  = X' beta_j + U_j,          U_j ~ N(0, sigma_j^2)        j = 0..J-1
#   S*   = Z' gamma + V,             V   ~ N(0, 1)
#   S    = j  iff  kappa_j < S* <= kappa_{j+1}
#          kappa_0 = -Inf, kappa_1 = 0, kappa_J = +Inf
#   (F_j(U_j), F_V(V)) ~ C_j(. , . ; tau_j)   bivariate copula, j = 0..J-1
#   U_j and U_k (j != k) are conditionally independent given V; their joint
#   dependence is never observed and is not identified (Heckman-Vytlacil).
#
# For J = 2 this is D = 1{Z'gamma + V > 0}: higher V = more inclined to
# adopt AI.  A positive tau_1 means people with high Y_1 shocks select into
# AI (positive selection on gains); a negative tau_0 combined with a positive
# tau_1 is the classical Roy comparative-advantage pattern.
#
# Copula families (parameterised by Kendall's tau in (-1, 1) so that all
# families are on a common, interpretable scale):
#   gaussian  rho = sin(pi tau / 2)
#   t         rho = sin(pi tau / 2), degrees of freedom nu > 2 (estimated)
#   clayton   theta = 2|tau| / (1 - |tau|); 90-degree rotation if tau < 0
#   gumbel    theta = 1 / (1 - |tau|);      90-degree rotation if tau < 0
#   frank     theta solves tau = 1 - 4/theta (1 - D_1(theta)); any sign
#
# Convention:  h(v | u) = dC(u, v) / du   (u = outcome uniform, v = selection
# uniform).  Everything is base R + numDeriv.
# =============================================================================

.roy_eps <- 1e-15

.clamp01 <- function(x, eps = .roy_eps) pmin(pmax(x, eps), 1 - eps)

roy_copula_families <- c("gaussian", "t", "clayton", "gumbel", "frank")

# ---- Frank: tau <-> theta ---------------------------------------------------
.debye1 <- function(theta) {
  vapply(theta, function(th) {
    if (abs(th) < 1e-8) return(1)
    stats::integrate(function(t) ifelse(t == 0, 1, t / expm1(t)),
                     0, th, rel.tol = 1e-10)$value / th
  }, numeric(1))
}

frank_tau <- function(theta) {
  ifelse(abs(theta) < 1e-8, 0, 1 - 4 / theta * (1 - .debye1(theta)))
}

# Monotone interpolation table, built once on first use.
.roy_cache <- new.env(parent = emptyenv())
.frank_theta_from_tau <- function(tau) {
  if (is.null(.roy_cache$frank)) {
    th <- c(-rev(exp(seq(log(1e-4), log(300), length.out = 600))), 0,
            exp(seq(log(1e-4), log(300), length.out = 600)))
    ta <- frank_tau(th)
    .roy_cache$frank <- stats::splinefun(ta, th, method = "monoH.FC")
  }
  .roy_cache$frank(tau)
}

#' Map Kendall's tau to native copula parameters
#'
#' @param family one of roy_copula_families
#' @param tau Kendall's tau in (-1, 1)
#' @param nu degrees of freedom (t copula only)
#' @return list(family, theta, rot, nu, tau)
#' @export
copula_par <- function(family, tau, nu = NULL) {
  family <- match.arg(family, roy_copula_families)
  tau <- max(min(tau, 0.995), -0.995)
  out <- list(family = family, tau = tau, rot = 0, nu = nu)
  out$theta <- switch(family,
    gaussian = sin(pi * tau / 2),
    t        = sin(pi * tau / 2),
    clayton  = 2 * abs(tau) / (1 - abs(tau)),
    gumbel   = 1 / (1 - abs(tau)),
    frank    = .frank_theta_from_tau(tau))
  if (family %in% c("clayton", "gumbel") && tau < 0) out$rot <- 90
  if (family == "t" && is.null(nu)) stop("t copula needs nu")
  out
}

# ---- base (unrotated, exchangeable) h-functions: dC_b(u, v)/du ---------------
.h_base <- function(v, u, cp) {
  u <- .clamp01(u); v <- .clamp01(v); th <- cp$theta
  h <- switch(cp$family,
    gaussian = {
      if (abs(th) < 1e-12) v else
        stats::pnorm((stats::qnorm(v) - th * stats::qnorm(u)) / sqrt(1 - th^2))
    },
    t = {
      nu <- cp$nu; x <- stats::qt(u, nu); y <- stats::qt(v, nu)
      stats::pt((y - th * x) / sqrt((nu + x^2) * (1 - th^2) / (nu + 1)), nu + 1)
    },
    clayton = {
      if (th < 1e-10) v else
        exp(-(1 + 1 / th) * log1p(u^th * expm1(-th * log(v))))
    },
    gumbel = {
      if (th - 1 < 1e-10) v else {
        x <- -log(u); y <- -log(v)
        A <- (x^th + y^th)^(1 / th)
        exp(-A + (1 - th) * log(A) + (th - 1) * log(x) + x)   # C * A^(1-th) * x^(th-1) / u
      }
    },
    frank = {
      if (abs(th) < 1e-8) v else {
        num <- exp(-th * u) * expm1(-th * v)
        den <- expm1(-th) + expm1(-th * u) * expm1(-th * v)
        num / den
      }
    })
  .clamp01(h, 1e-300)
}

#' Conditional copula CDF h(v | u) = dC(u, v)/du (handles rotation)
#' @param v selection uniform(s)
#' @param u outcome uniform(s)
#' @param cp copula parameter list from [copula_par()]
#' @export
copula_hfunc <- function(v, u, cp) {
  if (cp$rot == 90) .h_base(v, 1 - u, cp) else .h_base(v, u, cp)
}

#' Conditional copula CDF of u given v: dC(u, v)/dv (handles rotation)
#' @param u outcome uniform(s)
#' @param v selection uniform(s)
#' @param cp copula parameter list from [copula_par()]
#' @export
copula_hfunc2 <- function(u, v, cp) {
  if (cp$rot == 90) 1 - .h_base(1 - u, v, cp) else .h_base(u, v, cp)
}

# Vectorised bisection inverse of u -> copula_hfunc2(u, v): draw u | v.
.hinv2 <- function(t, v, cp, iter = 55) {
  lo <- rep(0, length(t)); hi <- rep(1, length(t))
  for (k in seq_len(iter)) {
    mid <- (lo + hi) / 2
    up  <- copula_hfunc2(mid, v, cp) < t
    lo[up] <- mid[up]; hi[!up] <- mid[!up]
  }
  .clamp01((lo + hi) / 2)
}

# ---- parameter vector bookkeeping ---------------------------------------------
# `sp` is the model spec: list(kx, kz, J, family, fix_nu, type, link)
#   type = "ordered": single latent index Z'gamma + V, cutpoints kappa
#   type = "mnl"    : Lee (1983) single-index copula switching with
#                     multinomial-logit selection (J >= 3 only; for J = 2 the
#                     binary logit model is the ordered model with link "logit")
.roy_layout <- function(sp) {
  J <- sp$J
  c(unlist(lapply(seq_len(J) - 1, function(j) paste0("beta", j, ":", sp$kx$names))),
    if (sp$type == "ordered") paste0("gamma:", sp$kz$names) else
      unlist(lapply(seq_len(J - 1), function(j) paste0("gamma", j, ":", sp$kz$names))),
    if (sp$type == "ordered" && J > 2) paste0("log_dkappa", 2:(J - 1)),
    paste0("log_sigma", seq_len(J) - 1),
    paste0("atanh_tau", seq_len(J) - 1),
    if (sp$family == "t" && is.null(sp$fix_nu)) "log_nu_m2")
}

.roy_unpack <- function(par, sp) {
  if (!is.null(sp$unpack)) return(sp$unpack(par))   # bridge from core roy_fit
  p <- sp$kx$k; q <- sp$kz$k; J <- sp$J; i <- 0
  take <- function(n) { r <- par[i + seq_len(n)]; i <<- i + n; r }
  beta  <- lapply(seq_len(J), function(j) take(p))
  if (sp$type == "ordered") {
    gamma <- take(q)
    kappa <- c(-Inf, 0)
    if (J > 2) kappa <- c(kappa, cumsum(exp(take(J - 2))))
    kappa <- c(kappa, Inf)
  } else {
    gamma <- matrix(take(q * (J - 1)), q, J - 1)
    kappa <- NULL
  }
  sigma <- exp(take(J))
  tau   <- tanh(take(J))
  nu    <- if (sp$family == "t") { if (is.null(sp$fix_nu)) 2 + exp(take(1)) else sp$fix_nu } else NULL
  list(beta = beta, gamma = gamma, kappa = kappa, sigma = sigma, tau = tau, nu = nu)
}

.linkfun <- function(link) if (link == "probit") stats::pnorm else stats::plogis

.mnl_prob <- function(Z, gamma) {
  e <- cbind(0, Z %*% gamma); e <- e - apply(e, 1, max)
  P <- exp(e); P / rowSums(P)
}

# ---- log-likelihood -----------------------------------------------------------
.roy_loglik_i <- function(par, y, X, Z, s, sp) {
  th <- .roy_unpack(par, sp)
  F <- .linkfun(sp$link)
  if (sp$type == "ordered") zg <- drop(Z %*% th$gamma) else P <- .mnl_prob(Z, th$gamma)
  ll <- numeric(length(y))
  for (j in seq_len(sp$J)) {
    idx <- which(s == j - 1)
    if (!length(idx)) next
    sg <- th$sigma[j]
    e  <- (y[idx] - drop(X[idx, , drop = FALSE] %*% th$beta[[j]])) / sg
    u  <- stats::pnorm(e)
    cp <- copula_par(sp$family, th$tau[j], th$nu)
    if (sp$type == "ordered") {
      # P(kappa_j < Z'g + V <= kappa_{j+1} | U_j) on the W = F(V) scale
      hhi <- if (j < sp$J) copula_hfunc(F(th$kappa[j + 1] - zg[idx]), u, cp) else 1
      hlo <- if (j > 1)    copula_hfunc(F(th$kappa[j]     - zg[idx]), u, cp) else 0
      pr  <- hhi - hlo
    } else {
      # Lee: S = j iff T_j > 1 - P_j, (U_j, T_j) ~ C_j  (tau_j > 0: high-U_j
      # people are more likely to choose regime j)
      pr <- 1 - copula_hfunc(1 - P[idx, j], u, cp)
    }
    ll[idx] <- stats::dnorm(e, log = TRUE) - log(sg) + log(pmax(pr, 1e-300))
  }
  ll
}

#' Fit the copula generalized Roy model (Model C) by maximum likelihood
#'
#' @param outcome formula for potential outcomes, e.g. y ~ x1 + x2
#'   (same regressors in every regime)
#' @param selection formula for selection, e.g. d ~ x1 + z1 + z2; the
#'   left-hand side is the observed regime, coded 0, 1, ..., J-1
#' @param data data.frame
#' @param family copula family: "gaussian", "t", "clayton", "gumbel", "frank"
#' @param type "ordered" (single latent AI-propensity index, the PDF's Model C
#'   with one U_D) or "mnl" (Lee 1983 copula switching with multinomial-logit
#'   selection; only (U_s, T_s) pairs are identified, so ATT-type effects
#'   across regimes are not available when J > 2)
#' @param link "probit" (V ~ N(0,1)) or "logit" for the ordered model
#' @param fix_nu fix the t-copula degrees of freedom (NULL = estimate)
#' @param start optional starting values (unconstrained scale)
#' @param multistart maximise from several starts and keep the best regular
#'   maximum. `NULL` (default) means `TRUE` for non-Gaussian families when
#'   `start` is not given. Starts: the two-step values, independence (all
#'   tau = 0) and the two-step taus with flipped signs. Copula likelihoods are
#'   often bimodal in the sign of tau and can spike as |tau| -> 1; a run is
#'   kept only if it is a regular maximum (positive definite Hessian, max
#'   |score| <= 1, every |tau| <= 0.98). Runs are stored in `fit$multistart`.
#' @param hessian compute the observed-information covariance matrix
#' @param control passed to optim()
#' @return object of class "roy_copula"
#' @export
roy_copula_fit <- function(outcome, selection, data,
                           family = c("gaussian", "t", "clayton", "gumbel", "frank"),
                           type = c("ordered", "mnl"), link = c("probit", "logit"),
                           fix_nu = NULL, start = NULL, multistart = NULL, hessian = TRUE,
                           control = list(maxit = 2000, reltol = 1e-10)) {
  family <- match.arg(family); type <- match.arg(type); link <- match.arg(link)
  cl <- match.call()
  y <- stats::model.response(stats::model.frame(outcome, data))
  X <- stats::model.matrix(outcome, data)
  mfz <- stats::model.frame(selection, data)
  s <- as.integer(as.character(stats::model.response(mfz)))
  Z <- stats::model.matrix(selection, data)
  J <- max(s) + 1L
  if (J < 2) stop("selection variable must take at least two values 0..J-1")
  if (!all(sort(unique(s)) == 0:(J - 1))) stop("regimes must be coded 0..J-1")
  if (type == "mnl") { link <- "logit"; if (J == 2) type <- "ordered" }
  sp <- list(kx = list(k = ncol(X), names = colnames(X)),
             kz = list(k = ncol(Z), names = colnames(Z)),
             J = J, family = family, fix_nu = fix_nu, type = type, link = link)
  nms <- .roy_layout(sp)

  user_start <- !is.null(start)
  if (is.null(multistart)) multistart <- !user_start && family != "gaussian"
  if (!user_start) start <- .roy_start(y, X, Z, s, sp)
  names(start) <- nms
  itau <- grep("^atanh_tau", nms)
  starts <- list(twostep = start)
  if (isTRUE(multistart)) {
    st <- start; st[itau] <- 0; starts$independence <- st
    st <- start; st[itau] <- -start[itau]; starts$flipped_tau <- st
  }

  nll <- function(p) {
    v <- -sum(.roy_loglik_i(p, y, X, Z, s, sp))
    if (!is.finite(v)) 1e100 else v
  }
  run1 <- function(st) {
    opt <- stats::optim(st, nll, method = "BFGS", control = control)
    if (opt$convergence != 0) {        # polish only when BFGS did not converge
      opt2 <- stats::nlminb(opt$par, nll, control = list(eval.max = 2000, iter.max = 1000))
      if (opt2$objective <= opt$value) {
        opt$par <- opt2$par; opt$value <- opt2$objective; opt$convergence <- opt2$convergence
      }
    }
    names(opt$par) <- nms
    opt
  }
  runs <- lapply(starts, run1)
  vals <- vapply(runs, function(o) o$value, numeric(1))

  # Regular interior maximum: PD Hessian, small score, |tau| away from 1.
  # Copula likelihoods can spike as |tau| -> 1, giving a higher but spurious
  # log-likelihood with a non-PD Hessian; such runs are rejected.
  regular <- function(o, H) {
    ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA_real_)
    g <- tryCatch(numDeriv::grad(nll, o$par), error = function(e) NA_real_)
    all(is.finite(ev)) && min(ev) > 0 && all(is.finite(g)) && max(abs(g)) <= 1 &&
      all(abs(tanh(o$par[itau])) <= 0.98)
  }
  ord <- order(vals); reg <- rep(NA, length(runs)); Hs <- vector("list", length(runs))
  pick <- NA_integer_
  check <- hessian || length(runs) > 1
  for (i in ord) {
    if (!check) { pick <- i; break }
    Hs[[i]] <- stats::optimHess(runs[[i]]$par, nll)
    reg[i] <- regular(runs[[i]], Hs[[i]])
    if (reg[i]) { pick <- i; break }
  }
  if (is.na(pick)) {
    pick <- ord[1]
    warning("no start reached a regular maximum (positive definite Hessian, |tau| <= 0.98); ",
            "the best run is returned, treat its standard errors with caution", call. = FALSE)
  } else if (pick != ord[1]) {
    warning(sprintf(paste("the highest log-likelihood (start '%s') is an irregular boundary point;",
                          "the best regular maximum (start '%s') is used"),
                    names(starts)[ord[1]], names(starts)[pick]), call. = FALSE)
  }
  opt <- runs[[pick]]

  V <- NULL
  if (hessian) {
    H <- Hs[[pick]]
    V <- tryCatch(solve(H), error = function(e) .ginv(H))
    dimnames(V) <- list(nms, nms)
  }
  ms <- data.frame(start = names(starts), loglik = -unname(vals), regular = reg,
                   chosen = seq_along(runs) == pick)
  structure(list(call = cl, family = family, J = J, sp = sp,
                 selection_type = type, link = link,
                 par = opt$par, vcov_raw = V,
                 theta = .roy_unpack(opt$par, sp), loglik = -opt$value, n = length(y),
                 npar = length(opt$par), convergence = opt$convergence,
                 outcome = outcome, selection = selection, multistart = ms,
                 X = X, Z = Z, s = s, y = y),
            class = "roy_copula")
}

.ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s <- svd(X); pos <- s$d > max(tol * s$d[1L], 0)
  s$v[, pos, drop = FALSE] %*% ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}

# ---- two-step (normal) estimator: starting values and MC benchmark -------------
#' Two-step Heckman / Lee estimator under joint normality
#'
#' Step 1: ordered probit (or MNL) for selection.  Step 2: OLS of Y on X and
#' the generalized inverse Mills ratio within each regime.  Returns beta,
#' sigma, the implied Gaussian-copula tau, and the selection parameters.
#' @param y outcome vector
#' @param X outcome design matrix
#' @param Z selection design matrix
#' @param s observed regime (0..J-1)
#' @param sp model spec list(kx, kz, J, family, fix_nu, type, link)
#' @return parameter vector on the unconstrained scale (Gaussian layout)
#' @export
roy_twostep <- function(y, X, Z, s, sp) {
  J <- sp$J; F <- .linkfun(sp$link)
  if (sp$type == "ordered") {
    nll_op <- function(p) {
      g <- p[seq_len(sp$kz$k)]
      kap <- c(-Inf, 0, if (J > 2) cumsum(exp(p[sp$kz$k + seq_len(J - 2)])), Inf)
      zg <- drop(Z %*% g)
      pr <- F(kap[s + 2] - zg) - F(kap[s + 1] - zg)
      -sum(log(pmax(pr, 1e-300)))
    }
    op <- stats::optim(c(rep(0, sp$kz$k), rep(0, max(J - 2, 0))), nll_op, method = "BFGS")$par
    g <- op[seq_len(sp$kz$k)]
    kap <- c(-Inf, 0, if (J > 2) cumsum(exp(op[sp$kz$k + seq_len(J - 2)])), Inf)
    zg <- drop(Z %*% g)
    selpar <- op
  } else {
    nll_mnl <- function(p) {
      P <- .mnl_prob(Z, matrix(p, sp$kz$k, J - 1))
      -sum(log(pmax(P[cbind(seq_along(s), s + 1)], 1e-300)))
    }
    selpar <- stats::optim(rep(0, sp$kz$k * (J - 1)), nll_mnl, method = "BFGS")$par
    P <- .mnl_prob(Z, matrix(selpar, sp$kz$k, J - 1))
  }
  b <- list(); lsig <- numeric(J); atau <- numeric(J)
  for (j in seq_len(J)) {
    idx <- s == j - 1
    if (sp$type == "ordered") {
      # E[V* | a < V <= b] with V* = Phi^{-1}(F(V)) (normal scores of W)
      a <- stats::qnorm(.clamp01(F(kap[j] - zg[idx])))
      bb <- stats::qnorm(.clamp01(F(kap[j + 1] - zg[idx])))
      lam <- (stats::dnorm(a) - stats::dnorm(bb)) /
             pmax(stats::pnorm(bb) - stats::pnorm(a), 1e-12)
    } else {
      # Lee: E[Phi^{-1}(T_j) | T_j > 1 - P_j] = phi(Phi^{-1}(P_j)) / P_j
      pj <- .clamp01(P[idx, j])
      lam <- stats::dnorm(stats::qnorm(pj)) / pj
    }
    fit <- stats::lm.fit(cbind(X[idx, , drop = FALSE], lam), y[idx])
    cf <- fit$coefficients; cf[is.na(cf)] <- 0
    b[[j]] <- cf[seq_len(sp$kx$k)]
    sg <- sqrt(max(mean(fit$residuals^2) + cf[sp$kx$k + 1]^2 * 0.5, 0.05))
    rho <- max(min(cf[sp$kx$k + 1] / sg, 0.9), -0.9)
    lsig[j] <- log(sg)
    atau[j] <- atanh(2 / pi * asin(rho))
  }
  c(unlist(b), selpar, lsig, atau)
}

.roy_start <- function(y, X, Z, s, sp) {
  c(roy_twostep(y, X, Z, s, sp),
    if (sp$family == "t" && is.null(sp$fix_nu)) log(8 - 2))
}

# ---- natural-scale coefficients and delta-method covariance -------------------
.roy_natural <- function(par, sp) {
  th <- .roy_unpack(par, sp); J <- sp$J
  ord <- sp$type == "ordered"
  out <- c(unlist(th$beta), as.vector(th$gamma),
           if (ord && J > 2) th$kappa[3:J], th$sigma, th$tau, th$nu)
  nm <- c(unlist(lapply(seq_len(J) - 1, function(j) paste0("beta", j, ":", sp$kx$names))),
          if (ord) paste0("gamma:", sp$kz$names) else
            unlist(lapply(seq_len(J - 1), function(j) paste0("gamma", j, ":", sp$kz$names))),
          if (ord && J > 2) paste0("kappa", 2:(J - 1)),
          paste0("sigma", seq_len(J) - 1), paste0("tau", seq_len(J) - 1),
          if (sp$family == "t") "nu")
  stats::setNames(out, nm)
}

#' @export
coef.roy_copula <- function(object, ...) .roy_natural(object$par, object$sp)

#' @export
vcov.roy_copula <- function(object, ...) {
  if (is.null(object$vcov_raw)) return(NULL)
  G <- numDeriv::jacobian(function(p) .roy_natural(p, object$sp), object$par)
  V <- G %*% object$vcov_raw %*% t(G)
  nm <- names(coef(object)); dimnames(V) <- list(nm, nm); V
}

#' @export
logLik.roy_copula <- function(object, ...) {
  structure(object$loglik, df = object$npar, nobs = object$n, class = "logLik")
}

#' @export
nobs.roy_copula <- function(object, ...) object$n

#' @export
print.roy_copula <- function(x, ...) {
  cat("Generalized Roy model with ", x$family, " copula dependence (Model C, ",
      x$selection_type, " selection, ", x$link, ")\n", sep = "")
  cat("Regimes:", x$J, " n =", x$n, " logLik =", format(x$loglik, digits = 8),
      " AIC =", format(stats::AIC(x), digits = 8), "\n")
  print(round(coef(x), 4))
  invisible(x)
}

#' @export
summary.roy_copula <- function(object, ...) {
  b <- coef(object); V <- vcov(object)
  se <- if (is.null(V)) NA else sqrt(pmax(diag(V), 0))
  tab <- cbind(Estimate = b, `Std. Error` = se, `z value` = b / se,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(b / se)))
  structure(list(fit = object, coefficients = tab), class = "summary.roy_copula")
}

#' @export
print.summary.roy_copula <- function(x, ...) {
  print(x$fit)
  cat("\nCoefficients (sigma, tau, nu on natural scale; delta-method SEs):\n")
  stats::printCoefmat(x$coefficients, digits = 4)
  invisible(x)
}

#' Fit every copula family and rank them by AIC/BIC
#' @param outcome,selection,data as in [roy_copula_fit()]
#' @param families copula families to compare
#' @param ... passed to [roy_copula_fit()]
#' @return list(table, fits, best)
#' @export
roy_copula_select <- function(outcome, selection, data,
                              families = roy_copula_families, ...) {
  fits <- lapply(families, function(f)
    tryCatch(roy_copula_fit(outcome, selection, data, family = f, ...),
             error = function(e) NULL))
  names(fits) <- families
  ok <- !vapply(fits, is.null, logical(1))
  tab <- data.frame(family = families[ok],
                    logLik = vapply(fits[ok], function(f) f$loglik, 0),
                    npar = vapply(fits[ok], function(f) f$npar, 0),
                    AIC = vapply(fits[ok], stats::AIC, 0),
                    BIC = vapply(fits[ok], stats::BIC, 0))
  tab <- tab[order(tab$AIC), ]
  list(table = tab, fits = fits, best = fits[[tab$family[1]]])
}

# ---- data generating process ------------------------------------------------
#' Simulate from the copula generalized Roy model (Model C)
#'
#' @param n sample size
#' @param beta list of J coefficient vectors (intercept first) on X = (1, x1, ...)
#' @param gamma selection coefficients on Z = (1, x1, z1, ...); for
#'   type = "mnl" a matrix with one column per regime 1..J-1
#' @param sigma outcome error sds (length J)
#' @param tau Kendall's tau between each U_j and its selection uniform (length J)
#' @param family copula family
#' @param nu t-copula df
#' @param kappa interior cutpoints beyond 0 (length J-2, increasing; ordered)
#' @param type "ordered" or "mnl" (Lee 1983 DGP with i.i.d. Gumbel utilities)
#' @param link "probit" or "logit" (ordered model)
#' @param covariates data.frame of covariates; default x1 ~ N(0,1),
#'   z1 ~ N(0,1), z2 ~ Bernoulli(0.5).  Columns are taken in order:
#'   X uses the first length(beta[[1]]) - 1 columns, Z the first nrow(gamma) - 1.
#' @return data.frame with y, s (observed regime), potential outcomes
#'   y0..y{J-1}, selection uniform(s) and covariates
#' @export
roy_copula_sim <- function(n, beta, gamma, sigma, tau, family = "gaussian",
                           nu = NULL, kappa = NULL, type = c("ordered", "mnl"),
                           link = c("probit", "logit"), covariates = NULL) {
  type <- match.arg(type); link <- match.arg(link)
  J <- length(beta)
  if (is.null(covariates)) {
    covariates <- data.frame(x1 = stats::rnorm(n), z1 = stats::rnorm(n),
                             z2 = stats::rbinom(n, 1, 0.5))
  }
  gamma <- as.matrix(gamma)
  X <- cbind(1, as.matrix(covariates[, seq_len(length(beta[[1]]) - 1), drop = FALSE]))
  Z <- cbind(1, as.matrix(covariates[, seq_len(nrow(gamma) - 1), drop = FALSE]))
  if (type == "ordered") {
    w <- stats::runif(n)                          # selection uniform F_V(V)
    v <- if (link == "probit") stats::qnorm(w) else stats::qlogis(w)
    kap <- c(-Inf, 0, kappa, Inf)
    s <- findInterval(drop(Z %*% gamma[, 1]) + v, kap, left.open = TRUE) - 1L
    W <- matrix(w, n, J)
  } else {
    mu <- cbind(0, Z %*% gamma)
    eps <- -log(-log(matrix(stats::runif(n * J), n, J)))      # Gumbel
    ut <- mu + eps
    s <- max.col(ut, ties.method = "first") - 1L
    # Lee transform: eta_s = max_{k != s}(mu_k + eps_k) - eps_s ~ logistic(L_s),
    # T_s = F(eta_s), S = s iff T_s <= P_s.  We use W_s = 1 - T_s.
    W <- sapply(seq_len(J), function(j) {
      L <- log(rowSums(exp(mu[, -j, drop = FALSE])))
      eta <- apply(ut[, -j, drop = FALSE], 1, max) - eps[, j]
      1 - stats::plogis(eta - L)
    })
  }
  ypot <- matrix(NA_real_, n, J)
  for (j in seq_len(J)) {
    cp <- copula_par(family, tau[j], nu)
    uj <- .hinv2(stats::runif(n), W[, j], cp)
    ypot[, j] <- drop(X %*% beta[[j]]) + sigma[j] * stats::qnorm(uj)
  }
  colnames(ypot) <- paste0("y", seq_len(J) - 1)
  y <- ypot[cbind(seq_len(n), s + 1L)]
  data.frame(y = y, s = s, ypot, covariates)
}

# ---- bridge to the RoyAI core (roy(), roy_sim_A(), roy_effects()) -------------
#' Copula families in the core plug-in format
#'
#' Builds [roy_dependence()] objects for the Clayton, Gumbel, Frank and
#' Student-t copulas so that `roy(..., dependence = "clayton")` estimates
#' Model A through the core.  The dependence parameter is Kendall's tau
#' (t copula: tau and nu), so families are directly comparable.  Negative
#' tau for Clayton/Gumbel uses the 90-degree rotation; the core simulator
#' assumes exchangeable copulas, so simulate rotated cases with
#' [roy_copula_sim()].
#' @return named list of roy_dependence objects
#' @export
roy_copula_dependence <- function() {
  mk <- function(fam) {
    if (fam == "t") {
      return(roy_dependence(
        name = "t", npar = 2L, par_names = c("tau", "nu"),
        h = function(v, u, theta) copula_hfunc(v, u, copula_par("t", theta[1], theta[2])),
        linkinv = function(eta) c(tanh(eta[1]), 2 + exp(eta[2])),
        link = function(theta) c(atanh(theta[1]), log(theta[2] - 2)),
        start = c(0, 8), tau = function(theta) theta[1]))
    }
    roy_dependence(
      name = fam, npar = 1L, par_names = "tau",
      h = function(v, u, theta) copula_hfunc(v, u, copula_par(fam, theta[1])),
      linkinv = function(eta) tanh(eta), link = function(theta) atanh(theta),
      start = 0, tau = function(theta) theta[1])
  }
  fams <- c("clayton", "gumbel", "frank", "t")
  stats::setNames(lapply(fams, mk), fams)
}

#' Register the copula families with the RoyAI core (called from .onLoad)
#' @export
roy_copula_register <- function() {
  if (!exists("roy_register_dependence", mode = "function")) return(invisible(FALSE))
  for (f in roy_copula_dependence()) roy_register_dependence(f, overwrite = TRUE)
  invisible(TRUE)
}

# Convert a core Model A fit (class "roy_fit_A") into the internal form
# used by roy_te(), roy_mte() and roy_late(); parameters stay on the core's
# natural scale and fit$vcov is used for the delta method.
.as_roy_copula <- function(fit) {
  if (inherits(fit, "roy_copula")) return(fit)
  if (!inherits(fit, "roy_fit_A"))
    stop("roy_te/roy_mte/roy_late need a roy_copula fit or a core Model A fit (roy_fit_A)")
  famname <- if (is.list(fit$family)) fit$family$name else as.character(fit$family)
  ours <- famname %in% c("clayton", "gumbel", "frank", "t")
  if (!ours && !famname %in% c("gaussian", "independence"))
    stop("dependence family '", famname, "' has no copula mapping")
  fam <- if (ours) famname else "gaussian"
  unpack <- function(p) {
    cl <- fit$relist(p)
    dep <- cl$dependence
    th0 <- if (length(dep)) unlist(dep[[1]]) else 0
    th1 <- if (length(dep)) unlist(dep[[2]]) else 0
    tau <- if (ours) c(th0[1], th1[1]) else 2 / pi * asin(c(th0[1], th1[1]))
    list(beta = list(unname(cl$outcome[[1]]), unname(cl$outcome[[2]])),
         gamma = unname(cl$selection), kappa = c(-Inf, 0, Inf),
         sigma = unname(cl$sigma), tau = unname(tau),
         nu = if (fam == "t") unname(th1[2]) else NULL)
  }
  sp <- list(J = 2, family = fam, type = "ordered", link = "probit", unpack = unpack)
  structure(list(par = fit$par, vcov_raw = fit$vcov, sp = sp, J = 2, family = fam,
                 link = "probit", selection_type = "ordered",
                 X = fit$spec$X, Z = fit$spec$Z,
                 outcome = fit$spec$outcome, selection = fit$spec$selection),
            class = "roy_copula")
}
