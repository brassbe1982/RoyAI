# =============================================================================
# roy_dynamic.R  --  Model E: dynamic Roy selection with learning about
#                    human-AI complementarity and learning by doing
# -----------------------------------------------------------------------------
# Monograph Section 5.5: Bellman equation (eq:bellman), Kalman learning
# (eq:kalman), learning by doing (eq:lbd), Hotz-Miller inversion (lem:hm).
#
# Two arrangements, s in {0 = human-only, 1 = AI}. Periods t = 1..T, discount
# factor beta (fixed: flow payoffs are not identified jointly with beta,
# Magnac and Thesmar 2002).
#   human-only earnings   y_it(0) = mu0 + eta_it,          eta ~ N(0, sigma0^2)
#   AI earnings           y_it(1) = mu0 + kappa_it + e_it, e   ~ N(0, sigma_e^2)
#   complementarity       kappa_i1 ~ N(mu_k, sd_k^2),
#                         kappa_i,t+1 = kappa_it + psi 1(s_it = 1)
#   beliefs               kappa_it | history ~ N(m_it, vs_n), n = # past AI periods
#                         (Kalman recursion; the worker's prior is the
#                         population distribution, so beliefs are computable
#                         from observed earnings histories)
#   flow utility          u_1 - u_0 = m_it - c_i,  c_i = X_i'theta_c (cost)
#   shocks                type-I extreme value with scale lambda
#
# Because beliefs depend only on observed histories, the likelihood factorises
# into an earnings part (a Kalman prediction-error likelihood) and a choice
# part (dynamic logit), and the value function depends on the state only
# through (t, n, m - c), so one backward induction serves every worker.
# =============================================================================

# Gauss-Hermite nodes and weights for a standard normal (Golub-Welsch).
.gh_normal <- function(k = 15) {
  if (!is.null(.roy_cache[[paste0("gh", k)]])) return(.roy_cache[[paste0("gh", k)]])
  J <- matrix(0, k, k)
  for (i in seq_len(k - 1)) J[i, i + 1] <- J[i + 1, i] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  out <- list(x = e$values, w = e$vectors[1, ]^2)
  .roy_cache[[paste0("gh", k)]] <- out
  out
}

# Linear interpolation with linear extrapolation (value functions grow
# linearly in beliefs far from the decision margin).
.lin_interp <- function(xg, yg, x) {
  i <- findInterval(x, xg, all.inside = TRUE)
  x0 <- xg[i]; x1 <- xg[i + 1]
  yg[i] + (yg[i + 1] - yg[i]) * (x - x0) / (x1 - x0)
}

.rd_postvar <- function(n, sd_k, sigma_e) 1 / (1 / sd_k^2 + n / sigma_e^2)

#' Solve Model E by backward induction
#'
#' Solves the finite-horizon dynamic Roy model of the monograph (Section 5.5)
#' for two arrangements with Bayesian learning about complementarity and
#' learning by doing. The state is `(t, n, m - c)`: period, number of past AI
#' periods (which fixes the posterior variance) and the expected net current
#' gain. Expectations over the next signal use Gauss-Hermite quadrature.
#'
#' @param par list with `sd_k`, `sigma_e`, `psi`, `lambda`
#'   (other elements are ignored)
#' @param T horizon
#' @param beta discount factor
#' @param grid grid of net expected gains `m - c`; default spans the prior
#' @param gh number of Gauss-Hermite nodes
#' @return list with `grid`, `dv` (array `T x (T+1) x length(grid)` of
#'   `v_1 - v_0`), `Vbar` (expected value functions) and inputs
#' @examples
#' sol <- roy_dynamic_solve(list(sd_k = 0.25, sigma_e = 0.3, psi = 0.05,
#'                               lambda = 0.1), T = 5)
#' dim(sol$dv)
#' @export
roy_dynamic_solve <- function(par, T, beta = 0.95, grid = NULL, gh = 15) {
  sd_k <- par$sd_k; se <- par$sigma_e; psi <- par$psi; lam <- par$lambda
  if (is.null(grid))
    grid <- seq(-1 - 6 * sd_k, 1 + 6 * sd_k + abs(psi) * T, length.out = 161)
  G <- length(grid); q <- .gh_normal(gh)
  dv <- array(NA_real_, c(T, T + 1, G))
  Vb <- array(0, c(T + 1, T + 1, G))                 # Vbar[T+1, , ] = 0
  for (t in T:1) {
    for (n in 0:(t - 1)) {
      vs <- .rd_postvar(n, sd_k, se)
      sdq <- vs / sqrt(vs + se^2)                     # sd of the belief update
      Vn0 <- Vb[t + 1, n + 1, ]
      v0 <- beta * Vn0
      v1 <- grid
      if (t < T) {
        Vn1 <- Vb[t + 1, n + 2, ]
        nxt <- outer(grid + psi, sdq * q$x, "+")
        EV <- matrix(.lin_interp(grid, Vn1, as.vector(nxt)), G) %*% q$w
        v1 <- v1 + beta * drop(EV)
      }
      mx <- pmax(v0, v1)
      Vb[t, n + 1, ] <- mx + lam * log(exp((v0 - mx) / lam) + exp((v1 - mx) / lam))
      dv[t, n + 1, ] <- v1 - v0
    }
  }
  list(grid = grid, dv = dv, Vbar = Vb, T = T, beta = beta, par = par)
}

# Beliefs (prior mean m_it and n_it before each period's choice) from observed
# choices and earnings; S and Y are N x T matrices.
.rd_beliefs <- function(S, Y, par) {
  N <- nrow(S); T <- ncol(S)
  M <- matrix(NA_real_, N, T); Nn <- matrix(0L, N, T)
  m <- rep(par$mu_k, N); n <- integer(N)
  for (t in seq_len(T)) {
    M[, t] <- m; Nn[, t] <- n
    vs <- .rd_postvar(n, par$sd_k, par$sigma_e)
    k <- vs / (vs + par$sigma_e^2)
    ai <- S[, t] == 1
    q <- Y[, t] - par$mu0
    m[ai] <- m[ai] + k[ai] * (q[ai] - m[ai]) + par$psi
    n[ai] <- n[ai] + 1L
  }
  list(M = M, N = Nn)
}

.rd_earn_ll <- function(S, Y, par, B = NULL) {
  if (is.null(B)) B <- .rd_beliefs(S, Y, par)
  vs <- .rd_postvar(B$N, par$sd_k, par$sigma_e)
  ifelse(S == 1,
         stats::dnorm(Y, par$mu0 + B$M, sqrt(vs + par$sigma_e^2), log = TRUE),
         stats::dnorm(Y, par$mu0, par$sigma0, log = TRUE))
}

.rd_dv_at <- function(sol, tt, nn, mt) {
  out <- numeric(length(tt))
  for (t in unique(tt)) for (n in unique(nn[tt == t])) {
    i <- which(tt == t & nn == n)
    out[i] <- .lin_interp(sol$grid, sol$dv[t, n + 1, ], mt[i])
  }
  out
}

.rd_names_e <- c("mu0", "log_sigma0", "mu_k", "log_sd_k", "log_sigma_e", "psi")

.rd_unpack <- function(p, cnames) {
  list(mu0 = p[[1]], sigma0 = exp(p[[2]]), mu_k = p[[3]], sd_k = exp(p[[4]]),
       sigma_e = exp(p[[5]]), psi = p[[6]],
       cost = unname(p[6 + seq_along(cnames)]),
       lambda = exp(p[[7 + length(cnames)]]))
}

.rd_grid <- function(B, cvec, par, T) {
  mt <- B$M - cvec
  lo <- min(mt) - 4 * par$sd_k - 0.5; hi <- max(mt) + 4 * par$sd_k + abs(par$psi) * T + 0.5
  seq(lo, hi, length.out = 161)
}

.rd_choice_ll <- function(S, Y, par, Xc, beta, B = NULL) {
  if (is.null(B)) B <- .rd_beliefs(S, Y, par)
  T <- ncol(S); cvec <- drop(Xc %*% par$cost)
  sol <- roy_dynamic_solve(par, T, beta, grid = .rd_grid(B, cvec, par, T))
  tt <- rep(seq_len(T), each = nrow(S))
  dv <- .rd_dv_at(sol, tt, as.vector(B$N), as.vector(B$M - cvec))
  eta <- dv / par$lambda
  ifelse(as.vector(S) == 1, stats::plogis(eta, log.p = TRUE), stats::plogis(-eta, log.p = TRUE))
}

.rd_panel <- function(data, id, time, choice, outcome) {
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
  ids <- unique(data[[id]]); tt <- sort(unique(data[[time]]))
  N <- length(ids); T <- length(tt)
  if (nrow(data) != N * T) stop("Model E needs a balanced panel (every id observed in every period)")
  S <- matrix(as.integer(data[[choice]]), N, T, byrow = TRUE)
  Y <- matrix(data[[outcome]], N, T, byrow = TRUE)
  if (!all(S %in% 0:1)) stop("choice must be coded 0 (human-only) / 1 (AI)")
  first <- data[!duplicated(data[[id]]), , drop = FALSE]
  list(S = S, Y = Y, ids = ids, times = tt, first = first, data = data)
}

#' Fit Model E: dynamic Roy model with learning about AI complementarity
#'
#' Estimates the two-arrangement dynamic Roy model of the monograph (Section
#' 5.5) from a balanced panel of choices (0 = human-only, 1 = AI) and
#' observed earnings. Workers learn their complementarity `kappa_i` from AI
#' earnings by Bayesian (Kalman) updating and acquire complementarity by
#' doing at rate `psi`; they choose by comparing expected net current gains
#' `m_it - c_i` plus continuation values, with logit shocks of scale `lambda`.
#' The discount factor is fixed (it is not identified jointly with flow
#' payoffs).
#'
#' Estimation is in two stages. Stage 1 maximises the earnings likelihood
#' (Kalman prediction errors) for `mu0, sigma0, mu_k, sd_k, sigma_e, psi`;
#' this is consistent on its own because choices depend on past earnings
#' only through beliefs. Stage 2 estimates the cost coefficients and
#' `lambda`:
#' * `method = "nfxp"` (default): nested fixed point. The model is solved by
#'   backward induction ([roy_dynamic_solve()]) inside the choice likelihood;
#'   the two stages are then polished jointly on the full likelihood, so the
#'   standard errors are full-information ones.
#' * `method = "ccp"`: conditional choice probabilities (Hotz and Miller 1993,
#'   monograph Lemma "lem:hm"). A flexible first-stage logit estimates
#'   `P(s = 1 | t, n, m, X)`; because staying human-only freezes the state, the
#'   continuation-value difference is
#'   `beta sum_{tau > t} beta^(tau-t-1) [ln P0_tau(n, m) - E ln P0_tau(n+1, m')]`,
#'   and the structural parameters come from one logit with this offset. Fast;
#'   stage-2 standard errors are conditional on stage 1.
#'
#' @param data long data.frame
#' @param id,time,choice,outcome column names
#' @param cost one-sided formula for the cost of AI use in terms of
#'   time-invariant worker covariates (first row per worker is used); `~ 1`
#'   by default
#' @param beta discount factor (fixed)
#' @param method `"nfxp"` or `"ccp"`
#' @param control passed to [stats::optim()]
#' @return object of class `"roy_dynamic"`
#' @examples
#' set.seed(2)
#' p <- roy_dynamic_sim(400, T = 6)
#' fit <- roy_dynamic(p, cost = ~ x, method = "ccp")
#' fit
#' @export
roy_dynamic <- function(data, id = "id", time = "t", choice = "s", outcome = "y",
                        cost = ~ 1, beta = 0.95, method = c("nfxp", "ccp"),
                        control = list(maxit = 500, reltol = 1e-10)) {
  method <- match.arg(method); cl <- match.call()
  pn <- .rd_panel(as.data.frame(data), id, time, choice, outcome)
  S <- pn$S; Y <- pn$Y; N <- nrow(S); T <- ncol(S)
  Xc <- stats::model.matrix(cost, pn$first); cn <- colnames(Xc)
  if (!any(S == 1) || !any(S == 0)) stop("need both arrangements in the data")

  # ---- stage 1: earnings ----
  y1 <- Y[S == 1]; y0 <- Y[S == 0]
  st_e <- c(mean(y0), log(stats::sd(y0)), mean(y1) - mean(y0),
            log(0.5 * stats::sd(y1)), log(0.7 * stats::sd(y1)), 0)
  nll_e <- function(p) {
    par <- .rd_unpack(c(p, rep(0, length(cn)), 0), cn)
    v <- -sum(.rd_earn_ll(S, Y, par)); if (is.finite(v)) v else 1e100
  }
  o1 <- stats::optim(st_e, nll_e, method = "BFGS", control = control)
  H1 <- stats::optimHess(o1$par, nll_e)
  par_e <- o1$par; names(par_e) <- .rd_names_e
  pe <- .rd_unpack(c(par_e, rep(0, length(cn)), 0), cn)
  B <- .rd_beliefs(S, Y, pe)

  # ---- stage 2 ----
  nms <- c(.rd_names_e, paste0("cost:", cn), "log_lambda")
  ccp <- .rd_ccp(S, B, Xc, pe, beta)
  if (method == "ccp") {
    par <- c(par_e, ccp$par)
    V <- matrix(0, length(par), length(par)); dimnames(V) <- list(nms, nms)
    V[1:6, 1:6] <- tryCatch(solve(H1), error = function(e) .ginv(H1))
    V[-(1:6), -(1:6)] <- ccp$vcov
    loglik <- -o1$value + sum(.rd_choice_ll(S, Y, .rd_unpack(par, cn), Xc, beta, B))
    conv <- o1$convergence
  } else {
    nll_c <- function(pc) {
      v <- -sum(.rd_choice_ll(S, Y, .rd_unpack(c(par_e, pc), cn), Xc, beta, B))
      if (is.finite(v)) v else 1e100
    }
    o2 <- stats::optim(ccp$par, nll_c, method = "BFGS", control = control)
    nll <- function(p) {
      par <- .rd_unpack(p, cn); Bp <- .rd_beliefs(S, Y, par)
      v <- -sum(.rd_earn_ll(S, Y, par, Bp)) - sum(.rd_choice_ll(S, Y, par, Xc, beta, Bp))
      if (is.finite(v)) v else 1e100
    }
    o3 <- stats::optim(c(par_e, o2$par), nll, method = "BFGS", control = control)
    par <- o3$par
    H <- stats::optimHess(par, nll)
    V <- tryCatch(solve(H), error = function(e) .ginv(H))
    loglik <- -o3$value; conv <- o3$convergence
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) warning("Hessian is not positive definite at the NFXP maximum", call. = FALSE)
  }
  names(par) <- nms; dimnames(V) <- list(nms, nms)
  th <- .rd_unpack(par, cn)
  structure(list(call = cl, method = method, par = par, vcov_raw = V, theta = th,
                 loglik = loglik, N = N, T = T, beta = beta, cost = cost, cnames = cn,
                 Xc = Xc, S = S, Y = Y, ids = pn$ids, times = pn$times,
                 ccp_first_stage = ccp$first, convergence = conv,
                 beliefs = .rd_beliefs(S, Y, th)),
            class = "roy_dynamic")
}

# CCP estimator of (cost, log lambda) given stage-1 earnings parameters.
.rd_ccp <- function(S, B, Xc, pe, beta) {
  N <- nrow(S); T <- ncol(S)
  xc <- Xc[, colnames(Xc) != "(Intercept)", drop = FALSE]
  feat <- function(t, n, m, xi) {
    d <- data.frame(t = t, n = n, m = m)
    if (ncol(xc)) d <- cbind(d, xi)
    d
  }
  xi_long <- xc[rep(seq_len(N), T), , drop = FALSE]
  tt <- rep(seq_len(T), each = N); nn <- as.vector(B$N); mm <- as.vector(B$M)
  d1 <- feat(tt, nn, mm, xi_long); d1$s <- as.vector(S)
  rhs <- "t + I(t^2) + n + I(n^2) + m + I(m^2) + m:t + n:t + m:n"
  if (ncol(xc)) rhs <- paste(rhs, "+", paste(colnames(xc), collapse = " + "),
                             "+", paste0("m:", colnames(xc), collapse = " + "))
  first <- suppressWarnings(stats::glm(stats::as.formula(paste("s ~", rhs)),
                                       family = stats::binomial(), data = d1))
  lnP0 <- function(t, n, m, xi) {
    eta <- stats::predict(first, newdata = feat(t, n, m, xi), type = "link")
    stats::plogis(-eta, log.p = TRUE)
  }
  q <- .gh_normal(15)
  vs <- .rd_postvar(nn, pe$sd_k, pe$sigma_e); sdq <- vs / sqrt(vs + pe$sigma_e^2)
  Gam <- numeric(length(tt))
  for (tau in 2:T) {
    i <- which(tt < tau)
    if (!length(i)) next
    disc <- beta^(tau - tt[i])                      # beta * beta^(tau - t - 1)
    a <- lnP0(rep(tau, length(i)), nn[i], mm[i], xi_long[i, , drop = FALSE])
    b <- 0
    for (j in seq_along(q$x))
      b <- b + q$w[j] * lnP0(rep(tau, length(i)), nn[i] + 1, mm[i] + pe$psi + sdq[i] * q$x[j],
                             xi_long[i, , drop = FALSE])
    Gam[i] <- Gam[i] + disc * (a - b)
  }
  Xl <- Xc[rep(seq_len(N), T), , drop = FALSE]
  g2 <- suppressWarnings(stats::glm.fit(cbind(m = mm, Xl), d1$s, family = stats::binomial(),
                                        offset = Gam))
  b <- g2$coefficients; Vb <- stats::summary.glm(structure(g2, class = c("glm", "lm")))$cov.unscaled
  bm <- b[1]
  if (!is.finite(bm) || bm <= 0) stop("CCP second stage: coefficient on beliefs is not positive")
  # natural -> (cost = -b_x / b_m, log lambda = -log b_m)
  par <- c(-b[-1] / bm, -log(bm))
  J <- matrix(0, length(par), length(b))
  k <- length(b) - 1
  J[seq_len(k), 1] <- b[-1] / bm^2
  J[seq_len(k), -1] <- diag(-1 / bm, k)
  J[k + 1, 1] <- -1 / bm
  list(par = unname(par), vcov = J %*% Vb %*% t(J), first = first, offset = Gam)
}

.rd_natural <- function(par, cn) {
  th <- .rd_unpack(par, cn)
  out <- c(th$mu0, th$sigma0, th$mu_k, th$sd_k, th$sigma_e, th$psi, th$cost, th$lambda)
  names(out) <- c("mu0", "sigma0", "mu_k", "sd_k", "sigma_e", "psi", paste0("cost:", cn), "lambda")
  out
}

#' Methods for Model E fits
#'
#' @param object,x a `"roy_dynamic"` fit
#' @param ... unused
#' @name roy_dynamic_methods
#' @return `coef()` a named numeric vector, `vcov()` a matrix, `logLik()` an object of class `logLik`, `nobs()` an integer, `summary()` an object of class `summary.roy_dynamic`; the `print` methods return their argument invisibly.
#' @examples
#' set.seed(2)
#' p <- roy_dynamic_sim(300, T = 5)
#' fit <- roy_dynamic(p, cost = ~ x, method = "ccp")
#' coef(fit)
#' summary(fit)
NULL

#' @rdname roy_dynamic_methods
#' @export
coef.roy_dynamic <- function(object, ...) .rd_natural(object$par, object$cnames)

#' @rdname roy_dynamic_methods
#' @export
vcov.roy_dynamic <- function(object, ...) {
  J <- numDeriv::jacobian(function(p) .rd_natural(p, object$cnames), object$par)
  V <- J %*% object$vcov_raw %*% t(J); nm <- names(coef(object))
  dimnames(V) <- list(nm, nm); V
}

#' @rdname roy_dynamic_methods
#' @export
logLik.roy_dynamic <- function(object, ...)
  structure(object$loglik, df = length(object$par), nobs = object$N * object$T, class = "logLik")

#' @rdname roy_dynamic_methods
#' @export
nobs.roy_dynamic <- function(object, ...) object$N * object$T

#' @rdname roy_dynamic_methods
#' @export
print.roy_dynamic <- function(x, ...) {
  cat("Model E: dynamic Roy model with learning (", toupper(x$method), ")\n", sep = "")
  cat(sprintf("N = %d workers, T = %d periods, beta = %.2f (fixed), AI share %.1f%%, logLik %.2f\n",
              x$N, x$T, x$beta, 100 * mean(x$S), x$loglik))
  print(round(coef(x), 4))
  invisible(x)
}

#' @rdname roy_dynamic_methods
#' @export
summary.roy_dynamic <- function(object, ...) {
  b <- coef(object); se <- sqrt(pmax(diag(vcov(object)), 0))
  tab <- cbind(Estimate = b, `Std. Error` = se, `z value` = b / se,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(b / se)))
  structure(list(fit = object, coefficients = tab), class = "summary.roy_dynamic")
}

#' @rdname roy_dynamic_methods
#' @export
print.summary.roy_dynamic <- function(x, ...) {
  print(x$fit); cat("\n"); stats::printCoefmat(x$coefficients, digits = 4)
  if (x$fit$method == "ccp") cat("\nStage-2 (CCP) standard errors are conditional on stage 1.\n")
  invisible(x)
}

#' Dynamic treatment effects for Model E
#'
#' By period, from a [roy_dynamic()] fit: the AI share; the average gain
#' `E[kappa_it]` for everyone (ATE), adopters (ATT) and non-adopters (ATU),
#' estimated by posterior means (unbiased under the model's rational
#' expectations); the naive cross-sectional earnings gap; the share of
#' adopters whose expected net current gain `m_it - c_i` is negative, who
#' adopt for the option value of learning and the human capital from doing
#' (monograph Propositions "prop:dyn"(b) and "prop:experiment"); and the
#' predicted AI share of the fitted forward-looking model versus a myopic
#' worker (`beta = 0`) at the same states.
#' @param fit a `"roy_dynamic"` fit
#' @return data.frame, one row per period
#' @examples
#' set.seed(2)
#' p <- roy_dynamic_sim(300, T = 5)
#' fit <- roy_dynamic(p, cost = ~ x, method = "ccp")
#' roy_dynamic_effects(fit)
#' @export
roy_dynamic_effects <- function(fit) {
  stopifnot(inherits(fit, "roy_dynamic"))
  th <- fit$theta; B <- fit$beliefs; S <- fit$S; Y <- fit$Y; T <- fit$T
  cvec <- drop(fit$Xc %*% th$cost); mt <- B$M - cvec
  sol <- roy_dynamic_solve(th, T, fit$beta, grid = .rd_grid(B, cvec, th, T))
  tt <- rep(seq_len(T), each = fit$N)
  dv <- matrix(.rd_dv_at(sol, tt, as.vector(B$N), as.vector(mt)), fit$N, T)
  data.frame(
    t = fit$times, share_ai = colMeans(S),
    ATE = colMeans(B$M),
    ATT = vapply(seq_len(T), function(t) mean(B$M[S[, t] == 1, t]), numeric(1)),
    ATU = vapply(seq_len(T), function(t) mean(B$M[S[, t] == 0, t]), numeric(1)),
    naive_gap = vapply(seq_len(T), function(t)
      mean(Y[S[, t] == 1, t]) - mean(Y[S[, t] == 0, t]), numeric(1)),
    option_value_adopters = vapply(seq_len(T), function(t)
      mean(mt[S[, t] == 1, t] < 0), numeric(1)),
    p_dynamic = colMeans(stats::plogis(dv / th$lambda)),
    p_myopic = colMeans(stats::plogis(mt / th$lambda)))
}

#' Simulate a panel from Model E
#'
#' Simulates `N` workers over `T` periods from the dynamic Roy model with
#' Kalman learning and learning by doing, solving the model by
#' [roy_dynamic_solve()]. The cost of AI use is `cost[1] + cost[2] * x` with
#' `x ~ N(0, 1)` a time-invariant worker covariate.
#' @param N workers
#' @param T periods
#' @param mu0,sigma0 human-only earnings mean and sd
#' @param mu_k,sd_k mean and sd of initial complementarity `kappa_i1`
#' @param sigma_e noise sd of AI earnings (signal noise)
#' @param psi learning-by-doing rate
#' @param cost cost coefficients on (1, x)
#' @param lambda logit scale
#' @param beta discount factor
#' @return long data.frame with `id`, `t`, `s`, `y`, `x`, potential outcomes
#'   `y0`, `y1`, true `kappa` and beliefs `m`, `n`
#' @examples
#' set.seed(1)
#' p <- roy_dynamic_sim(200, T = 5)
#' head(p)
#' tapply(p$s, p$t, mean)   # adoption rate by period
#' @export
roy_dynamic_sim <- function(N, T = 10, mu0 = 1, sigma0 = 0.3, mu_k = 0, sd_k = 0.25,
                            sigma_e = 0.3, psi = 0.05, cost = c(0.1, 0.1),
                            lambda = 0.1, beta = 0.95) {
  par <- list(mu0 = mu0, sigma0 = sigma0, mu_k = mu_k, sd_k = sd_k, sigma_e = sigma_e,
              psi = psi, cost = cost, lambda = lambda)
  x <- stats::rnorm(N); cvec <- cost[1] + cost[2] * x
  sol <- roy_dynamic_solve(par, T, beta,
                           grid = seq(mu_k - max(cvec) - 6 * sd_k - 1,
                                      mu_k - min(cvec) + 6 * sd_k + abs(psi) * T + 1,
                                      length.out = 201))
  kap <- stats::rnorm(N, mu_k, sd_k); m <- rep(mu_k, N); n <- integer(N)
  out <- vector("list", T)
  for (t in seq_len(T)) {
    dv <- .rd_dv_at(sol, rep(t, N), n, m - cvec)
    s <- as.integer(stats::runif(N) < stats::plogis(dv / lambda))
    y0 <- mu0 + stats::rnorm(N, 0, sigma0)
    y1 <- mu0 + kap + stats::rnorm(N, 0, sigma_e)
    out[[t]] <- data.frame(id = seq_len(N), t = t, s = s, y = ifelse(s == 1, y1, y0),
                           x = x, y0 = y0, y1 = y1, kappa = kap, m = m, n = n)
    vs <- .rd_postvar(n, sd_k, sigma_e); k <- vs / (vs + sigma_e^2)
    ai <- s == 1
    m[ai] <- m[ai] + k[ai] * (y1[ai] - mu0 - m[ai]) + psi
    kap[ai] <- kap[ai] + psi
    n[ai] <- n[ai] + 1L
  }
  d <- do.call(rbind, out)
  d[order(d$id, d$t), , drop = FALSE]
}
