# ---------------------------------------------------------------------------
# Dependence families: the plug-in point for roy_copula
# ---------------------------------------------------------------------------
#
# In Model A every outcome disturbance U_s (s = 0, 1) is paired with the
# selection disturbance U_D through a bivariate dependence family. The core
# likelihood only ever needs the conditional distribution ("h-function")
#
#     h(v | u; theta) = P(V <= v | U = u),   U = F_s(U_s),  V = F_D(U_D)
#
# so any bivariate copula with a computable h-function can be plugged in
# without touching the estimator. Gaussian and independence families ship
# with the core; roy_copula registers further families (Clayton, Frank,
# Gumbel, Student-t, rotations, ...) via roy_register_dependence().

.roy_registry <- new.env(parent = emptyenv())
.roy_registry$dependence <- list()
.roy_registry$methods <- list()

#' Construct a dependence family
#'
#' @param name Family name, e.g. `"gaussian"`.
#' @param h Function `h(v, u, theta)` returning \eqn{P(V \le v \mid U = u)}
#'   for vectors of uniforms `u`, `v` and a parameter vector `theta`.
#' @param npar Number of dependence parameters per outcome-selection pair.
#' @param linkinv Map from unconstrained real vector `eta` to `theta`.
#' @param link Inverse of `linkinv`.
#' @param start Default starting value of `theta`.
#' @param tau Optional function `tau(theta)` giving Kendall's tau, used for
#'   reporting and for comparing families on a common scale.
#' @param cond_mean Optional closed form for \eqn{E[U_s \mid U_D > c]}; when
#'   absent the core integrates numerically (see [roy_potential()]).
#' @param par_names Names of the elements of `theta`.
#' @return An object of class `roy_dependence`.
#' @examples
#' # a copy of the Gaussian family under a new name, built from its pieces
#' g <- roy_dep_gaussian()
#' g2 <- roy_dependence("gaussian_copy", h = g$h, npar = 1, linkinv = g$linkinv,
#'                      link = g$link, start = g$start, tau = g$tau)
#' g2
#' @export
roy_dependence <- function(name, h, npar, linkinv, link, start,
                           tau = NULL, cond_mean = NULL,
                           par_names = if (npar == 1) "theta" else paste0("theta", seq_len(npar))) {
  stopifnot(is.character(name), length(name) == 1L, is.function(h),
            is.function(linkinv), is.function(link))
  structure(list(name = name, h = h, npar = as.integer(npar),
                 linkinv = linkinv, link = link, start = start, tau = tau,
                 cond_mean = cond_mean, par_names = par_names),
            class = "roy_dependence")
}

#' @export
print.roy_dependence <- function(x, ...) {
  cat("<roy_dependence>", x$name, "with", x$npar, "parameter(s):",
      paste(x$par_names, collapse = ", "), "\n")
  invisible(x)
}

#' Register or look up a dependence family
#'
#' `roy_register_dependence()` adds a family to the package registry so that
#' `roy(..., dependence = "name")` can use it. This is how the copula module
#' plugs into the core estimators.
#'
#' @param family A [roy_dependence()] object.
#' @param overwrite Replace an existing family of the same name?
#' @param name Family name to look up.
#' @return `roy_register_dependence()` returns the family invisibly;
#'   `roy_get_dependence()` returns the family; `roy_dependence_families()`
#'   returns the registered names.
#' @export
roy_register_dependence <- function(family, overwrite = FALSE) {
  stopifnot(inherits(family, "roy_dependence"))
  if (!overwrite && !is.null(.roy_registry$dependence[[family$name]]))
    stop("dependence family '", family$name, "' already registered; use overwrite = TRUE")
  .roy_registry$dependence[[family$name]] <- family
  invisible(family)
}

#' @rdname roy_register_dependence
#' @examples
#' g <- roy_dep_gaussian()
#' g$name <- "gaussian_copy"
#' roy_register_dependence(g, overwrite = TRUE)
#' "gaussian_copy" %in% roy_dependence_families()
#' roy_get_dependence("gaussian_copy")
#' @export
roy_get_dependence <- function(name) {
  if (inherits(name, "roy_dependence")) return(name)
  fam <- .roy_registry$dependence[[name]]
  if (is.null(fam))
    stop("unknown dependence family '", name, "'. Registered: ",
         paste(names(.roy_registry$dependence), collapse = ", "))
  fam
}

#' @rdname roy_register_dependence
#' @export
roy_dependence_families <- function() names(.roy_registry$dependence)

#' Gaussian dependence (the classical normal Roy / switching model)
#'
#' With Gaussian margins this reproduces the textbook normal endogenous
#' switching likelihood exactly; `theta` is the correlation \eqn{\rho}.
#' @return A [roy_dependence()] object.
#' @examples
#' g <- roy_dep_gaussian()
#' g
#' g$h(0.5, 0.3, 0.4)
#' @export
roy_dep_gaussian <- function() {
  roy_dependence(
    name = "gaussian", npar = 1L, par_names = "rho",
    h = function(v, u, theta) {
      rho <- theta[1]
      stats::pnorm((stats::qnorm(v) - rho * stats::qnorm(u)) / sqrt(1 - rho^2))
    },
    linkinv = function(eta) tanh(eta),
    link = function(theta) atanh(theta),
    start = 0,
    tau = function(theta) 2 / pi * asin(theta),
    # E[e_s | U_D > -c] for standard normal e_s, U_D: rho * phi(c) / Phi(c)
    cond_mean = function(c, theta, select = TRUE) {
      if (select) theta * stats::dnorm(c) / stats::pnorm(c)
      else -theta * stats::dnorm(c) / stats::pnorm(-c)
    }
  )
}

#' Independence (no selection on unobservables)
#'
#' Useful as the null model in diagnostics: with independence, selection is
#' on observables only and OLS within each regime is consistent.
#' @return A [roy_dependence()] object.
#' @examples
#' roy_dep_independence()
#' @export
roy_dep_independence <- function() {
  roy_dependence(
    name = "independence", npar = 0L, par_names = character(0),
    h = function(v, u, theta) v,
    linkinv = function(eta) numeric(0),
    link = function(theta) numeric(0),
    start = numeric(0),
    tau = function(theta) 0,
    cond_mean = function(c, theta, select = TRUE) rep(0, length(c))
  )
}

# Registration of built-ins happens in .onLoad (see zzz.R).
.register_builtin_dependence <- function() {
  roy_register_dependence(roy_dep_gaussian(), overwrite = TRUE)
  roy_register_dependence(roy_dep_independence(), overwrite = TRUE)
}

#' Register an estimation method
#'
#' Extension modules (e.g. a copula likelihood for Model B) can register an
#' estimator under a `(model, method)` key. The fitter receives a
#' [roy_spec()] object plus `...` and must return an object inheriting from
#' `roy_fit` (see `PLUGIN_CONTRACT.md` in the package sources).
#'
#' @param model `"A"` or `"B"`.
#' @param method Method name, e.g. `"copula_ml"`.
#' @param fitter Function `fitter(spec, ...)`.
#' @param overwrite Replace an existing registration?
#' @return The fitter, invisibly.
#' @export
roy_register_method <- function(model, method, fitter, overwrite = FALSE) {
  key <- paste(model, method, sep = ":")
  if (!overwrite && !is.null(.roy_registry$methods[[key]]))
    stop("method '", key, "' already registered; use overwrite = TRUE")
  .roy_registry$methods[[key]] <- fitter
  invisible(fitter)
}

#' @rdname roy_register_method
#' @examples
#' # a toy method that refits the specification with the built-in two-step
#' # estimator; a real plug-in would implement its own estimator here
#' my_twostep <- function(spec, ...)
#'   roy(spec$outcome, spec$selection, data = spec$data, method = "twostep")
#' roy_register_method("A", "my_twostep", my_twostep, overwrite = TRUE)
#' roy_methods()
#' dat <- roy_sim_A(n = 800, seed = 1)
#' fit <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, method = "my_twostep")
#' coef(fit)[1:3]
#' @export
roy_methods <- function() names(.roy_registry$methods)
