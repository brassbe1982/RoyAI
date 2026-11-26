#' @keywords internal
#' @importFrom stats coef vcov logLik nobs predict
#' @rawNamespace export(roy_copula_families)
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  .register_builtin_dependence()
  # Extension modules shipped in the package (e.g. roy_copula) register
  # their dependence families and estimators here.
  ns <- environment(.onLoad)
  if (exists("roy_copula_register", envir = ns, inherits = FALSE)) {
    tryCatch(get("roy_copula_register", envir = ns)(), error = function(e)
      warning("RoyAI: copula registration failed: ", conditionMessage(e), call. = FALSE))
  }
  invisible()
}

#' Copula families shipped with roy_copula
#'
#' Character vector of the dependence families implemented in `roy_copula`
#' (registered for [roy()] at load time by `roy_copula_register()`).
#' @name roy_copula_families
#' @format A character vector.
#' @examples
#' roy_copula_families
NULL
