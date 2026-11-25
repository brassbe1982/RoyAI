# ---------------------------------------------------------------------------
# roy_spec: model specification shared by every estimator and plug-in
# ---------------------------------------------------------------------------

#' Specify a generalized Roy model
#'
#' Builds the design matrices shared by all RoyAI estimators. The number of
#' distinct values of the selection variable decides the model: two values
#' give Model A (binary AI adoption), three or more give Model B
#' (multinomial production arrangements; by convention 0 = human-only,
#' 1 = AI-augmented, 2 = AI-dominant).
#'
#' @param outcome Formula for the potential outcomes, e.g. `y ~ x1 + x2`.
#'   The same regressors enter every regime equation (\eqn{Y_s = X'\beta_s + U_s}).
#' @param selection Formula for selection, e.g. `d ~ x1 + x2 + z1`. Variables
#'   that appear here but not in `outcome` are the exclusion restrictions.
#' @param data A data frame.
#' @param regimes Optional labels for the regimes, in order.
#' @return An object of class `roy_spec`.
#' @export
roy_spec <- function(outcome, selection, data, regimes = NULL) {
  stopifnot(inherits(outcome, "formula"), inherits(selection, "formula"),
            is.data.frame(data))
  vars <- unique(c(all.vars(outcome), all.vars(selection)))
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("variables not found in data: ", paste(miss, collapse = ", "))
  cc <- stats::complete.cases(data[, vars, drop = FALSE])
  if (any(!cc)) message(sum(!cc), " row(s) with missing values dropped")
  data <- data[cc, , drop = FALSE]

  mf_y <- stats::model.frame(outcome, data)
  y <- stats::model.response(mf_y)
  X <- stats::model.matrix(outcome, mf_y)
  mf_d <- stats::model.frame(selection, data)
  draw <- stats::model.response(mf_d)
  Z <- stats::model.matrix(selection, mf_d)

  if (is.factor(draw)) {
    lev <- levels(draw); d <- as.integer(draw) - 1L
  } else if (is.logical(draw)) {
    lev <- c("FALSE", "TRUE"); d <- as.integer(draw)
  } else {
    lev <- sort(unique(draw)); d <- match(draw, lev) - 1L
    lev <- as.character(lev)
  }
  J <- length(lev)
  if (J < 2) stop("the selection variable has a single value")
  if (is.null(regimes)) {
    regimes <- if (J == 3) c("human_only", "ai_augmented", "ai_dominant")
               else if (J == 2) c("no_ai", "ai_adopter") else lev
  }
  stopifnot(length(regimes) == J)
  tab <- tabulate(d + 1L, nbins = J)
  if (any(tab < ncol(X) + 2))
    stop("too few observations in regime(s): ", paste(regimes[tab < ncol(X) + 2], collapse = ", "))

  excl <- setdiff(colnames(Z), colnames(X))
  if (!length(excl))
    warning("no exclusion restriction: every selection regressor also enters the outcome; ",
            "identification then rests on functional form only", call. = FALSE)

  structure(list(
    y = as.numeric(y), d = d, X = X, Z = Z, J = J, levels = lev, regimes = regimes,
    model = if (J == 2) "A" else "B",
    outcome = outcome, selection = selection,
    terms_x = stats::delete.response(stats::terms(mf_y)),
    terms_z = stats::delete.response(stats::terms(mf_d)),
    xlev_x = stats::.getXlevels(stats::terms(mf_y), mf_y),
    xlev_z = stats::.getXlevels(stats::terms(mf_d), mf_d),
    excluded = excl, n = length(d), data = data),
    class = "roy_spec")
}

#' @export
print.roy_spec <- function(x, ...) {
  cat("<roy_spec> Model", x$model, "|", x$n, "obs |", x$J, "regimes:",
      paste(sprintf("%s (n=%d)", x$regimes, tabulate(x$d + 1L, x$J)), collapse = ", "), "\n")
  cat("  outcome:  ", deparse(x$outcome), "\n")
  cat("  selection:", deparse(x$selection), "\n")
  cat("  excluded from outcome:", if (length(x$excluded)) paste(x$excluded, collapse = ", ") else "<none>", "\n")
  invisible(x)
}

# Design matrices for new data (used by predict / roy_potential)
.design_new <- function(spec, newdata) {
  if (is.null(newdata)) return(list(X = spec$X, Z = spec$Z))
  X <- stats::model.matrix(spec$terms_x, stats::model.frame(spec$terms_x, newdata, xlev = spec$xlev_x))
  Z <- stats::model.matrix(spec$terms_z, stats::model.frame(spec$terms_z, newdata, xlev = spec$xlev_z))
  list(X = X, Z = Z)
}
