# =============================================================================
# Monte Carlo study for Model C (copula generalized Roy model) -- Paper 1
# -----------------------------------------------------------------------------
# Usage:  Rscript mc_copula.R <design> [cores]
#   design = main_small | main_large | misspec | j3 | weakiv | all
# Raw replications are checkpointed to results/raw/<cell>.rds, so an
# interrupted run resumes where it stopped.  Tables: mc_tables.R
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
design <- if (exists("MC_LOAD_ONLY")) "none" else if (length(args)) args[1] else "main_small"
cores  <- if (length(args) > 1) as.integer(args[2]) else parallel::detectCores()

here <- if (exists("MC_HERE")) MC_HERE else local({
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
pkg_root <- normalizePath(file.path(here, "..", ".."))
source(file.path(pkg_root, "R", "roy_copula.R"))
source(file.path(pkg_root, "R", "roy_mte.R"))
out_dir <- file.path(here, "results", "raw"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG")

# ---- DGP constants (state these in the paper) --------------------------------
DGP <- list(
  beta  = list(c(0.5, 0.3), c(1.0, 0.5)),        # Y0, Y1 on (1, x1)
  gamma = c(0.2, 0.3, 0.8, 0.6),                 # selection on (1, x1, z1, z2)
  gamma_weak = c(0.2, 0.3, 0.15, 0.10),          # weak-instrument variant
  sigma = c(0.8, 1.0),
  nu = 5,                                        # t-copula df
  beta3 = list(c(0.5, 0.3), c(1.0, 0.5), c(1.3, 0.8)),
  sigma3 = c(0.8, 1.0, 1.2), kappa3 = 1.2)
# x1 ~ N(0,1) (in X and Z), z1 ~ N(0,1) and z2 ~ Bernoulli(0.5) excluded from X.
# LATE: z2 switched from 0 to 1 ("AI access").
# (tau_0, tau_1) = (-tau, +tau): Roy sorting on comparative advantage.
FAMS <- c("gaussian", "t", "clayton", "gumbel", "frank")
UGRID <- seq(0.05, 0.95, by = 0.05)

draw_cov <- function(n) data.frame(x1 = stats::rnorm(n), z1 = stats::rnorm(n),
                                   z2 = stats::rbinom(n, 1, 0.5))

truth_cache <- new.env()
get_truth <- function(fam, tau, gamma = DGP$gamma, J = 2) {
  key <- paste(fam, tau, paste(gamma, collapse = "_"), J)
  if (!is.null(truth_cache[[key]])) return(truth_cache[[key]])
  set.seed(20261124)
  cov <- draw_cov(1e6)
  if (J == 2) {
    tv <- roy_te_truth(DGP$beta, gamma, DGP$sigma, c(-tau, tau), fam, nu = DGP$nu,
                       covariates = cov, instrument = "z2")
    th <- list(beta = DGP$beta, gamma = gamma, kappa = c(-Inf, 0, Inf),
               sigma = DGP$sigma, tau = c(-tau, tau), nu = DGP$nu)
    mg <- .roy_mgrid(th, fam, 2)
    xb <- sum(c(1, 0) * (DGP$beta[[2]] - DGP$beta[[1]]))      # at x1 = 0
    mte <- xb + .interp_m(mg[[2]], 1 - UGRID) - .interp_m(mg[[1]], 1 - UGRID)
    res <- list(te = tv, mte = mte,
                par = c(unlist(DGP$beta), gamma, DGP$sigma, -tau, tau,
                        if (fam == "t") DGP$nu))
  } else {
    taus <- c(-tau, tau / 2, tau)
    t10 <- roy_te_truth(DGP$beta3, gamma, DGP$sigma3, taus, fam, nu = DGP$nu,
                        kappa = DGP$kappa3, covariates = cov, to = 1, from = 0)
    t20 <- roy_te_truth(DGP$beta3, gamma, DGP$sigma3, taus, fam, nu = DGP$nu,
                        kappa = DGP$kappa3, covariates = cov, to = 2, from = 0)
    res <- list(te = c(ATE10 = t10[["ATE"]], ATT10 = t10[["TT_s1"]],
                       ATE20 = t20[["ATE"]], ATT20 = t20[["TT_s2"]]),
                par = c(unlist(DGP$beta3), gamma, DGP$kappa3, DGP$sigma3, taus,
                        if (fam == "t") DGP$nu))
  }
  truth_cache[[key]] <- res
  res
}

# ---- estimators ------------------------------------------------------------
# Treatment effects (and LATE) for a parameter vector on the sample.
te_fun <- function(par, sp, X, Z, Z0, Z1) {
  th <- .roy_unpack(par, sp)
  mg <- .roy_mgrid(th, sp$family, 2)
  zg <- drop(Z %*% th$gamma); P <- stats::pnorm(zg)
  xdb <- drop(X %*% (th$beta[[2]] - th$beta[[1]]))
  Gf <- function(P) -.interp_M(mg[[2]], 1 - P) + .interp_M(mg[[1]], 1 - P)
  P0 <- stats::pnorm(drop(Z0 %*% th$gamma)); P1 <- stats::pnorm(drop(Z1 %*% th$gamma))
  c(ATE = mean(xdb),
    ATT = sum(P * xdb + Gf(P)) / sum(P),
    ATU = sum((1 - P) * xdb - Gf(P)) / sum(1 - P),
    LATE = sum((P1 - P0) * xdb + Gf(P1) - Gf(P0)) / sum(P1 - P0))
}

mte_fun <- function(par, sp) {
  th <- .roy_unpack(par, sp); mg <- .roy_mgrid(th, sp$family, sp$J)
  sum(c(1, 0) * (th$beta[[2]] - th$beta[[1]])) +
    .interp_m(mg[[2]], 1 - UGRID) - .interp_m(mg[[1]], 1 - UGRID)
}

ols_est <- function(d, X) {
  out <- list(); V <- list()
  for (j in 0:1) {
    idx <- d$s == j
    f <- stats::lm.fit(X[idx, ], d$y[idx])
    out[[j + 1]] <- f$coefficients
    s2 <- sum(f$residuals^2) / (sum(idx) - ncol(X))
    V[[j + 1]] <- s2 * chol2inv(qr.R(f$qr))
  }
  xb <- colMeans(X); x1 <- colMeans(X[d$s == 1, ]); x0 <- colMeans(X[d$s == 0, ])
  db <- out[[2]] - out[[1]]; Vd <- V[[1]] + V[[2]]
  list(beta = c(out[[1]], out[[2]]),
       te = c(ATE = sum(xb * db), ATT = sum(x1 * db), ATU = sum(x0 * db), LATE = NA),
       te_se = sqrt(c(drop(t(xb) %*% Vd %*% xb), drop(t(x1) %*% Vd %*% x1),
                      drop(t(x0) %*% Vd %*% x0), NA)))
}

one_rep <- function(n, fam, tau, gamma, fit_fams = fam, hess_fams = fam,
                    want_te = TRUE) {
  cov <- draw_cov(n)
  d <- roy_copula_sim(n, DGP$beta, gamma, DGP$sigma, c(-tau, tau), fam,
                      nu = DGP$nu, covariates = cov)
  X <- cbind(1, d$x1); Z <- cbind(1, d$x1, d$z1, d$z2)
  Z0 <- Z; Z0[, 4] <- 0; Z1 <- Z; Z1[, 4] <- 1
  res <- list(share1 = mean(d$s))
  # OLS by regime and normal two-step
  res$ols <- ols_est(d, X)
  spg <- list(kx = list(k = 2, names = c("(Intercept)", "x1")),
              kz = list(k = 4, names = c("(Intercept)", "x1", "z1", "z2")),
              J = 2, family = "gaussian", fix_nu = NULL, type = "ordered", link = "probit")
  ts <- roy_twostep(d$y, X, Z, d$s, spg)
  res$twostep <- list(par = .roy_natural(ts, spg),
                      te = if (want_te) te_fun(ts, spg, X, Z, Z0, Z1))
  # Copula FIML
  res$fiml <- list()
  for (f in fit_fams) {
    hs <- f %in% hess_fams
    fit <- tryCatch(roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = f,
                                   hessian = hs), error = function(e) NULL)
    if (is.null(fit)) { res$fiml[[f]] <- NULL; next }
    r <- list(par = coef(fit), loglik = fit$loglik, npar = fit$npar,
              conv = fit$convergence)
    if (hs) {
      r$se <- sqrt(pmax(diag(vcov(fit)), 0))
      if (want_te) {
        g <- function(p) te_fun(p, fit$sp, X, Z, Z0, Z1)
        r$te <- g(fit$par)
        G <- numDeriv::jacobian(g, fit$par, method = "simple")
        r$te_se <- sqrt(pmax(diag(G %*% fit$vcov_raw %*% t(G)), 0))
        r$mte <- mte_fun(fit$par, fit$sp)
      }
    }
    res$fiml[[f]] <- r
  }
  res
}

one_rep_j3 <- function(n, fam, tau) {
  cov <- draw_cov(n)
  taus <- c(-tau, tau / 2, tau)
  d <- roy_copula_sim(n, DGP$beta3, DGP$gamma, DGP$sigma3, taus, fam, nu = DGP$nu,
                      kappa = DGP$kappa3, covariates = cov)
  fit <- tryCatch(roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = fam),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  t10 <- roy_te(fit, to = 1, from = 0); t20 <- roy_te(fit, to = 2, from = 0)
  list(shares = tabulate(d$s + 1, 3) / n, par = coef(fit),
       se = sqrt(pmax(diag(vcov(fit)), 0)), conv = fit$convergence,
       te = c(ATE10 = t10$estimate[1], ATT10 = t10$estimate[3],
              ATE20 = t20$estimate[1], ATT20 = t20$estimate[4]),
       te_se = c(t10$se[1], t10$se[3], t20$se[1], t20$se[4]))
}

run_cell <- function(cell, R, FUN) {
  file <- file.path(out_dir, paste0(cell, ".rds"))
  if (file.exists(file)) { message("skip ", cell); return(invisible()) }
  t0 <- Sys.time()
  seed <- sum(utf8ToInt(cell) * seq_along(utf8ToInt(cell)))
  set.seed(seed)
  reps <- parallel::mclapply(seq_len(R), function(r) tryCatch(FUN(), error = function(e) NULL),
                             mc.cores = cores, mc.set.seed = TRUE)
  saveRDS(list(cell = cell, R = R, reps = reps, seed = seed,
               time = as.numeric(Sys.time() - t0, units = "secs")), file)
  mirror <- Sys.getenv("MC_MIRROR")          # optional durable copy of each cell
  if (nzchar(mirror)) file.copy(file, file.path(mirror, basename(file)), overwrite = TRUE)
  message(sprintf("%s done: R=%d, %.0f s, %d failed", cell, R,
                  as.numeric(Sys.time() - t0, units = "secs"),
                  sum(vapply(reps, is.null, TRUE))))
}

# ---- designs -----------------------------------------------------------------
taus <- c(0, 0.3, 0.6)
if (design %in% c("main_small", "main_n2000", "main_n500", "all")) {
  ns <- switch(design, main_n2000 = 2000, main_n500 = 500, c(2000, 500))
  for (n in ns) for (tau in rev(taus)) for (f in FAMS) {
    R <- if (f == "t") 250 else 500
    run_cell(sprintf("main_%s_tau%.1f_n%d", f, tau, n), R,
             function() one_rep(n, f, tau, DGP$gamma))
  }
}
if (design %in% c("misspec", "all")) {
  # every family fitted to every DGP: AIC/BIC hit rates + Gaussian misspecification
  for (tau in c(0.6, 0.3)) for (f in FAMS) {
    run_cell(sprintf("misspec_%s_tau%.1f_n2000", f, tau), 100,
             function() one_rep(2000, f, tau, DGP$gamma, fit_fams = FAMS,
                                hess_fams = unique(c("gaussian", f))))
  }
}
if (design %in% c("weakiv", "all")) {
  for (f in c("gaussian", "clayton", "gumbel")) {
    run_cell(sprintf("weakiv_%s_tau0.6_n2000", f), 300,
             function() one_rep(2000, f, 0.6, DGP$gamma_weak))
  }
}
if (design %in% c("j3", "all")) {
  for (f in FAMS) {
    run_cell(sprintf("j3_%s_tau0.6_n3000", f), if (f == "t") 100 else 200,
             function() one_rep_j3(3000, f, 0.6))
  }
}
if (design %in% c("main_large", "all")) {
  for (tau in rev(taus)) for (f in FAMS) {
    R <- if (f == "t") 50 else 100
    run_cell(sprintf("main_%s_tau%.1f_n10000", f, tau), R,
             function() one_rep(10000, f, tau, DGP$gamma))
  }
}
