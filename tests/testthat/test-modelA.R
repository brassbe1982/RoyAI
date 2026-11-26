dat <- roy_sim_A(4000, seed = 11)
tr <- attr(dat, "truth")

test_that("ML recovers Model A parameters", {
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
  expect_s3_class(f, "roy_fit_A")
  expect_true(f$converged)
  se <- sqrt(diag(vcov(f)))[names(tr)]
  expect_true(all(abs(f$par[names(tr)] - tr) < 4 * se))
})

test_that("two-step agrees with ML and has corrected standard errors", {
  f1 <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
  f2 <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, method = "twostep")
  b <- grep("^y[01]:", names(f1$par), value = TRUE)
  expect_equal(f2$par[b], f1$par[b], tolerance = 0.1)
  expect_true(all(is.finite(diag(f2$vcov)[b])))
})

test_that("Gaussian family via the generic h-function path matches the fast path", {
  g <- roy_dep_gaussian(); g$name <- "gaussian_generic"
  roy_register_dependence(g, overwrite = TRUE)
  small <- dat[1:1500, ]
  f1 <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = small)
  f2 <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = small, dependence = "gaussian_generic")
  expect_equal(f1$loglik, f2$loglik, tolerance = 1e-4)
  # numerical-integration potential outcomes agree with the closed form
  co <- f2$coef; co$cov_sd <- NULL
  q <- roy_potential(f2, 1, 0, coef_list = co)
  expect_equal(q, roy_potential(f1, 1, 0), tolerance = 1e-3)
})

test_that("effects are close to the simulated truth", {
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
  e <- roy_effects(f)
  truth <- c(mean(dat$y1 - dat$y0), mean((dat$y1 - dat$y0)[dat$d == 0]), mean((dat$y1 - dat$y0)[dat$d == 1]))
  expect_true(all(abs(e$estimate - truth) < 4 * e$se + 0.05))
})

test_that("independence family equals regime-wise OLS plus probit", {
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, dependence = "independence")
  o1 <- coef(lm(y ~ x1 + x2, dat, subset = d == 1))
  expect_equal(unname(f$par[c("y1:(Intercept)", "y1:x1", "y1:x2")]), unname(o1), tolerance = 1e-3)
})

test_that("classical Roy is a restriction of the generalized model", {
  dc <- roy_sim_classical(3000, seed = 5)
  fc <- roy(y ~ x1 + x2, d ~ x1 + x2, data = dc, type = "classical")
  fg <- suppressWarnings(roy(y ~ x1 + x2, d ~ x1 + x2, data = dc))
  expect_lte(fc$loglik, fg$loglik + 1e-6)
  lr <- roy_lrtest(fc, fg)
  expect_gt(lr$p.value, 0.001)
  b <- grep("^y[01]:", names(fc$par), value = TRUE)
  tc <- attr(dc, "truth")
  expect_true(all(abs(fc$par[b] - tc[b]) < 4 * sqrt(diag(fc$vcov)[b])))
  # extended Roy with an observed cost shifter
  de <- roy_sim_classical(3000, pi = 0.6, seed = 6)
  fe <- roy(y ~ x1 + x2, d ~ x1 + x2 + w, data = de, type = "extended")
  expect_equal(unname(fe$par["cost:w"]), 0.6, tolerance = 0.2)
})

test_that("diagnostics run and flag missing exclusion restrictions", {
  skip_on_cran()  # slow (> 3 s); runs when NOT_CRAN = "true"
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
  dg <- roy_diagnostics(f)
  expect_s3_class(dg, "roy_diagnostics")
  expect_lt(dg$instruments$p_value, 1e-6)
  expect_lt(dg$selection$p_value, 0.01)
  f0 <- suppressWarnings(roy(y ~ x1 + x2, d ~ x1 + x2, data = dat))
  expect_true(any(grepl("exclusion", roy_diagnostics(f0)$flags)))
  expect_output(print(dg), "RoyAI diagnostics")
  expect_output(print(summary(f)), "Model A")
})

test_that("predict and propensity have the right shape", {
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, method = "twostep")
  expect_equal(dim(predict(f)), c(nrow(dat), 2))
  expect_length(predict(f, type = "observed"), nrow(dat))
  nd <- dat[1:5, ]
  expect_equal(dim(predict(f, nd, type = "propensity")), c(5, 2))
})
