test_that("copula log-densities integrate the h-function", {
  for (f in roy_copula_families) for (tau in c(-0.5, 0.3, 0.7)) {
    cp <- copula_par(f, tau, nu = 5); u <- c(0.1, 0.4, 0.8); v <- c(0.3, 0.6, 0.9)
    num <- (copula_hfunc(v + 1e-6, u, cp) - copula_hfunc(v - 1e-6, u, cp)) / 2e-6
    expect_equal(exp(copula_logdens(u, v, cp)), num, tolerance = 1e-5)
  }
})

test_that("Model D simulator produces censored intensity in [0, 1)", {
  set.seed(1)
  d <- roy_continuous_sim(3000, tau = 0.4)
  expect_true(all(d$a >= 0 & d$a < 1))
  expect_gt(mean(d$a == 0), 0.15); expect_lt(mean(d$a == 0), 0.45)
  expect_lt(abs(cor(d$eps, d$v, method = "kendall") - 0.4), 0.05)
})

test_that("Model D Gaussian ML recovers the parameters and effects", {
  set.seed(11)
  d <- roy_continuous_sim(3000, tau = 0.4)
  fit <- roy_continuous(y ~ x1, a ~ x1 + z1 + z2, d)
  expect_true(fit$regular)
  b <- coef(fit)
  expect_lt(abs(b[["tau"]] - 0.4), 0.06)
  expect_lt(abs(b[["gamma:z1"]] - 0.8), 0.08)
  expect_lt(abs(b[["alpha:a"]] - 0.5), 0.15)
  e <- roy_continuous_effects(fit, a = 0.5, grid = c(0, 0.5), se = TRUE)
  tr <- roy_continuous_truth(tau = 0.4, covariates = d, a = 0.5, grid = c(0, 0.5))
  z <- (e$estimate - tr$truth) / e$se
  expect_true(all(abs(z) < 4))
  expect_s3_class(summary(fit), "summary.roy_continuous")
  expect_equal(nobs(fit), 3000)
})

test_that("Model D non-Gaussian fit runs the multistart safeguard", {
  set.seed(4)
  d <- roy_continuous_sim(1500, tau = 0.4, family = "frank")
  fit <- roy_continuous(y ~ x1, a ~ x1 + z1 + z2, d, family = "frank", scale = ~ 1)
  expect_equal(nrow(fit$multistart), 3)
  expect_true(any(fit$multistart$chosen))
  expect_gt(coef(fit)[["tau"]], 0.2)
})

test_that("Model D validates its inputs", {
  d <- data.frame(y = rnorm(10), a = c(rep(0, 5), rep(1.2, 5)), z = rnorm(10))
  expect_error(roy_continuous(y ~ 1, a ~ z, d), "\\[0, 1\\)")
})
