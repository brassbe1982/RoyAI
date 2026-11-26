test_that("h-functions are proper conditional CDFs for every family", {
  u <- (1:4000 - 0.5) / 4000
  for (f in roy_copula_families) for (tau in c(-0.5, 0.3, 0.7)) {
    cp <- copula_par(f, tau, nu = 5)
    # E_u[h(v | u)] = v  (uniform margins)
    expect_equal(mean(copula_hfunc(0.3, u, cp)), 0.3, tolerance = 1e-3)
    # monotone in v
    h <- copula_hfunc(c(0.1, 0.5, 0.9), rep(0.4, 3), cp)
    expect_true(all(diff(h) > 0))
  }
})

test_that("Kendall tau is reproduced by simulation", {
  set.seed(1)
  for (f in roy_copula_families) {
    cp <- copula_par(f, -0.4, nu = 5)
    w <- runif(3000); u <- RoyAI:::.hinv2(runif(3000), w, cp)
    expect_lt(abs(cor(u, w, method = "kendall") + 0.4), 0.05)
  }
})

test_that("Gaussian conditional mean matches the closed form", {
  cp <- copula_par("gaussian", 0.4); w <- c(0.01, 0.3, 0.5, 0.95)
  expect_equal(copula_cond_mean(w, cp, 1.3), 1.3 * cp$theta * qnorm(w), tolerance = 1e-6)
})

test_that("analytic ATE/ATT/ATU/LATE match a large simulation", {
  skip_on_cran()  # slow (> 3 s); runs when NOT_CRAN = "true"
  set.seed(5)
  B <- list(c(0.5, 0.3), c(1, 0.5)); G <- c(0.2, 0.3, 0.8, 0.6)
  n <- 2e5
  cov <- data.frame(x1 = rnorm(n), z1 = rnorm(n), z2 = rbinom(n, 1, 0.5))
  for (f in c("clayton", "gumbel")) {
    tr <- roy_te_truth(B, G, c(0.8, 1), c(-0.3, 0.4), f, covariates = cov,
                       instrument = "z2")
    d <- roy_copula_sim(n, B, G, c(0.8, 1), c(-0.3, 0.4), f, covariates = cov)
    g <- d$y1 - d$y0
    expect_lt(abs(tr[["ATT"]] - mean(g[d$s == 1])), 0.02)
    expect_lt(abs(tr[["ATU"]] - mean(g[d$s == 0])), 0.02)
  }
})

test_that("FIML recovers Model C parameters (binary and ordered J = 3)", {
  skip_on_cran()  # slow (> 3 s); runs when NOT_CRAN = "true"
  set.seed(2)
  d <- roy_copula_sim(3000, list(c(0.5, 0.3), c(1, 0.5)), c(0.2, 0.3, 0.8, 0.6),
                      c(0.8, 1), c(-0.4, 0.5), "clayton")
  fit <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "clayton")
  b <- coef(fit)
  expect_equal(unname(b[c("tau0", "tau1")]), c(-0.4, 0.5), tolerance = 0.1)
  expect_equal(unname(b[c("beta1:(Intercept)", "beta1:x1")]), c(1, 0.5), tolerance = 0.1)
  te <- roy_te(fit)
  expect_equal(te$parameter, c("ATE", "ATT", "ATU"))
  expect_true(all(te$se > 0))
  m <- roy_mte(fit, u = c(0.1, 0.5, 0.9))
  expect_true(m$mte[1] > m$mte[3])           # Roy sorting: gains fall with resistance
  la <- roy_late(fit, "z2", data = d)
  expect_true(is.finite(la$estimate[1]))

  d3 <- roy_copula_sim(3000, list(c(0.5, 0.3), c(1, 0.5), c(1.3, 0.8)),
                       c(0.2, 0.3, 0.8, 0.6), c(0.8, 1, 1.2), c(-0.3, 0.3, 0.5),
                       "frank", kappa = 1.2)
  f3 <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d3, family = "frank")
  expect_equal(unname(coef(f3)["kappa2"]), 1.2, tolerance = 0.1)
  expect_equal(nrow(roy_te(f3, to = 2, from = 0)), 4)
})

test_that("copula families plug into the core roy() estimator", {
  skip_on_cran()  # slow (> 3 s); runs when NOT_CRAN = "true"
  skip_if_not(exists("roy", mode = "function"))
  expect_true(all(c("clayton", "gumbel", "frank", "t") %in% roy_dependence_families()))
  set.seed(3)
  d <- roy_copula_sim(2000, list(c(0.5, 0.3), c(1, 0.5)), c(0.2, 0.3, 0.8, 0.6),
                      c(0.8, 1), c(0.3, 0.5), "gumbel")
  d$d <- d$s
  f <- roy(y ~ x1, d ~ x1 + z1 + z2, data = d, dependence = "gumbel")
  own <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "gumbel")
  expect_equal(f$loglik, own$loglik, tolerance = 1e-3)
  expect_equal(roy_te(f)$estimate, roy_te(own)$estimate, tolerance = 1e-2)
  expect_equal(nrow(roy_mte(f, u = c(0.2, 0.8))), 2)
})

test_that("multistart rejects an irregular boundary maximum (Frank on Clayton data)", {
  skip_on_cran()  # slow (> 3 s); runs when NOT_CRAN = "true"
  # Clayton DGP, (tau0, tau1) = (-0.3, 0.7). From the flipped-tau start the Frank
  # likelihood climbs to tau1 -> 1 with a non-PD Hessian and a spuriously high
  # log-likelihood; roy() must return the regular maximum instead.
  set.seed(1)
  d <- roy_copula_sim(3000, beta = list(c(0.5, 0.3), c(1, 0.8)), gamma = c(0.2, 0.3, 0.7, 0.7),
                      sigma = c(0.8, 1), tau = c(-0.3, 0.7), family = "clayton")
  fit <- suppressWarnings(roy(y ~ x1, s ~ x1 + z1 + z2, data = d, dependence = "frank"))
  expect_true("regular" %in% names(fit$multistart))
  expect_lt(abs(fit$par[["tau1"]]), 0.95)
  ev <- eigen(fit$optim$hessian, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), 0)
  expect_true(isTRUE(fit$multistart$regular[fit$multistart$chosen]))
})

test_that("roy_copula_fit rejects an edge start and keeps the best regular maximum", {
  skip_on_cran()  # slow (> 3 s); runs when NOT_CRAN = "true"
  set.seed(1)
  d <- roy_copula_sim(3000, beta = list(c(0.5, 0.3), c(1, 0.8)), gamma = c(0.2, 0.3, 0.7, 0.7),
                      sigma = c(0.8, 1), tau = c(-0.3, 0.7), family = "clayton")
  ref <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "frank")
  expect_setequal(ref$multistart$start, c("twostep", "independence", "flipped_tau"))
  expect_true(isTRUE(ref$multistart$regular[ref$multistart$chosen]))
  # a start pushed to the boundary (tau1 = 0.999) must not be the one kept
  st <- ref$par; st["atanh_tau1"] <- atanh(0.999)
  fit <- suppressWarnings(roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "frank",
                                         start = st, multistart = TRUE))
  expect_true(all(abs(fit$theta$tau) <= 0.98))
  expect_gt(min(eigen(solve(fit$vcov_raw), symmetric = TRUE, only.values = TRUE)$values), 0)
  # Gaussian family: single start by default
  g <- roy_copula_fit(y ~ x1, s ~ x1 + z1 + z2, d, family = "gaussian")
  expect_equal(nrow(g$multistart), 1)
})
