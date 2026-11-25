test_that("DMF terms equal the conditional means of Gumbel errors", {
  set.seed(1); n <- 4e5; v <- c(0, 0.5, -0.3)
  eps <- -log(-log(matrix(runif(n * 3), n)))
  s <- max.col(sweep(eps, 2, v, "+")) - 1L
  P <- matrix(exp(v) / sum(exp(v)), 1)
  for (k in 0:2) {
    emp <- colMeans(eps[s == k, ]) + digamma(1)
    expect_equal(as.vector(roy_dmf_terms(P, k)), emp, tolerance = 0.02)
  }
})

test_that("probit and multinomial logit stages recover parameters", {
  dat <- roy_sim_B(6000, seed = 2)
  sp <- roy_spec(y ~ x1 + x2, s ~ x1 + x2 + z1 + z2, dat)
  fs <- roy_select(sp)
  tr <- attr(dat, "truth")
  expect_equal(as.vector(fs$coef), unname(tr[grep("^sel", names(tr))]), tolerance = 0.15)
  expect_equal(unname(rowSums(fs$prob)), rep(1, sp$n))
  da <- roy_sim_A(4000, seed = 3)
  pa <- roy_select(roy_spec(y ~ x1 + x2, d ~ x1 + x2 + z, da))
  expect_equal(unname(pa$coef), c(-0.2, 0.4, 0.2, 0.8), tolerance = 0.15)
})

test_that("inverse Mills ratios have the right sign", {
  expect_true(all(roy_imr(c(-1, 0, 1), 1) > 0))
  expect_true(all(roy_imr(c(-1, 0, 1), 0) < 0))
})
