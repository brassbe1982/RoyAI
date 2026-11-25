test_that("Model E solver: value increases in beliefs and dv is monotone", {
  par <- list(sd_k = 0.25, sigma_e = 0.3, psi = 0.05, lambda = 0.1)
  sol <- roy_dynamic_solve(par, T = 5)
  expect_true(all(diff(sol$dv[1, 1, ]) > 0))
  expect_true(all(diff(sol$Vbar[1, 1, ]) >= -1e-10))
  # in the last period the decision is myopic
  expect_equal(sol$dv[5, 1, ], sol$grid)
  # option value: dv exceeds the current gain early on (Prop. experiment)
  mid <- which.min(abs(sol$grid))
  expect_gt(sol$dv[1, 1, mid], sol$grid[mid])
})

test_that("Kalman beliefs match the closed-form posterior when psi = 0", {
  set.seed(1)
  p <- roy_dynamic_sim(300, T = 4, psi = 0)
  last <- p[p$t == 4, ]
  par <- list(mu0 = 1, mu_k = 0, sd_k = 0.25, sigma_e = 0.3, psi = 0)
  S <- matrix(p$s, 300, 4, byrow = TRUE); Y <- matrix(p$y, 300, 4, byrow = TRUE)
  B <- RoyAI:::.rd_beliefs(S, Y, par)
  expect_equal(B$M[, 4], last$m, tolerance = 1e-10)
  i <- which(rowSums(S[, 1:3]) == 3)[1]
  post <- sum(Y[i, 1:3] - 1) / 0.3^2 / (1 / 0.25^2 + 3 / 0.3^2)
  expect_equal(B$M[i, 4], post, tolerance = 1e-10)
})

test_that("Model E CCP and NFXP recover the structural parameters", {
  set.seed(3)
  p <- roy_dynamic_sim(800, T = 8)
  f1 <- roy_dynamic(p, cost = ~ x, method = "ccp")
  b <- coef(f1)
  expect_lt(abs(b[["psi"]] - 0.05), 0.01)
  expect_lt(abs(b[["sigma_e"]] - 0.3), 0.02)
  expect_lt(abs(b[["lambda"]] - 0.1), 0.03)
  skip_on_cran()
  f2 <- roy_dynamic(p, cost = ~ x, method = "nfxp")
  b2 <- coef(f2); se <- sqrt(diag(vcov(f2)))
  truth <- c(mu0 = 1, sigma0 = 0.3, mu_k = 0, sd_k = 0.25, sigma_e = 0.3, psi = 0.05,
             `cost:(Intercept)` = 0.1, `cost:x` = 0.1, lambda = 0.1)
  expect_true(all(abs((b2 - truth) / se) < 4))
  ef <- roy_dynamic_effects(f2)
  expect_equal(nrow(ef), 8)
  expect_true(all(ef$p_dynamic >= ef$p_myopic - 1e-8))
})

test_that("Model E rejects unbalanced panels", {
  p <- roy_dynamic_sim(20, T = 3)
  expect_error(roy_dynamic(p[-1, ]), "balanced")
})
