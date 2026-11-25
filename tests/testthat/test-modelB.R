dat <- roy_sim_B(8000, seed = 21)
tr <- attr(dat, "truth")

test_that("Model B two-step recovers outcome equations", {
  f <- roy(y ~ x1 + x2, s ~ x1 + x2 + z1 + z2, data = dat)
  expect_s3_class(f, "roy_fit_B")
  b <- grep("^y[0-2]:(\\(Intercept\\)|x1|x2)$", names(tr), value = TRUE)
  se <- sqrt(diag(f$vcov))[b]
  expect_true(all(abs(f$par[b] - tr[b]) < 4 * se))
})

test_that("Model B effects and diagnostics", {
  f <- roy(y ~ x1 + x2, s ~ x1 + x2 + z1 + z2, data = dat)
  e <- roy_effects(f)
  expect_equal(nrow(e), 12)
  Y <- as.matrix(dat[, c("y0", "y1", "y2")])
  att <- mean((Y[, 2] - Y[, 1])[dat$s == 1])
  est <- e$estimate[e$contrast == "ai_augmented vs human_only" & e$population == "S = ai_augmented"]
  sev <- e$se[e$contrast == "ai_augmented vs human_only" & e$population == "S = ai_augmented"]
  expect_lt(abs(est - att), 4 * sev)
  dg <- roy_diagnostics(f)
  expect_equal(nrow(dg$overlap$summary), 3)
})

test_that("Lee correction fits and gives NA counterfactuals", {
  f <- roy(y ~ x1 + x2, s ~ x1 + x2 + z1 + z2, data = dat, correction = "lee")
  expect_true(all(is.na(roy_potential(f, 1, 0))))
  expect_true(all(is.finite(roy_potential(f, 1, 1))))
})

test_that("ML is refused for Model B in the core", {
  expect_error(roy(y ~ x1 + x2, s ~ x1 + x2 + z1 + z2, data = dat, method = "ml"), "twostep")
})
