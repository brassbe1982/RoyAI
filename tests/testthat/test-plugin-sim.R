test_that("a new dependence family plugs into ML, simulation and effects", {
  # Clayton copula (exchangeable), theta > 0
  clayton <- roy_dependence(
    name = "clayton_test", npar = 1L, par_names = "theta",
    h = function(v, u, theta) (u^(-theta) * (u^(-theta) + v^(-theta) - 1)^(-1 - 1 / theta)) * u^(-1) * u,
    linkinv = function(eta) exp(eta), link = function(theta) log(theta), start = 1,
    tau = function(theta) theta / (theta + 2))
  # h(v|u) = u^(-theta-1) (u^-theta + v^-theta - 1)^(-1-1/theta)
  clayton$h <- function(v, u, theta) u^(-theta - 1) * (u^(-theta) + v^(-theta) - 1)^(-1 - 1 / theta)
  roy_register_dependence(clayton, overwrite = TRUE)
  expect_true("clayton_test" %in% roy_dependence_families())
  dat <- roy_sim_A(2500, dependence = "clayton_test", theta0 = 1, theta1 = 2, seed = 31)
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, dependence = "clayton_test")
  expect_true(f$converged)
  expect_equal(unname(f$par[c("theta0", "theta1")]), c(1, 2), tolerance = 0.5)
  e <- roy_effects(f)
  expect_true(all(is.finite(e$estimate)))
  fg <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
  expect_equal(nrow(roy_ic(clayton = f, gaussian = fg)), 2)
})

test_that("registered estimation methods are dispatched", {
  roy_register_method("A", "dummy", function(spec, ...) {
    f <- RoyAI:::.fit_A_twostep(spec); f$method <- "dummy"; f
  }, overwrite = TRUE)
  dat <- roy_sim_A(800, seed = 1)
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, method = "dummy")
  expect_equal(f$method, "dummy")
})

test_that("roy_mc summarizes bias and coverage", {
  mc <- roy_mc(function() roy_sim_A(800),
               function(d) roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = d, method = "twostep"),
               R = 8, seed = 1)
  expect_s3_class(mc, "roy_mc")
  expect_true(all(c("bias", "rmse", "coverage") %in% names(mc$summary)))
})

test_that("bootstrap standard errors", {
  dat <- roy_sim_A(600, seed = 4)
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat, method = "twostep", se = "bootstrap", B = 15)
  expect_match(f$se_type, "bootstrap")
  expect_true(all(is.finite(roy_effects(f)$se)))
})

test_that("PIAAC preparation pipeline runs on synthetic data", {
  ex <- piaac_example(2000)
  d <- piaac_onet_prep(ex$piaac, ex$exposure)
  expect_true(all(c("log_wage", "d", "s", "ai_exposure") %in% names(d)))
  expect_setequal(unique(d$s), 0:2)
  expect_true(length(attr(d, "notes")) >= 3)
  # robustness exposure measure and alternative variable names
  d2 <- piaac_onet_prep(ex$piaac, ex$exposure, exposure_col = "gpt_beta")
  expect_false(isTRUE(all.equal(d$ai_exposure, d2$ai_exposure)))
  p2 <- ex$piaac; names(p2)[names(p2) == "EARNHRBONUSPPPC2"] <- "EARNHRBONUSPPP"
  d3 <- piaac_onet_prep(p2, ex$exposure)
  expect_equal(unname(attr(d3, "variables")["earnings"]), "EARNHRBONUSPPP")
  # PUF missing codes become NA; AGEG10LFS replaces a suppressed AGE_R
  p4 <- ex$piaac; p4$EARNHRBONUSPPPC2[1:10] <- 999999999996; p4$ISCO2C[11:20] <- 999
  p4$AGEG10LFS <- 3; p4$AGE_R <- NULL
  d4 <- piaac_onet_prep(p4, ex$exposure)
  expect_false(any(d4$id %in% p4$SEQID[1:20]))
  expect_true(all(d4$age == 39.5))
  p2$EARNHRBONUSPPP <- NULL
  expect_error(piaac_onet_prep(p2, ex$exposure), "earnings")
  xw <- data.frame(soc = c("15-1252", "43-4051"), isco08 = c("2512", "4222"))
  exs <- data.frame(soc = c("15-1252.00", "43-4051.00"), aioe = c(1.2, 0.8))
  i2 <- onet_to_isco2(exs, xw, "aioe")
  expect_equal(i2$isco08_2d, c("25", "42"))
  ci <- isco_soc_crosswalk()
  expect_true(all(c("isco08", "soc") %in% names(ci)))
  expect_gt(length(unique(substr(ci$isco08, 1, 2))), 40)
  expect_equal(nrow(piaac_analysis_countries()), 6)
})

test_that("O*NET complementarity and SOC exposure assembly", {
  items <- c("Face-to-Face Discussions with Individuals and Within Teams", "Public Speaking",
             "Health and Safety of Other Workers", "Work Outcomes and Results of Other Workers",
             "Outdoors, Exposed to All Weather Conditions", "Physical Proximity",
             "Consequence of Error", "Freedom to Make Decisions", "Frequency of Decision Making",
             "Determine Tasks, Priorities and Goals", "Degree of Automation",
             "Importance of Repeating Same Tasks")
  wc <- expand.grid(code = c("11-1011.00", "43-9041.00"), item = items, stringsAsFactors = FALSE)
  wc <- data.frame(`O*NET-SOC Code` = wc$code, `Element Name` = wc$item, `Scale ID` = "CX",
                   `Data Value` = ifelse(wc$code == "11-1011.00", 5, 1), check.names = FALSE)
  jz <- data.frame(`O*NET-SOC Code` = c("11-1011.00", "43-9041.00"), `Job Zone` = c(5, 2), check.names = FALSE)
  cm <- onet_complementarity(wc, jz)
  expect_gt(cm$complementarity[cm$soc == "11-1011"], cm$complementarity[cm$soc == "43-9041"])
  aioe <- data.frame(`SOC Code` = c("11-1011", "43-9041"), AIOE = c(1.3, 1.2), check.names = FALSE)
  gpt <- data.frame(`O*NET-SOC Code` = c("11-1011.00", "11-1011.03", "43-9041.00"),
                    dv_rating_beta = c(0.4, 0.6, 0.8), human_rating_beta = c(0.3, 0.3, 0.7),
                    check.names = FALSE)
  ex <- onet_exposure_soc(aioe, gpt, cm)
  expect_equal(ex$gpt_beta[ex$soc == "11-1011"], 0.5)
  expect_equal(nrow(ex), 2)
})

test_that("copula ML fits use multistart and accept dependence starts", {
  skip_if_not("frank" %in% roy_dependence_families())
  set.seed(11)
  d <- roy_sim_A(n = 600, rho0 = -0.5, rho1 = -0.4)
  f <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = d, dependence = "frank")
  expect_equal(f$multistart$start, c("default", "twostep_tau", "flipped_tau"))
  expect_equal(sum(f$multistart$chosen), 1)
  expect_equal(f$loglik, max(f$multistart$loglik), tolerance = 1e-3)
  f1 <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = d, dependence = "frank",
            start = list(-0.3, -0.3), multistart = FALSE)
  expect_equal(f1$multistart$start, c("user", "default"))
})
