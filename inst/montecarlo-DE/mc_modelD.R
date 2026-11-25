# Model D Monte Carlo (monograph Chapter 7, design "Model D").
# Rscript mc_modelD.R [R] [cores]
# DGP: roy_continuous_sim() defaults: mu = 1 + 0.5 x1 + 0.8 a - 0.4 a^2,
# log sigma(a) = log(0.6) + 0.5 a, gamma = (0.6, 0.3, 0.8, 0.6) (~30% non-users),
# tau = 0.4, families Gaussian and Frank, n in {1000, 3000}.
# Estimators: copula ML (correct family), Gaussian ML on Frank data
# (misspecified), control-function two-step (Gaussian, homoskedastic), naive OLS.
suppressMessages(library(RoyAI))
args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 100
cores <- if (length(args) >= 2) as.integer(args[2]) else 1
out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "results")
set.seed(20261125)
big <- data.frame(x1 = rnorm(2e5), z1 = rnorm(2e5), z2 = rbinom(2e5, 1, 0.5))
A <- c(0.2, 0.5, 0.8)
truth <- lapply(c(gaussian = "gaussian", frank = "frank"), function(f)
  roy_continuous_truth(tau = 0.4, family = f, covariates = big, a = A, grid = c(0, 0.5, 0.95)))

one <- function(r, n, dgp) {
  set.seed(1e5 * match(dgp, c("gaussian", "frank")) + 10 * n + r)
  d <- roy_continuous_sim(n, tau = 0.4, family = dgp)
  res <- list()
  fits <- list(ml = dgp)
  if (dgp != "gaussian") fits$ml_gauss <- "gaussian"
  for (k in names(fits)) {
    fit <- tryCatch(suppressWarnings(roy_continuous(y ~ x1, a ~ x1 + z1 + z2, d, family = fits[[k]])),
                    error = function(e) NULL)
    if (is.null(fit)) next
    b <- coef(fit); se <- sqrt(diag(vcov(fit)))
    e <- roy_continuous_effects(fit, a = A, grid = c(0, 0.5, 0.95))
    res[[k]] <- data.frame(est = k, param = c("beta_a", "beta_a2", "alpha_a", "tau",
                                              paste(e$estimand, e$a)),
                           value = c(b[c("beta:a", "beta:I(a^2)", "alpha:a", "tau")], e$estimate),
                           se = c(se[c("beta:a", "beta:I(a^2)", "alpha:a", "tau")], e$se),
                           regular = fit$regular)
  }
  ols <- lm(y ~ x1 + a + I(a^2), d)
  cf <- RoyAI:::.rc_twostep(list(y = d$y, a = d$a, tA = ifelse(d$a > 0, -log1p(-d$a), 0),
                                 X = cbind(1, d$x1), D = cbind(d$a, d$a^2),
                                 Z = cbind(1, d$x1, d$z1, d$z2)), NULL)
  for (k in c("ols", "cf")) {
    bb <- if (k == "ols") coef(ols)[3:4] else cf$beta[3:4]
    ate <- bb[1] * 0.1 + bb[2] * ((A + 0.1)^2 - A^2)
    res[[k]] <- data.frame(est = k, param = c("beta_a", "beta_a2", paste("ATE", A)),
                           value = c(bb, ate), se = NA, regular = NA)
  }
  out <- do.call(rbind, res); out$rep <- r; out$n <- n; out$dgp <- dgp
  rownames(out) <- NULL; out
}

cells <- expand.grid(n = c(1000, 3000), dgp = c("gaussian", "frank"), stringsAsFactors = FALSE)
raw <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  t0 <- Sys.time()
  r <- do.call(rbind, parallel::mclapply(seq_len(R), one, n = cells$n[i], dgp = cells$dgp[i],
                                         mc.cores = cores))
  message(cells$dgp[i], " n=", cells$n[i], ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  r
}))
saveRDS(raw, file.path(out, "mc_modelD_raw.rds"))

tv <- function(dgp, param) {
  tr <- truth[[dgp]]; key <- paste(tr$estimand, tr$a)
  base <- c(beta_a = 0.8, beta_a2 = -0.4, alpha_a = 0.5, tau = 0.4)
  ifelse(param %in% names(base), base[param], tr$truth[match(param, key)])
}
raw$truth <- tv(raw$dgp[1], raw$param)
for (g in unique(raw$dgp)) raw$truth[raw$dgp == g] <- tv(g, raw$param[raw$dgp == g])
key <- unique(raw[, c("dgp", "n", "est", "param")])
tab <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  k <- key[i, ]
  s <- raw[raw$dgp == k$dgp & raw$n == k$n & raw$est == k$est & raw$param == k$param, ]
  data.frame(k, truth = s$truth[1], mean = mean(s$value), bias = mean(s$value - s$truth),
             rmse = sqrt(mean((s$value - s$truth)^2)),
             coverage = if (all(is.na(s$se))) NA else mean(abs(s$value - s$truth) <= 1.96 * s$se),
             share_regular = if (all(is.na(s$regular))) NA else mean(s$regular), R = nrow(s))
}))
write.csv(tab, file.path(out, "mc_modelD_summary.csv"), row.names = FALSE)
print(tab, digits = 3)
