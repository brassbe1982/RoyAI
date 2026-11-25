# Model E Monte Carlo (monograph Chapter 7, design "Model E").
# Rscript mc_modelE.R [R] [cores]
# DGP: roy_dynamic_sim(N = 1000, T = 10), mu_k = 0, sd_k = 0.25, sigma_e = 0.3,
# cost = 0.1 + 0.1 x, lambda = 0.1, beta = 0.95, psi in {0, 0.1}.
# Estimators: CCP (Hotz-Miller) and NFXP. Static comparison on the final-period
# cross-section: naive earnings gap and the Gaussian Roy model (Heckman
# two-step, x as excluded cost shifter) against the true ATT/ATE of period T.
suppressMessages(library(RoyAI))
args <- commandArgs(TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 100
cores <- if (length(args) >= 2) as.integer(args[2]) else 1
out <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "results")
truth0 <- c(mu0 = 1, sigma0 = 0.3, mu_k = 0, sd_k = 0.25, sigma_e = 0.3, psi = NA,
            `cost:(Intercept)` = 0.1, `cost:x` = 0.1, lambda = 0.1)

one <- function(r, psi) {
  set.seed(7e5 + 1e4 * psi + r)
  p <- roy_dynamic_sim(1000, T = 10, psi = psi)
  tr <- truth0; tr["psi"] <- psi
  res <- list()
  for (m in c("ccp", "nfxp")) {
    t0 <- Sys.time()
    f <- tryCatch(suppressWarnings(roy_dynamic(p, cost = ~ x, method = m)), error = function(e) NULL)
    if (is.null(f)) next
    b <- coef(f); se <- sqrt(pmax(diag(vcov(f)), 0))
    res[[m]] <- data.frame(est = m, param = names(b), value = b, se = se, truth = tr[names(b)],
                           secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  fin <- p[p$t == 10, ]
  att <- mean(fin$kappa[fin$s == 1]); ate <- mean(fin$kappa)
  naive <- mean(fin$y[fin$s == 1]) - mean(fin$y[fin$s == 0])
  st <- tryCatch({
    rf <- suppressWarnings(roy(y ~ 1, s ~ x, data = fin, method = "ml"))
    ef <- roy_effects(rf, se = FALSE)
    c(ATE = ef$estimate[grepl("ATE", ef$population)][1], ATT = ef$estimate[grepl("ATT", ef$population)][1])
  }, error = function(e) c(ATE = NA, ATT = NA))
  res$static <- data.frame(est = c("naive_gap", "static_roy", "static_roy"),
                           param = c("ATT_T", "ATT_T", "ATE_T"),
                           value = c(naive, st[["ATT"]], st[["ATE"]]), se = NA,
                           truth = c(att, att, ate), secs = NA)
  o <- do.call(rbind, res); o$rep <- r; o$psi <- psi; rownames(o) <- NULL; o
}

raw <- do.call(rbind, lapply(c(0, 0.1), function(psi) {
  t0 <- Sys.time()
  r <- do.call(rbind, parallel::mclapply(seq_len(R), one, psi = psi, mc.cores = cores))
  message("psi=", psi, ": ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  r
}))
saveRDS(raw, file.path(out, "mc_modelE_raw.rds"))
key <- unique(raw[, c("psi", "est", "param")])
tab <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  k <- key[i, ]
  s <- raw[raw$psi == k$psi & raw$est == k$est & raw$param == k$param & is.finite(raw$value), ]
  data.frame(k, truth = mean(s$truth), mean = mean(s$value), bias = mean(s$value - s$truth),
             rmse = sqrt(mean((s$value - s$truth)^2)),
             coverage = if (all(is.na(s$se))) NA else mean(abs(s$value - s$truth) <= 1.96 * s$se),
             secs = if (all(is.na(s$secs))) NA else mean(s$secs), R = nrow(s))
}))
write.csv(tab, file.path(out, "mc_modelE_summary.csv"), row.names = FALSE)
print(tab, digits = 3)
