# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# Copula Roy (Model A) on PIAAC with the exclusion restriction relaxed:
# outcome = log wage - delta * ict_home, delta calibrated on non-exposed occupations
# (zero-first-stage group, van Kippersluis & Rietveld 2018). PV1.
source("./prep2.R")
out <- file.path(ROYAI_ROOT, "RoyAI-data/results/credible/cells")
d <- prep_main(1)
deltas <- list(base = c(0, 0.0180, 0.0301), job = c(0, 0.0171, 0.0234))
cells <- expand.grid(family = c("frank", "gaussian"), controls = c("base", "job"), k = 1:3, stringsAsFactors = FALSE)
cells <- cells[order(cells$controls == "job", cells$family == "gaussian"), ]
run <- function(i) {
  cc <- cells[i, ]; dl <- deltas[[cc$controls]][cc$k]
  f <- file.path(out, sprintf("%s_%s_d%.4f.rds", cc$family, cc$controls, dl))
  if (file.exists(f)) return(f)
  ctrl <- if (cc$controls == "base") CTRL else paste(CTRL, JOB)
  dd <- d; dd$y_adj <- dd$log_wage - dl * dd$ict_home
  t0 <- Sys.time()
  fit <- tryCatch(roy_copula_fit(as.formula(paste("y_adj ~", ctrl)), as.formula(paste("d ~", ctrl, "+ ict_home")),
                                 dd, family = cc$family), error = function(e) e)
  r <- if (inherits(fit, "error")) list(error = conditionMessage(fit)) else
    list(te = roy_te(fit), coef = coef(fit)[grepl("tau|sigma", names(coef(fit)))], loglik = fit$loglik,
         aic = AIC(fit), conv = fit$convergence, multistart = fit$multistart, n = nrow(dd))
  r$family <- cc$family; r$controls <- cc$controls; r$delta <- dl; r$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(r, f); f
}
invisible(parallel::mclapply(seq_len(nrow(cells)), run, mc.cores = 4, mc.preschedule = FALSE))
cat("done\n")
