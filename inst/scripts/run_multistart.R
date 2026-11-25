# Multistart check of the copula family comparison (main specification, PV1).
# Each family is refitted from the optima of the other families and from the
# positive-dependence mode, and the best log-likelihood is kept.
suppressPackageStartupMessages(library(RoyAI))
src <- readLines(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "run_piaac_labour.R"))
end <- grep("^# ---- family selection", src); eval(parse(text = src[1:(end - 1)]))
sel <- readRDS(file.path(out_dir, "sel_pv1.rds")); mm <- readRDS(file.path(out_dir, "multimodal.rds"))
base <- function(p) p[!grepl("log_nu", names(p))]
starts <- list(frank = base(sel$fits$frank$par), gumbel = base(sel$fits$gumbel$par), t = base(sel$fits$t$par),
               posmode = base(mm$fr0_from03$par))
fams <- c("t", "gumbel", "frank", "clayton", "gaussian")
jobs <- expand.grid(family = fams, start = names(starts), nu0 = NA, stringsAsFactors = FALSE)
jobs <- jobs[jobs$family != jobs$start, ]
jobs <- rbind(jobs, data.frame(family = "t", start = c("frank", "gumbel"), nu0 = base::log(0.05)))
fits <- cache("multistart", pmap(seq_len(nrow(jobs)), function(k) {
  fam <- jobs$family[k]; st <- starts[[jobs$start[k]]]
  if (fam == "t") st <- c(st, log_nu_m2 = if (is.na(jobs$nu0[k])) base::log(3) else jobs$nu0[k])
  f <- tryCatch(roy_copula_fit(f_out(), f_sel(), d1, family = fam, start = unname(st)), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  te <- roy_te(f); cf <- coef(f)
  list(family = fam, start = jobs$start[k], nu0 = jobs$nu0[k], loglik = f$loglik, AIC = AIC(f), conv = f$convergence,
       tau0 = cf[["tau0"]], tau1 = cf[["tau1"]], nu = if (fam == "t") cf[["nu"]] else NA,
       te = te, par = f$par)
}))
rows <- do.call(rbind, lapply(fits, function(x) if (!is.null(x)) data.frame(family = x$family, start = x$start, nu0 = x$nu0,
  loglik = x$loglik, AIC = x$AIC, conv = x$conv, tau0 = x$tau0, tau1 = x$tau1, nu = x$nu,
  ATE = x$te$estimate[1], ATT = x$te$estimate[2], ATU = x$te$estimate[3])))
# add the default-start fits
for (f in fams) { x <- sel$fits[[f]]; te <- res$sel_te[[f]] <- roy_te(x); cf <- coef(x)
  rows <- rbind(rows, data.frame(family = f, start = "default", nu0 = NA, loglik = x$loglik, AIC = AIC(x), conv = x$convergence,
    tau0 = cf[["tau0"]], tau1 = cf[["tau1"]], nu = if (f == "t") cf[["nu"]] else NA,
    ATE = te$estimate[1], ATT = te$estimate[2], ATU = te$estimate[3])) }
rows <- rows[order(rows$family, rows$AIC), ]
write.csv(rows, file.path(out_dir, "multistart_all.csv"), row.names = FALSE)
best <- do.call(rbind, lapply(split(rows, rows$family), function(r) r[1, ]))
best <- best[order(best$AIC), ]; best$dAIC <- best$AIC - min(best$AIC)
best$weight <- exp(-best$dAIC / 2) / sum(exp(-best$dAIC / 2))
write.csv(best, file.path(out_dir, "multistart_best.csv"), row.names = FALSE)
print(rows, digits = 4); cat("\nBest per family:\n"); print(best, digits = 4)
