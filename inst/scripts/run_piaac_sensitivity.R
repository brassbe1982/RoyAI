# Sensitivity of the copula Roy estimates to a direct effect of ICT use at
# home on earnings (a "plausibly exogenous" analysis in the spirit of
# Conley, Hansen and Rossi, 2012). For a direct effect delta, the outcome is
# log_wage - delta * ict_home, and the model is re-estimated with ict_home
# still excluded. delta = 0 is the main specification.
# Usage: Rscript run_piaac_sensitivity.R <puf_dir> <exposure.csv> <out_dir> [cores]
suppressPackageStartupMessages(library(RoyAI))
src <- readLines(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "run_piaac_labour.R"))
end <- grep("^# ---- family selection", src)
eval(parse(text = src[1:(end - 1)]))   # data, formulas, helpers, d1

deltas <- c(0, 0.01, 0.02, 0.03, 0.04)
fams <- c("frank", "gaussian", "clayton")
grid <- expand.grid(delta = deltas, family = fams, stringsAsFactors = FALSE)
out <- cache("sens_delta", pmap(seq_len(nrow(grid)), function(k) {
  dd <- d1; dd$log_wage <- dd$log_wage - grid$delta[k] * dd$ict_home
  fit <- tryCatch(roy_copula_fit(f_out(), f_sel(), dd, family = grid$family[k]), error = function(e) e)
  if (inherits(fit, "error")) return(NULL)
  te <- roy_te(fit); cf <- coef(fit)
  data.frame(delta = grid$delta[k], family = grid$family[k], AIC = AIC(fit),
             tau0 = cf[["tau0"]], tau1 = cf[["tau1"]],
             ATE = te$estimate[1], ATT = te$estimate[2], ATU = te$estimate[3],
             se_ATE = te$se[1], se_ATT = te$se[2], se_ATU = te$se[3])
}))
tab <- do.call(rbind, out)
# OLS benchmark at each delta
tab$OLS <- sapply(tab$delta, function(dl) {
  dd <- d1; dd$log_wage <- dd$log_wage - dl * dd$ict_home
  coef(lm(update(f_out(), . ~ . + d), data = dd))[["d"]]
})
saveRDS(tab, file.path(out_dir, "sensitivity_delta.rds"))
write.csv(tab, file.path(out_dir, "sensitivity_delta.csv"), row.names = FALSE)
print(tab, digits = 3)
