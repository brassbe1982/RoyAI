# Check for multiple local optima: refit the delta = 0 model starting from
# the delta = 0.03 solution (and vice versa), and compare log-likelihoods.
suppressPackageStartupMessages(library(RoyAI))
src <- readLines(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "run_piaac_labour.R"))
end <- grep("^# ---- family selection", src); eval(parse(text = src[1:(end - 1)]))
sel <- readRDS(file.path(out_dir, "sel_pv1.rds"))
adj <- function(dl) { dd <- d1; dd$log_wage <- dd$log_wage - dl * dd$ict_home; dd }
out <- cache("multimodal", {
  fr03 <- roy_copula_fit(f_out(), f_sel(), adj(0.03), family = "frank")
  # same parameters except the intercept/ict_home shift do not apply: start from fr03 raw par
  fr0_from03 <- roy_copula_fit(f_out(), f_sel(), adj(0), family = "frank", start = fr03$par)
  fr03_from0 <- roy_copula_fit(f_out(), f_sel(), adj(0.03), family = "frank", start = sel$fits$frank$par)
  list(fr03 = fr03, fr0_from03 = fr0_from03, fr03_from0 = fr03_from0)
})
show <- function(nm, f) { te <- roy_te(f, se = FALSE); cf <- coef(f)
  cat(sprintf("%-22s logLik %.2f tau0 %.2f tau1 %.2f ATE %.3f ATT %.3f ATU %.3f\n", nm, f$loglik, cf[["tau0"]], cf[["tau1"]],
              te$estimate[1], te$estimate[2], te$estimate[3])) }
show("delta0 default start", sel$fits$frank); show("delta0 from delta.03", out$fr0_from03)
show("delta.03 default", out$fr03); show("delta.03 from delta0", out$fr03_from0)
