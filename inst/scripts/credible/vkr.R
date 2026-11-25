# van Kippersluis-Rietveld plug-in: LATE(delta_hat) with SE combining the two independent samples
source("./prep2.R")
d <- prep_main(1); dz <- d[!is.na(d$z_indict), ]
dn <- d0; dn$skill_lit <- dn$lit1; dn <- droplevels(dn[dn$exposed == 0 & !is.na(dn$z_indict) & !is.na(dn$ict_home), ])
out <- list()
for (cc in c("base","job")) for (pl in c("cond","uncond")) {
  ctrl <- if (cc == "base") CTRL else paste(CTRL, JOB)
  rhs <- if (pl == "cond") "ict_work + ict_home +" else "ict_home +"
  mp <- lm(as.formula(paste("log_wage ~", rhs, ctrl)), dn)
  dh <- coef(mp)["ict_home"]; dse <- sqrt(vcovCL(mp, cluster = dn$ind)["ict_home","ict_home"])
  dz$y_adj <- dz$log_wage - dh * dz$ict_home
  m <- ivreg(as.formula(paste("y_adj ~ d +", ctrl, "| ict_home +", ctrl)), data = dz)
  fs <- lm(as.formula(paste("d ~", ctrl, "+ ict_home")), dz); pi <- coef(fs)["ict_home"]
  se <- sqrt(vcovCL(m, cluster = dz$ind)["d","d"] + (dse / pi)^2)
  out[[paste(cc, pl)]] <- data.frame(controls = cc, placebo = pl, delta = dh, delta_se = dse, pi = pi,
                                     late = coef(m)["d"], se = se, lo = coef(m)["d"] - 1.96 * se, hi = coef(m)["d"] + 1.96 * se)
}
out <- do.call(rbind, out); rownames(out) <- NULL; print(out, digits = 3); saveRDS(out, "vkr.rds")
