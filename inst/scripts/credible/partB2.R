source("./prep2.R")
res <- readRDS("partB.rds")
d <- prep_main(1)
med <- ave(ifelse(d$d == 1, d$comp_ind, NA), d$country, FUN = function(x) median(x, na.rm = TRUE))
d$s <- ifelse(d$d == 0, 0L, ifelse(d$comp_ind >= med, 1L, 2L))
dz <- d[!is.na(d$z_indict), ]
JOBL <- "+ firmsize_f + sector_f + contract_f + isic1_f + tenure"
out <- list()
for (cc in c("base","joblite")) for (dl in c(0, 0.0180)) {
  ctrl <- if (cc == "base") CTRL else paste(CTRL, JOBL)
  dz$y_adj <- dz$log_wage - dl * dz$ict_home
  fb <- roy(as.formula(paste("y_adj ~", ctrl)), as.formula(paste("s ~", ctrl, "+ ict_home + z_indict + z_indmeth")), dz, correction = "dmf")
  e <- roy_effects(fb); e$controls <- cc; e$delta <- dl
  # collinearity of correction terms within each regime
  cf <- coef(fb); mt <- grep(":m[0-9]$", names(cf), value = TRUE)
  e$mterms_absmax <- max(abs(cf[mt]))
  out[[paste(cc, dl)]] <- e
  # OLS benchmark
  m <- lm(as.formula(paste("log_wage ~ factor(s) +", ctrl)), dz)
  print(c(cc, dl)); print(coeftest(m, vcov = vcovCL(m, cluster = dz$ind))[c("factor(s)1","factor(s)2"),1:2])
  print(e[, setdiff(names(e), "mterms_absmax")], digits = 3)
}
res$B_effects <- out
saveRDS(res, "partB.rds")
