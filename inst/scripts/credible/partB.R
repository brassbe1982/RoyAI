source("./prep2.R")
d <- prep_main(1)
# worker-level arrangement among adopters: split by individual complementarity (within country, among adopters)
med <- ave(ifelse(d$d == 1, d$comp_ind, NA), d$country, FUN = function(x) median(x, na.rm = TRUE))
d$s <- ifelse(d$d == 0, 0L, ifelse(d$comp_ind >= med, 1L, 2L))
d$s_occ <- ifelse(d$d == 0, 0L, ifelse(d$complementarity > median(d$complementarity), 1L, 2L))
res <- list()
res$shares <- table(d$s)
# Prop 5.3 check: share of adopter variation in s that is within occupation
ad <- d[d$d == 1, ]
res$within_occ_var <- c(worker_level = 1 - summary(lm((s == 2) ~ occ2, ad))$r.squared,
                        occ_level = 1 - summary(lm((s_occ == 2) ~ occ2, ad))$r.squared)
res$occ_both <- sum(apply(table(ad$occ2, ad$s), 1, min) >= 10); res$n_occ <- nlevels(ad$occ2)
print(res[c("shares","within_occ_var","occ_both","n_occ")])
dz <- d[!is.na(d$z_indict), ]
cl <- dz$ind
# descriptive premia with occupation effects (OLS)
m <- lm(as.formula(paste("log_wage ~ factor(s) +", CTRL)), dz)
res$ols <- coeftest(m, vcov = vcovCL(m, cluster = cl))[c("factor(s)1","factor(s)2"), 1:2]; print(res$ols)
# multinomial first stage: instrument relevance for each margin (LR tests)
library(nnet)
fs_full <- multinom(as.formula(paste("factor(s) ~", CTRL, "+ ict_home + z_indict + z_indmeth")), dz, trace = FALSE, maxit = 500)
fs_noz <- multinom(as.formula(paste("factor(s) ~", CTRL)), dz, trace = FALSE, maxit = 500)
res$fs_LR <- c(LR = 2 * (logLik(fs_full) - logLik(fs_noz)))
res$fs_coef <- coef(fs_full)[, c("ict_home","z_indict","z_indmeth")]
res$fs_se <- summary(fs_full)$standard.errors[, c("ict_home","z_indict","z_indmeth")]
print(res$fs_coef); print(res$fs_se)
# relevance for augmented-vs-dominant margin among adopters
lp <- glm(as.formula(paste("(s==1) ~", CTRL, "+ ict_home + z_indict + z_indmeth")), dz[dz$d==1,], family = binomial("probit"))
res$aug_margin <- coeftest(lp)[c("ict_home","z_indict","z_indmeth"), ]; print(res$aug_margin)
# Model B two-step DMF via the package, Z = ict_home + z_indict + z_indmeth
for (dl in c(0, 0.0180)) {
  dz$y_adj <- dz$log_wage - dl * dz$ict_home
  fb <- tryCatch(roy(as.formula(paste("y_adj ~", CTRL)), as.formula(paste("s ~", CTRL, "+ ict_home + z_indict + z_indmeth")),
                     dz, correction = "dmf"), error = function(e) e)
  res[[paste0("B_dmf_d", dl)]] <- fb
  if (!inherits(fb, "error")) { print(dl); print(summary(fb)); print(tryCatch(roy_te(fb), error=function(e) conditionMessage(e))) } else print(fb)
}
saveRDS(res, "partB.rds")
