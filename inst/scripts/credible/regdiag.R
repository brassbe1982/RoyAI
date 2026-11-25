# Why are the copula Roy maxima flagged irregular? Hessian eigenvalues and score at the optimum.
source("./prep2.R")
d <- prep_main(1)
out <- list()
for (dl in c(0, 0.0180)) {
  dd <- d; dd$y_adj <- dd$log_wage - dl * dd$ict_home
  fit <- roy_copula_fit(as.formula(paste("y_adj ~", CTRL)), as.formula(paste("d ~", CTRL, "+ ict_home")), dd, family = "frank")
  H <- solve(fit$vcov_raw); ev <- eigen(H, symmetric = TRUE)
  nll <- function(p) -fit$loglik  # placeholder
  out[[as.character(dl)]] <- list(min_eig = min(ev$values), n_neg = sum(ev$values <= 0),
     worst_dir = head(sort(abs(ev$vectors[, which.min(ev$values)]), decreasing = TRUE), 6),
     worst_names = names(fit$par)[order(-abs(ev$vectors[, which.min(ev$values)]))][1:6],
     cond = max(ev$values) / max(min(ev$values), 1e-12), ms = fit$multistart, te = roy_te(fit), loglik = fit$loglik)
  print(out[[as.character(dl)]])
}
saveRDS(out, "regdiag.rds")
