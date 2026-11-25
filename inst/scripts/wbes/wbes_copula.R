# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# WBES copula Roy model (Model A): outcome = relative log labour productivity (lp) or wage (lw);
# selection on the extended controls + leave-one-out peer-adoption instrument z.
# Random subsample of firms with a valid instrument (n = NSUB) to keep each fit tractable.
# Usage: Rscript wbes_copula.R <family> <outcome> <delta>
suppressPackageStartupMessages(library(RoyAI))
source(file.path(ROYAI_ROOT, "RoyAI-data/results/wbes/scripts/wbes_prep.R"))
a <- commandArgs(TRUE); fam <- a[1]; y <- a[2]; dl <- as.numeric(a[3]); NSUB <- as.integer(Sys.getenv("NSUB", "20000"))
OUT <- file.path(ROYAI_ROOT, "RoyAI-data/results/wbes/copula"); dir.create(OUT, showWarnings = FALSE)
fn <- file.path(OUT, sprintf("%s_%s_d%.2f.rds", fam, y, dl)); if (file.exists(fn)) quit()
d <- wbes_prep(); d <- d[!is.na(z)]; set.seed(20261125); d <- as.data.frame(d[sample.int(nrow(d), min(NSUB, nrow(d)))])
d$ya <- d[[y]] - dl * d$z
cc <- paste(CTRL_W, EXT_W)
t0 <- Sys.time()
fit <- tryCatch(roy_copula_fit(as.formula(paste("ya ~", cc)), as.formula(paste("D ~", cc, "+ z")), d, family = fam), error = function(e) e)
r <- if (inherits(fit, "error")) list(error = conditionMessage(fit)) else
  list(te = roy_te(fit), coef = coef(fit)[grepl("^(tau|sigma|nu)", names(coef(fit)))], se = sqrt(diag(vcov(fit)))[grepl("^(tau|sigma|nu)", names(coef(fit)))],
       loglik = fit$loglik, aic = AIC(fit), conv = fit$convergence, multistart = fit$multistart, n = nrow(d))
r$family <- fam; r$y <- y; r$delta <- dl; r$minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(r, fn); print(r[c("te","coef","aic","minutes")])
