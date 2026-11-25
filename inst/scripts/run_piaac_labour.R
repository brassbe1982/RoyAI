# RoyAI, Labour Economics application: copula generalized Roy model on
# PIAAC Cycle 2 (2023) public-use files linked to occupational AI exposure.
#
# Usage: Rscript run_piaac_labour.R <puf_dir> <exposure_isco08_2d.csv> <out_dir> [cores]
#
# Main design (binary, J = 2): workers in AI-exposed occupations (AIOE above
# the pooled median). Treatment d = 1 if ICT use at work is in the top
# within-country tercile ("AI-ready work"). Outcome: log hourly earnings incl.
# bonuses (PPP USD). Controls: literacy (PV), schooling, age, age^2, female,
# country and 2-digit occupation effects. Excluded from the outcome equation:
# ICT use at home (ICTHOMEC2), an assumption that the paper discusses.
# Secondary: ordered 3-regime model on all occupations (human-only <
# AI-augmented < AI-dominant, from the exposure x complementarity split).
suppressPackageStartupMessages(library(RoyAI))
args <- commandArgs(trailingOnly = TRUE)
puf_dir <- args[1]; exposure_csv <- args[2]; out_dir <- args[3]
cores <- if (length(args) >= 4) as.integer(args[4]) else 1L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console() }
pmap <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
cache <- function(name, expr) {
  f <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  v <- expr; saveRDS(v, f); v
}

keep <- c("SEQID", "CNTRYID", "SPFWT0", "EARNHRBONUSPPPC2", "ISCO2C", "AGE_R", "AGEG10LFS",
          "GENDER_R", "YRSQUALC2", "ICTWORKC2", "ICTHOMEC2",
          paste0("PVLIT", 1:10), paste0("PVNUM", 1:10), paste0("PVAPS", 1:10))
puf <- read_piaac_puf(puf_dir, countries = piaac_analysis_countries()$iso3, columns = keep)
ex <- read.csv(exposure_csv, colClasses = c(isco08_2d = "character"))
log("PUF rows:", nrow(puf))

# ---- data ---------------------------------------------------------------
prep_all <- function(pv, exposure_col = "aioe", ict_q = 2 / 3) {
  d <- piaac_onet_prep(puf, ex, pv = pv, exposure_col = exposure_col,
                       thresholds = list(ict = ict_q), extra = c(ict_home = "ICTHOMEC2"))
  d <- d[stats::complete.cases(d[, c("log_wage", "skill_lit", "skill_num", "age", "schooling",
                                     "ict_work", "ict_home", "ai_exposure")]), ]
  d$country <- factor(d$country)
  d$occ2 <- factor(d$isco08_2d)
  d$exposed <- as.integer(d$ai_exposure > stats::quantile(d$ai_exposure, 0.5))
  d$ict_top <- as.integer(d$ict_work > stats::ave(d$ict_work, d$country, FUN = function(x)
    stats::quantile(x, ict_q)))
  d
}
# exposed-occupation sample; drop 2-digit groups without both regimes
prep_main <- function(pv, ..., min_cell = 20) {
  d <- prep_all(pv, ...)
  d <- d[d$exposed == 1, ]
  d$d <- d$ict_top
  tab <- table(d$occ2, d$d)
  ok <- rownames(tab)[apply(tab, 1, min) >= min_cell]
  droplevels(d[d$occ2 %in% ok, ])
}

ctrl <- "skill_lit + schooling + age + age2 + female + country"
f_out <- function(extra = "+ occ2", skill = "skill_lit")
  as.formula(paste("log_wage ~", sub("skill_lit", skill, ctrl), extra))
f_sel <- function(extra = "+ occ2", skill = "skill_lit", lhs = "d")
  as.formula(paste(lhs, "~", sub("skill_lit", skill, ctrl), extra, "+ ict_home"))

rubin <- function(tabs) {
  M <- length(tabs); est <- sapply(tabs, `[[`, "estimate"); se <- sapply(tabs, `[[`, "se")
  if (is.null(dim(est))) { est <- matrix(est, nrow = 1); se <- matrix(se, nrow = 1) }
  qbar <- rowMeans(est); W <- rowMeans(se^2); B <- apply(est, 1, stats::var)
  data.frame(parameter = tabs[[1]]$parameter, estimate = qbar, se = sqrt(W + (1 + 1 / M) * B),
             between_share = (1 + 1 / M) * B / (W + (1 + 1 / M) * B))
}
tau_of <- function(fit) { cf <- coef(fit); cf[grepl("tau|nu", names(cf))] }

res <- list()
d1 <- prep_main(1)
res$sample <- list(n = nrow(d1), occ2 = levels(d1$occ2), by_country = table(d1$country, d1$d),
                   share_d = mean(d1$d), notes = attr(d1, "notes"))
log("main sample:", nrow(d1), "occupations:", nlevels(d1$occ2), "share d:", round(mean(d1$d), 3))

# descriptives by regime (PV1)
vars <- c("log_wage", "skill_lit", "skill_num", "schooling", "age", "female", "ict_work", "ict_home", "ai_exposure")
res$desc <- data.frame(variable = vars,
                       mean_d0 = sapply(vars, function(v) mean(d1[[v]][d1$d == 0])),
                       mean_d1 = sapply(vars, function(v) mean(d1[[v]][d1$d == 1])),
                       sd_all = sapply(vars, function(v) sd(d1[[v]])))

# OLS benchmark and first stage strength
ols <- lm(update(f_out(), . ~ . + d), data = d1)
res$ols <- summary(ols)$coefficients["d", ]
fs <- glm(f_sel(), data = d1, family = binomial(link = "probit"))
fs0 <- glm(update(f_sel(), . ~ . - ict_home), data = d1, family = binomial(link = "probit"))
res$first_stage <- list(coef = summary(fs)$coefficients["ict_home", ],
                        lr = 2 * (logLik(fs) - logLik(fs0)))
# pseudo outcome check: ict_home in the outcome equation of each regime (OLS)
res$excl_check <- sapply(0:1, function(k)
  summary(lm(update(f_out(), . ~ . + ict_home), data = d1[d1$d == k, ]))$coefficients["ict_home", ])

# ---- family selection (PV1) ----------------------------------------------
log("family selection, PV1")
sel <- cache("sel_pv1", {
  fams <- roy_copula_families
  fits <- pmap(fams, function(f) tryCatch(roy_copula_fit(f_out(), f_sel(), d1, family = f),
                                          error = function(e) NULL))
  names(fits) <- fams
  ok <- !vapply(fits, is.null, logical(1))
  tab <- data.frame(family = fams[ok], logLik = sapply(fits[ok], `[[`, "loglik"),
                    npar = sapply(fits[ok], `[[`, "npar"), AIC = sapply(fits[ok], AIC),
                    BIC = sapply(fits[ok], BIC), conv = sapply(fits[ok], `[[`, "convergence"))
  list(table = tab[order(tab$AIC), ], fits = fits)
})
res$sel_table <- sel$table
best <- sel$table$family[1]
log("AIC best family:", best)
res$sel_te <- lapply(sel$fits[!vapply(sel$fits, is.null, TRUE)], function(f) roy_te(f))
res$sel_dep <- lapply(sel$fits[!vapply(sel$fits, is.null, TRUE)], tau_of)

# ---- PV1-PV10 for the AIC-best family and the Gaussian ---------------------
fam_pv <- unique(c(best, "gaussian"))
res$pv <- list()
for (fam in fam_pv) {
  log("PV loop:", fam)
  fits <- cache(paste0("pv_", fam), pmap(1:10, function(pv) {
    if (pv == 1) return(sel$fits[[fam]])
    roy_copula_fit(f_out(), f_sel(), prep_main(pv), family = fam)
  }))
  res$pv[[fam]] <- list(te = rubin(lapply(fits, roy_te)),
                        dep = rowMeans(sapply(fits, tau_of)),
                        conv = sapply(fits, `[[`, "convergence"))
  if (fam == best) fits_best <- fits
}

fb <- sel$fits[[best]]
# MTE curve (PV1, best family) and LATE for an ICT-at-home shift p25 -> p75
res$mte <- roy_mte(fb, u = seq(0.02, 0.98, by = 0.02))
q <- stats::quantile(d1$ict_home, c(.25, .75))
res$late <- roy_late(fb, "ict_home", from = q[1], to = q[2], data = d1)
res$late_q <- q
# propensity overlap
Fz <- pnorm(drop(fb$Z %*% fb$theta$gamma))
res$overlap <- list(quant_d0 = quantile(Fz[d1$d == 0], c(0, .01, .05, .5, .95, .99, 1)),
                    quant_d1 = quantile(Fz[d1$d == 1], c(0, .01, .05, .5, .95, .99, 1)),
                    share_out = mean(Fz < 0.02 | Fz > 0.98))
res$pscore <- data.frame(p = Fz, d = d1$d)

# ---- heterogeneity (PV1, best family): ATT/ATU on subgroups ----------------
log("heterogeneity")
grp <- list(women = d1$female == 1, men = d1$female == 0,
            age_u40 = d1$age < 40, age_40p = d1$age >= 40,
            tertiary = d1$schooling >= 15, non_tertiary = d1$schooling < 15,
            lit_top_half = d1$skill_lit >= median(d1$skill_lit), lit_bottom_half = d1$skill_lit < median(d1$skill_lit))
for (cn in levels(d1$country)) grp[[paste0("country_", cn)]] <- d1$country == cn
res$het <- do.call(rbind, lapply(names(grp), function(g) {
  t <- roy_te(fb, newdata = d1[grp[[g]], ]); t$group <- g; t$n <- sum(grp[[g]]); t
}))

# ---- robustness (PV1, best family unless noted) ----------------------------
log("robustness")
rob_specs <- list(
  no_occ_fe      = list(data = function() d1, fo = f_out(""), fs = f_sel("")),
  numeracy       = list(data = function() d1, fo = f_out(skill = "skill_num"), fs = f_sel(skill = "skill_num")),
  gpt_exposure   = list(data = function() prep_main(1, exposure_col = "gpt_beta"), fo = f_out(), fs = f_sel()),
  ict_top_half   = list(data = function() prep_main(1, ict_q = 0.5), fo = f_out(), fs = f_sel()),
  t_copula       = list(data = function() d1, fo = f_out(), fs = f_sel(), fam = "t"),
  no_exclusion   = list(data = function() d1, fo = f_out(), fs = update(f_sel(), . ~ . - ict_home))
)
res$rob <- cache("rob", pmap(names(rob_specs), function(nm) {
  r <- rob_specs[[nm]]; dd <- r$data(); fam <- if (is.null(r$fam)) best else r$fam
  fit <- tryCatch(roy_copula_fit(r$fo, r$fs, dd, family = fam), error = function(e) e)
  if (inherits(fit, "error")) return(list(spec = nm, error = conditionMessage(fit)))
  list(spec = nm, family = fam, n = nrow(dd), te = roy_te(fit), dep = tau_of(fit),
       conv = fit$convergence, loglik = fit$loglik)
}))
names(res$rob) <- names(rob_specs)

# ---- secondary: ordered 3-regime model, all occupations (PV1) --------------
log("3-regime ordered model")
d3 <- prep_all(1)
tab3 <- table(d3$occ2, d3$s)
res$J3_n <- nrow(d3); res$J3_shares <- prop.table(table(d3$s))
res$J3 <- cache("J3", pmap(unique(c("gaussian", best)), function(fam) {
  fit <- tryCatch(roy_copula_fit(f_out(""), f_sel("", lhs = "s"), d3, family = fam), error = function(e) e)
  if (inherits(fit, "error")) return(list(family = fam, error = conditionMessage(fit)))
  list(family = fam, AIC = AIC(fit), loglik = fit$loglik, conv = fit$convergence, dep = tau_of(fit),
       te10 = roy_te(fit, to = 1, from = 0), te20 = roy_te(fit, to = 2, from = 0),
       te21 = roy_te(fit, to = 2, from = 1))
}))

saveRDS(res, file.path(out_dir, "piaac_labour_results.rds"))
sink(file.path(out_dir, "piaac_labour_results.txt"))
op <- options(width = 140, digits = 4)
cat("RoyAI Labour Economics application: PIAAC Cycle 2 x AI exposure\n\n")
str(res$sample[c("n", "occ2", "share_d")]); print(res$sample$by_country); print(res$sample$notes)
cat("\n== Descriptives (PV1) ==\n"); print(res$desc)
cat("\n== OLS coefficient on d ==\n"); print(res$ols)
cat("\n== First stage (probit): ict_home ==\n"); print(res$first_stage)
cat("\n== ict_home in regime outcome equations (OLS) ==\n"); print(res$excl_check)
cat("\n== Copula family selection (PV1) ==\n"); print(res$sel_table)
for (f in names(res$sel_te)) { cat("\n--", f, "--\n"); print(res$sel_te[[f]]); print(res$sel_dep[[f]]) }
for (f in names(res$pv)) { cat("\n== Rubin over 10 PVs:", f, "==\n"); print(res$pv[[f]]) }
cat("\n== LATE (ict_home p25 -> p75) ==\n"); print(res$late_q); print(res$late)
cat("\n== Overlap ==\n"); print(res$overlap)
cat("\n== MTE (selected u) ==\n"); print(res$mte[res$mte$u %in% c(.02, .1, .26, .5, .74, .9, .98), ])
cat("\n== Heterogeneity ==\n"); print(res$het)
cat("\n== Robustness ==\n"); print(res$rob)
cat("\n== 3-regime ordered (n =", res$J3_n, ") ==\n"); print(res$J3_shares); print(res$J3)
options(op); sink()
log("done")
