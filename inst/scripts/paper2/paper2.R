# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# Paper 2 (applied): who gains from AI-ready work? PIAAC Cycle 2, AI-exposed occupations.
# Selection-on-observables bounds, heterogeneity by human capital, worker-level three-arrangement
# effects (multinomial AIPW) and a local-IV MTE with calibrated exclusion. Needs prep2.R + d_all.rds.
source("./prep2.R")
suppressPackageStartupMessages(library(nnet))
set.seed(20261125)
OUT <- file.path(ROYAI_ROOT, "RoyAI-data/results/paper2")
B <- 199
rubin <- function(est, se) { m <- length(est); qb <- mean(est); W <- mean(se^2); Bv <- var(est)
  c(est = qb, se = sqrt(W + (1 + 1/m) * Bv)) }
addarr <- function(d) { med <- ave(ifelse(d$d == 1, d$comp_ind, NA), d$country, FUN = function(x) median(x, na.rm = TRUE))
  d$arr <- factor(ifelse(d$d == 0, "human", ifelse(d$comp_ind >= med, "augmented", "dominant")), levels = c("human","augmented","dominant")); d }
# outcome regressions on a common design matrix (robust to factor levels missing in a subgroup)
mfit <- function(d, ctrl, idx) { X <- model.matrix(as.formula(paste("~", ctrl)), d); y <- d$log_wage
  b <- lm.fit(X[idx, , drop = FALSE], y[idx])$coefficients; b[is.na(b)] <- 0; drop(X %*% b) }
res <- list()
d1 <- addarr(prep_main(1)); d1 <- droplevels(d1[!is.na(d1$ind), ])
# ---- A. descriptives by arrangement
vars <- c("log_wage","skill_lit","skill_num","schooling","age","female","ict_home","comp_ind","taskdisc","influence","tenure")
res$desc <- sapply(split(d1[, vars], d1$arr), function(x) colMeans(x, na.rm = TRUE))
res$desc <- cbind(res$desc, all = colMeans(d1[, vars], na.rm = TRUE)); res$n_arr <- table(d1$arr)
# ---- B. OLS premium over 10 PVs (base, job), cluster country x ISIC2
ols_pv <- function(ctrl) { r <- sapply(1:10, function(k) { d <- prep_main(k); d <- droplevels(d[!is.na(d$ind), ])
  m <- lm(as.formula(paste("log_wage ~ d +", ctrl)), d); c(coef(m)["d"], sqrt(vcovCL(m, cluster = ~ind)["d","d"])) })
  rubin(r[1,], r[2,]) }
res$ols <- rbind(base = ols_pv(CTRL), job = ols_pv(paste(CTRL, JOB)))
# ---- C. Heterogeneity in the premium (job controls, PV1)
FJ <- paste(CTRL, JOB)
d1$lit_t <- cut(d1$skill_lit, quantile(d1$skill_lit, 0:3/3), include.lowest = TRUE, labels = c("low","mid","high"))
d1$tertiary <- as.integer(d1$schooling >= 15)
het <- function(term) { m <- lm(as.formula(paste("log_wage ~ d *", term, "+", FJ)), d1); V <- vcovCL(m, cluster = ~ind)
  k <- grep("^d$|^d:", names(coef(m)), value = TRUE); cbind(est = coef(m)[k], se = sqrt(diag(V)[k])) }
res$het <- list(lit_t = het("lit_t"), comp = het("comp_ind"), tertiary = het("tertiary"), female = het("female"),
                skill_num = het("skill_num"))
# ---- D. AIPW binary: ATE / ATT / ATU and the sorting test ATT - ATU (selection on observables)
aipw_bin <- function(d, ctrl) {
  ps <- glm(as.formula(paste("d ~", ctrl)), binomial, d); p <- pmin(pmax(fitted(ps), 0.01), 0.99)
  m1 <- mfit(d, ctrl, d$d == 1); m0 <- mfit(d, ctrl, d$d == 0); y <- d$log_wage; D <- d$d
  ate <- mean(m1 - m0 + D * (y - m1) / p - (1 - D) * (y - m0) / (1 - p))
  att <- sum(D * (y - m0) - (1 - D) * p / (1 - p) * (y - m0)) / sum(D)
  atu <- sum((1 - D) * (m1 - y) + D * (1 - p) / p * (y - m1)) / sum(1 - D)
  c(ATE = ate, ATT = att, ATU = atu, ATT_ATU = att - atu, trim_share = mean(fitted(ps) < .01 | fitted(ps) > .99)) }
cboot <- function(d, fun, B) { idx <- split(seq_len(nrow(d)), as.character(d$ind)); t(replicate(B, { s <- sample(names(idx), replace = TRUE)
  db <- droplevels(d[unlist(idx[s]), , drop = FALSE])
  tryCatch(fun(db), error = function(e) rep(NA, length(fun(d)))) })) }
dd <- d1[!is.na(d1$ind), ]; dd <- droplevels(dd)
ctrlB <- gsub("occ2", "occ2", FJ)
res$aipw <- lapply(c(base = CTRL, job = FJ), function(cc) { est <- aipw_bin(dd, cc); bt <- cboot(dd, function(x) aipw_bin(x, cc), B)
  rbind(est = est, se = apply(bt, 2, sd, na.rm = TRUE)) })
res$aipw_pv <- sapply(1:10, function(k) { d <- prep_main(k); d <- droplevels(d[!is.na(d$ind), ]); aipw_bin(d, FJ) })
# ---- E. Three arrangements: multinomial AIPW (generalised propensity score)
aipw_mn <- function(d, ctrl) {
  ps <- multinom(as.formula(paste("arr ~", ctrl)), d, trace = FALSE, MaxNWts = 5000); P <- pmax(fitted(ps), 0.01); P <- P / rowSums(P)
  y <- d$log_wage; mu <- sapply(levels(d$arr), function(a) mfit(d, ctrl, d$arr == a))
  psi <- sapply(1:3, function(j) mu[, j] + (d$arr == levels(d$arr)[j]) * (y - mu[, j]) / P[, j])
  att <- function(j, k) { T <- d$arr == levels(d$arr)[j]; K <- d$arr == levels(d$arr)[k]
    (sum(T * (y - mu[, k])) - sum(K * P[, j] / P[, k] * (y - mu[, k]))) / sum(T) }
  c(ATE_aug_hum = mean(psi[,2] - psi[,1]), ATE_dom_hum = mean(psi[,3] - psi[,1]), ATE_aug_dom = mean(psi[,2] - psi[,3]),
    ATT_aug = att(2, 1), ATT_dom = att(3, 1)) }
res$mn <- lapply(c(base = CTRL, job = FJ), function(cc) { est <- aipw_mn(dd, cc); bt <- cboot(dd, function(x) aipw_mn(x, cc), B)
  rbind(est = est, se = apply(bt, 2, sd, na.rm = TRUE)) })
res$mn_ols <- lapply(c(base = CTRL, job = FJ), function(cc) { m <- lm(as.formula(paste("log_wage ~ arr +", cc)), dd)
  k <- c("arraugmented","arrdominant"); cbind(est = coef(m)[k], se = sqrt(diag(vcovCL(m, cluster = ~ind))[k])) })
# ---- F. Local-IV MTE (Heckman-Vytlacil), instrument ict_home, outcome adjusted by delta * ict_home
mte_fit <- function(d, ctrl, delta, grid = seq(0.05, 0.95, by = 0.05)) {
  pr <- glm(as.formula(paste("d ~", ctrl, "+ ict_home")), binomial(link = "probit"), d); P <- fitted(pr)
  d$P <- P; d$P2 <- P^2; d$P3 <- P^3; d$ya <- d$log_wage - delta * d$ict_home
  X <- c("skill_lit","schooling","female","comp_ind")
  for (x in X) d[[paste0("Px_", x)]] <- P * (d[[x]] - mean(d[[x]]))
  m <- lm(as.formula(paste("ya ~", ctrl, "+ P + P2 + P3 +", paste0("Px_", X, collapse = "+"))), d)
  b <- coef(m); mte <- b["P"] + 2 * b["P2"] * grid + 3 * b["P3"] * grid^2
  supp <- c(q01 = unname(quantile(P[d$d == 1], .01)), q99 = unname(quantile(P[d$d == 0], .99)))
  c(setNames(mte, paste0("u", grid)), slope_avg = unname(mte[length(mte)] - mte[1]), supp) }
res$mte <- lapply(c(d0 = 0, d017 = 0.017, d023 = 0.023), function(dl) {
  est <- mte_fit(dd, FJ, dl); bt <- cboot(dd, function(x) mte_fit(x, FJ, dl), B)
  rbind(est = est, lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE)) })
# ---- G. results from the credible-result step (instrument diagnostics, Oster, vKR) carried over
res$partA <- readRDS(file.path(ROYAI_ROOT, "RoyAI-data/results/credible/partA.rds"))$res[c("ols","oster","iv","placebo","placebo_uncond")]
res$vkr <- readRDS(file.path(ROYAI_ROOT, "RoyAI-data/results/credible/vkr.rds"))
res$copula_cells <- lapply(list.files(file.path(ROYAI_ROOT, "RoyAI-data/results/credible/cells"), full.names = TRUE), readRDS)
saveRDS(res, file.path(OUT, "paper2_results.rds"))
sink(file.path(OUT, "paper2_results.txt")); print(res[setdiff(names(res), c("copula_cells","partA"))], digits = 3); sink()
cat("done\n")
