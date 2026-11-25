# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# Paper 3 (distributional): AI-ready work and the wage distribution. PIAAC Cycle 2, six countries,
# all wage earners (AI-ready work is possible only in AI-exposed occupations).
# Counterfactual wage distributions under selection on observables with two polar assumptions on the
# dependence between potential-outcome residuals (rank invariance vs independence), a Shapley
# decomposition of the change in Var(log wage), unconditional QTT, mean-gain sensitivity to the
# credible ATT bounds, and a calibrated CES equilibrium component.
source("./prep2.R")
set.seed(20261125)
OUT <- file.path(ROYAI_ROOT, "RoyAI-data/results/paper3")
FJ <- paste(CTRL, JOB)
gini <- function(y) { x <- sort(exp(y)); n <- length(x); 2 * sum(x * seq_len(n)) / (n * sum(x)) - (n + 1) / n }
ineq <- function(y) { q <- quantile(y, c(.1, .5, .9)); c(var = var(y), gini = gini(y), p90p10 = unname(q[3] - q[1]),
  p50p10 = unname(q[2] - q[1]), p90p50 = unname(q[3] - q[2]), mean = mean(y)) }
# full analysis frame: all wage earners with complete data; exposed sample defines where D can be 1
frame_all <- function(pv = 1) { d <- d0; d$skill_lit <- d[[paste0("lit", pv)]]
  d <- d[!is.na(d$ict_home) & !is.na(d$skill_lit) & !is.na(d$ind), ]
  d$d <- ifelse(d$exposed == 1, d$ict_top, 0L)
  ok <- levels(prep_main(pv)$occ2)                    # exposed occupations with common support (as in Paper 2)
  d <- d[d$exposed == 0 | d$occ2 %in% ok, ]; d$occ2 <- droplevels(d$occ2); d }
core <- function(da, mean_gain = NULL, resid = c("rank","indep"), sigma = Inf, expand = 0.20) {
  resid <- match.arg(resid)
  de <- da[da$exposed == 1, ]; de$occ2 <- droplevels(de$occ2)
  f1 <- lm(as.formula(paste("log_wage ~", FJ)), de[de$d == 1, ]); f0 <- lm(as.formula(paste("log_wage ~", FJ)), de[de$d == 0, ])
  m1 <- predict(f1, de); m0 <- predict(f0, de); cate <- m1 - m0
  if (!is.null(mean_gain)) cate <- cate - mean(cate[de$d == 1]) + mean_gain     # re-centre ATT on a bound
  y <- de$log_wage; D <- de$d
  e <- ifelse(D == 1, y - m1, y - m0)
  if (resid == "indep") { e1 <- sample(e[D == 1], length(e), TRUE); e0 <- sample(e[D == 0], length(e), TRUE)
    Y0 <- ifelse(D == 0, y, m1 - cate + e0); Y1 <- ifelse(D == 1, y, m0 + cate + e1)
  } else { Y0 <- ifelse(D == 0, y, y - cate); Y1 <- ifelse(D == 1, y, y + cate) }
  gain <- Y1 - Y0
  # human-capital and task-complementarity gradients of the gain (projection among adopters)
  g <- lm(gain ~ I(skill_lit) + schooling + skill_num + comp_ind, data = cbind(de, gain = gain), subset = D == 1)
  a <- mean(gain[D == 1]); Xh <- cbind(de$skill_lit, de$schooling, de$skill_num); bh <- coef(g)[2:4]; bc <- coef(g)[5]
  hc <- drop(sweep(Xh, 2, colMeans(Xh[D == 1, ])) %*% bh); tc <- (de$comp_ind - mean(de$comp_ind[D == 1])) * bc
  rest <- gain - a - hc - tc
  ne <- da$log_wage[da$exposed == 0]                       # non-exposed workers: unchanged
  V <- function(yex) var(c(ne, yex))
  # Shapley decomposition of Var(observed) - Var(no AI-ready work)
  comps <- list(mean = D * a, skill = D * hc, task = D * (tc + rest))
  perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  sh <- rowMeans(sapply(perms, function(p) { cur <- Y0; out <- numeric(3)
    for (k in p) { nx <- cur + comps[[k]]; out[k] <- V(nx) - V(cur); cur <- nx }; out }))
  names(sh) <- c("reallocation_mean_gain", "skill_gradient", "task_complementarity_and_residual")
  # counterfactual scenarios
  sc <- list(no_ai = Y0, observed = ifelse(D == 1, Y1, Y0))
  p <- fitted(glm(as.formula(paste("d ~", FJ)), binomial, de))
  nadd <- round(expand * nrow(de)); cand <- which(D == 0); add <- cand[order(-p[cand])][seq_len(min(nadd, length(cand)))]
  Dx <- D; Dx[add] <- 1L
  sc$expand_pe <- ifelse(Dx == 1, Y1, Y0)
  # CES equilibrium: relative price of adopter vs non-adopter efficiency units, elasticity sigma
  if (is.finite(sigma)) { s1 <- mean(D); s1x <- mean(Dx); dln <- log(s1x / (1 - s1x)) - log(s1 / (1 - s1))
    dw <- -(1 / sigma) * dln; sc$expand_ge <- ifelse(Dx == 1, Y1 + (1 - s1x) * dw, Y0 - s1x * dw) } else sc$expand_ge <- sc$expand_pe
  sc$universal <- Y1
  tab <- t(sapply(sc, function(yy) ineq(c(ne, yy))))
  list(tab = tab, shapley = sh, dvar = V(sc$observed) - V(Y0), att = a, grad = coef(g)[-1],
       gain_sd = sd(gain[D == 1]), frac_neg = mean(gain[D == 1] < 0), n = nrow(da), n_exp = nrow(de))
}
res <- list(); da <- frame_all(1)
res$obs_by_group <- sapply(split(da$log_wage, interaction(da$exposed, da$d, drop = TRUE)), ineq)
res$obs_all <- ineq(da$log_wage)
grid <- expand.grid(resid = c("rank","indep"), mg = c("soo", "lo", "hi", "copula"), sigma = c(Inf, 4, 2), stringsAsFactors = FALSE)
mgv <- c(soo = NA, lo = -0.01, hi = 0.06, copula = 0.75)
res$cells <- lapply(seq_len(nrow(grid)), function(i) { g <- grid[i, ]; mg <- if (g$mg == "soo") NULL else mgv[[g$mg]]
  r <- core(da, mg, g$resid, g$sigma); r$cell <- g; r })
# cluster bootstrap for the headline cell (SOO mean, rank invariance, PE)
cb <- function(B = 99) { idx <- split(seq_len(nrow(da)), da$ind)
  t(replicate(B, { s <- sample(names(idx), replace = TRUE); db <- da[unlist(idx[s]), ]; db$occ2 <- droplevels(db$occ2)
    r <- tryCatch(core(db, NULL, "rank"), error = function(e) NULL); if (is.null(r)) rep(NA, 9) else
      c(r$shapley, dvar = r$dvar, att = r$att, var_noai = r$tab["no_ai","var"], var_obs = r$tab["observed","var"],
        dvar_exp = r$tab["expand_pe","var"] - r$tab["observed","var"], gini_d = r$tab["observed","gini"] - r$tab["no_ai","gini"]) })) }
bt <- cb(); res$boot <- rbind(est = apply(bt, 2, mean, na.rm = TRUE), se = apply(bt, 2, sd, na.rm = TRUE),
                              lo = apply(bt, 2, quantile, .025, na.rm = TRUE), hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
# PV robustness of the headline numbers
res$pv <- t(sapply(1:10, function(k) { r <- core(frame_all(k), NULL, "rank"); c(r$shapley, dvar = r$dvar, att = r$att) }))
# unconditional quantile effects on the treated (exposed sample, IPW, job controls)
de <- da[da$exposed == 1, ]; de$occ2 <- droplevels(de$occ2); p <- pmin(pmax(fitted(glm(as.formula(paste("d ~", FJ)), binomial, de)), .01), .99)
w0 <- ifelse(de$d == 0, p / (1 - p), 0); taus <- seq(.1, .9, .1)
wq <- function(y, w, t) { o <- order(y); cw <- cumsum(w[o]) / sum(w); y[o][which(cw >= t)[1]] }
qtt <- function(de, p) { w0 <- ifelse(de$d == 0, p / (1 - p), 0)
  sapply(taus, function(t) quantile(de$log_wage[de$d == 1], t) - wq(de$log_wage, w0, t)) }
res$qtt <- qtt(de, p); idx <- split(seq_len(nrow(de)), de$ind)
bq <- t(replicate(99, { s <- sample(names(idx), replace = TRUE); db <- de[unlist(idx[s]), ]; db$occ2 <- droplevels(db$occ2)
  pb <- pmin(pmax(fitted(glm(as.formula(paste("d ~", FJ)), binomial, db)), .01), .99); qtt(db, pb) }))
res$qtt_ci <- rbind(est = res$qtt, lo = apply(bq, 2, quantile, .025), hi = apply(bq, 2, quantile, .975)); colnames(res$qtt_ci) <- taus
res$raw_qd <- sapply(taus, function(t) quantile(de$log_wage[de$d == 1], t) - quantile(de$log_wage[de$d == 0], t))
saveRDS(res, file.path(OUT, "paper3_results.rds"))
sink(file.path(OUT, "paper3_results.txt")); print(res[c("obs_all","obs_by_group","boot","pv","qtt_ci","raw_qd")], digits = 3)
for (r in res$cells) { cat("\n--", paste(names(r$cell), r$cell, collapse = " "), " ATT", round(r$att, 3), " dVar", round(r$dvar, 4), "\n"); print(round(r$tab, 4)); print(round(r$shapley, 4)) }
sink(); cat("done\n")
