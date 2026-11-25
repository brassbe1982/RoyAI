# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# WBES firm-level Roy extension: OLS, Oster bounds, 2SLS with a leave-one-out peer-adoption instrument,
# Conley/plausibly-exogenous sensitivity, AIPW ATT vs ATU, heterogeneity by size and region.
source(file.path(ROYAI_ROOT, "RoyAI-data/results/wbes/scripts/wbes_prep.R"))
suppressPackageStartupMessages({library(sandwich); library(lmtest); library(AER)})
set.seed(20261125)
OUT <- file.path(ROYAI_ROOT, "RoyAI-data/results/wbes")
d <- wbes_prep(); res <- list(flow = attr(d, "flow"))
dz <- d[!is.na(z)]; res$n <- c(all = nrow(d), with_iv = nrow(dz), countries = uniqueN(dz$CNTname), cy = uniqueN(dz$cy))
f <- function(...) as.formula(paste0(...))
ct <- function(m, v, data) { V <- vcovCL(m, cluster = data$clus); c(est = unname(coef(m)[v]), se = sqrt(V[v, v])) }
res$desc <- rbind(dz[, lapply(.SD, mean), by = D, .SDcols = c("lp","lw","lsh","lsize","lage","mexp","female","manuf","exporter","listed","foreign_tech","training","cert","credit","z")][order(D)],
                  dz[, c(D = 9, lapply(.SD, mean)), .SDcols = c("lp","lw","lsh","lsize","lage","mexp","female","manuf","exporter","listed","foreign_tech","training","cert","credit","z")])
res$share_D <- mean(dz$D)
outs <- c(lp = "lp", lw = "lw", lsh = "lsh")
ctrls <- c(base = CTRL_W, ext = paste(CTRL_W, EXT_W))
# OLS + Oster (delta = 1, Rmax = 1.3 R2)
res$ols <- do.call(rbind, lapply(names(outs), function(y) do.call(rbind, lapply(names(ctrls), function(cc) {
  m <- lm(f(y, " ~ D + ", ctrls[cc]), dz); m0 <- lm(f(y, " ~ D"), dz)
  b1 <- coef(m)["D"]; r1 <- summary(m)$r.squared; b0 <- coef(m0)["D"]; r0 <- summary(m0)$r.squared; rmax <- min(1, 1.3 * r1)
  data.frame(y = y, controls = cc, t(ct(m, "D", dz)), raw = unname(b0), oster = unname(b1 - (b0 - b1) * (rmax - r1) / (r1 - r0)), r2 = r1) }))))
# first stage and 2SLS
res$iv <- do.call(rbind, lapply(names(outs), function(y) do.call(rbind, lapply(names(ctrls), function(cc) {
  fs <- lm(f("D ~ z + ", ctrls[cc]), dz); fz <- ct(fs, "z", dz)
  m <- ivreg(f(y, " ~ D + ", ctrls[cc], " | z + ", ctrls[cc]), data = dz)
  rf <- lm(f(y, " ~ z + ", ctrls[cc]), dz)
  data.frame(y = y, controls = cc, fs = fz[1], fs_se = fz[2], F = (fz[1] / fz[2])^2, t(ct(m, "D", dz)), rf = ct(rf, "z", dz)[1]) }))))
# exclusion check: instrument's association with the outcome within each regime
res$within <- do.call(rbind, lapply(names(outs), function(y) do.call(rbind, lapply(0:1, function(k) {
  dk <- dz[D == k]; m <- lm(f(y, " ~ z + ", ctrls["ext"]), dk); data.frame(y = y, D = k, t(ct(m, "z", dk))) }))))
# plausibly exogenous (Conley et al. 2012): y - delta*z, LATE as a function of delta; delta that zeroes the LATE
dgrid <- seq(0, 0.5, by = 0.05)
res$conley <- do.call(rbind, lapply(names(outs), function(y) do.call(rbind, lapply(dgrid, function(dl) {
  dz$ya <- dz[[y]] - dl * dz$z; m <- ivreg(f("ya ~ D + ", ctrls["ext"], " | z + ", ctrls["ext"]), data = dz)
  data.frame(y = y, delta = dl, t(ct(m, "D", dz))) }))))
# AIPW ATE/ATT/ATU (selection on observables, extended controls)
aipw <- function(dd, y, cc) { ps <- glm(f("D ~ ", cc), binomial, dd); p <- pmin(pmax(fitted(ps), .01), .99)
  m1 <- predict(lm(f(y, " ~ ", cc), dd[D == 1]), dd); m0 <- predict(lm(f(y, " ~ ", cc), dd[D == 0]), dd); Y <- dd[[y]]; T <- dd$D
  att <- sum(T * (Y - m0) - (1 - T) * p / (1 - p) * (Y - m0)) / sum(T); atu <- sum((1 - T) * (m1 - Y) + T * (1 - p) / p * (Y - m1)) / sum(1 - T)
  c(ATE = mean(m1 - m0 + T * (Y - m1) / p - (1 - T) * (Y - m0) / (1 - p)), ATT = att, ATU = atu, ATT_ATU = att - atu) }
cb <- function(dd, fun, B = 99) { cl <- unique(dd$clus); idx <- split(seq_len(nrow(dd)), dd$clus)
  t(replicate(B, { s <- sample(names(idx), replace = TRUE); fun(dd[unlist(idx[s])]) })) }
res$aipw <- lapply(names(outs), function(y) { e <- aipw(d, y, ctrls["ext"]); bt <- cb(d, function(x) aipw(x, y, ctrls["ext"]))
  rbind(est = e, se = apply(bt, 2, sd)) }); names(res$aipw) <- names(outs)
# heterogeneity: by size class and region (OLS ext, 2SLS ext)
d[, sizec := cut(size_num, c(0, 19, 99, Inf), labels = c("small","medium","large"))]; dz <- d[!is.na(z)]
res$het <- do.call(rbind, lapply(c("sizec","region"), function(g) do.call(rbind, lapply(levels(factor(dz[[g]])), function(l) {
  dl <- dz[get(g) == l]; cc <- if (g == "region") gsub("region \\+ ", "", ctrls["ext"]) else ctrls["ext"]
  o <- lm(f("lp ~ D + ", cc), dl); iv <- tryCatch(ivreg(f("lp ~ D + ", cc, " | z + ", cc), data = dl), error = function(e) NULL)
  fs <- lm(f("D ~ z + ", cc), dl); fz <- ct(fs, "z", dl)
  data.frame(group = g, level = l, n = nrow(dl), shareD = mean(dl$D), t(ct(o, "D", dl)), F = (fz[1]/fz[2])^2,
             iv = if (is.null(iv)) NA else coef(iv)["D"], iv_se = if (is.null(iv)) NA else sqrt(vcovCL(iv, cluster = dl$clus)["D","D"])) }))))
saveRDS(res, file.path(OUT, "wbes_results.rds"))
sink(file.path(OUT, "wbes_results.txt")); print(res, digits = 3); sink(); cat("done\n")
