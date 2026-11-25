source("./prep2.R")
d <- prep_main(1); cat("n", nrow(d), "occ", nlevels(d$occ2), "share d", mean(d$d), "\n")
dz <- d[!is.na(d$z_indict), ]; cat("n with industry IV", nrow(dz), "\n")
cl <- ~ ind
f <- function(rhs) as.formula(rhs)
res <- list()
# OLS
o1 <- lm(f(paste("log_wage ~ d +", CTRL)), dz); o2 <- lm(f(paste("log_wage ~ d +", CTRL, JOB)), dz)
ct <- function(m, v="d") coeftest(m, vcov = vcovCL(m, cluster = cl))[v, 1:2]
res$ols <- rbind(base = ct(o1), job = ct(o2)); print(res$ols)
# Oster delta=1 bound with Rmax = 1.3 R2
r0 <- summary(lm(log_wage ~ d, dz))$r.squared; b0 <- coef(lm(log_wage ~ d, dz))["d"]
osterb <- function(m) { b1 <- coef(m)["d"]; r1 <- summary(m)$r.squared; rmax <- min(1, 1.3 * r1)
  c(beta_tilde = unname(b1), beta_star_delta1 = unname(b1 - (b0 - b1) * (rmax - r1) / (r1 - r0))) }
res$oster <- rbind(base = osterb(o1), job = osterb(o2)); print(res$oster)
# first stages & 2SLS for each candidate instrument
ivs <- list(ict_home = "ict_home", ind_newict = "z_indict", fath_exp = "z_fathexp + fath_miss",
            ind_newict_meth = "z_indict + z_indmeth", all = "ict_home + z_indict + z_fathexp + fath_miss")
fsF <- function(ctrl, z) { m1 <- lm(f(paste("d ~", ctrl, "+", z)), dz)
  zz <- trimws(strsplit(z, "\\+")[[1]]); zz <- setdiff(zz, "fath_miss")
  V <- vcovCL(m1, cluster = cl)[zz, zz, drop=FALSE]; b <- coef(m1)[zz]
  c(F = drop(t(b) %*% solve(V) %*% b) / length(zz), coef1 = unname(b[1]), se1 = sqrt(V[1,1])) }
iv1 <- function(ctrl, z) { m <- ivreg(f(paste("log_wage ~ d +", ctrl, "|", ctrl, "+", z)), data = dz)
  s <- summary(m, vcov = vcovCL(m, cluster = cl), diagnostics = TRUE)
  c(late = coef(m)["d"], se = sqrt(vcovCL(m, cluster = cl)["d","d"]),
    sargan_p = if (nrow(s$diagnostics) >= 3) s$diagnostics[3, 4] else NA) }
tab <- do.call(rbind, lapply(names(ivs), function(nm) do.call(rbind, lapply(c("base","job"), function(cc) {
  ctrl <- if (cc == "base") CTRL else paste(CTRL, JOB)
  data.frame(iv = nm, controls = cc, t(fsF(ctrl, ivs[[nm]])), t(iv1(ctrl, ivs[[nm]])))}))))
res$iv <- tab; print(tab, digits = 3)
# placebo: reduced form of each instrument on wages in NON-exposed occupations, conditional on own ICT use
dn <- d0; dn$skill_lit <- dn$lit1; dn <- dn[dn$exposed == 0 & !is.na(dn$z_indict) & !is.na(dn$ict_home), ]
dn <- droplevels(dn)
pl <- lapply(c(CTRL, paste(CTRL, JOB)), function(ctrl) {
  m <- lm(f(paste("log_wage ~ ict_work + ict_home + z_indict + z_fathexp + fath_miss +", ctrl)), dn)
  coeftest(m, vcov = vcovCL(m, cluster = ~ind))[c("ict_home","z_indict","z_fathexp"), 1:2] })
names(pl) <- c("base","job"); res$placebo <- pl; print(pl)
# within-regime direct-effect check in exposed sample (as before)
res$excl <- lapply(c("ict_home","z_indict"), function(z) sapply(0:1, function(k) {
  dk <- dz[dz$d == k, ]; m <- lm(f(paste("log_wage ~", z, "+", CTRL, JOB)), dk); coeftest(m, vcov=vcovCL(m, cluster=dk["ind"]))[z,1:2]}))
names(res$excl) <- c("ict_home","z_indict"); print(res$excl)

# placebo without conditioning on own ICT use at work
res$placebo_uncond <- lapply(c(CTRL, paste(CTRL, JOB)), function(ctrl) {
  m <- lm(f(paste("log_wage ~ ict_home +", ctrl)), dn); coeftest(m, vcov = vcovCL(m, cluster = ~ind))["ict_home", 1:2] })
print(res$placebo_uncond)
# van Kippersluis-Rietveld / Conley et al.: 2SLS with Y - delta*ict_home over a delta grid
dgrid <- seq(0, 0.035, by = 0.0025)
res$conley <- do.call(rbind, lapply(c("base","job"), function(cc) { ctrl <- if (cc=="base") CTRL else paste(CTRL, JOB)
  do.call(rbind, lapply(dgrid, function(dl) { dz$y_adj <- dz$log_wage - dl * dz$ict_home
    m <- ivreg(f(paste("y_adj ~ d +", ctrl, "| ict_home +", ctrl)), data = dz)
    data.frame(controls = cc, delta = dl, late = unname(coef(m)["d"]), se = sqrt(vcovCL(m, cluster = dz$ind)["d","d"])) })) }))
print(res$conley, digits = 3)
saveRDS(list(res = res), "partA.rds")
