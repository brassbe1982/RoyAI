# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# shared prep: analysis frames for the "credible applied result" step
suppressPackageStartupMessages({library(RoyAI); library(AER); library(sandwich); library(lmtest)})
d0 <- readRDS("./d_all.rds")
ex <- read.csv(file.path(ROYAI_ROOT, "RoyAI-data/exposure_isco08_2d.csv"), colClasses=c(isco08_2d="character"))
d0 <- d0[stats::complete.cases(d0[, c("log_wage","skill_lit","skill_num","age","schooling","ict_work","ict_home","ai_exposure")]), ]
d0$country <- factor(d0$country); d0$occ2 <- factor(d0$isco08_2d)
mk01 <- function(x) ifelse(is.na(x), NA, as.integer(x == 1))
d0$newict01 <- mk01(d0$newict); d0$newmeth01 <- mk01(d0$newmeth); d0$newmach01 <- mk01(d0$newmach)
d0$repet01 <- mk01(d0$repet)
d0$isic2[d0$isic2 %in% c("997","998","999","9996","9997","9998","9999",".",".n",".v")] <- NA
d0$isic1[d0$isic1 %in% c("997","998","999")] <- NA
d0$ind <- interaction(d0$country, d0$isic2, drop = TRUE)
# leave-one-out industry x country shares (all wage earners, valid answers)
loo <- function(x, g) { ok <- !is.na(x) & !is.na(g)
  s <- ave(ifelse(ok, x, 0), g, FUN = sum); n <- ave(as.numeric(ok), g, FUN = sum)
  out <- (s - ifelse(ok, x, 0)) / (n - ok); out[is.na(g) | (n - ok) < 10] <- NA; out }
d0$z_indict  <- loo(d0$newict01, d0$ind)
d0$z_indmeth <- loo(d0$newmeth01, d0$ind)
d0$z_indmach <- loo(d0$newmach01, d0$ind)
d0$ind_n <- ave(rep(1, nrow(d0)), d0$ind, FUN = length)
# father's occupation AI exposure
f2 <- sprintf("%02d", as.integer(ifelse(d0$fath2 %in% 0:99, d0$fath2, NA)))
d0$z_fathexp <- ex$aioe[match(f2, ex$isco08_2d)]
d0$fath_miss <- as.integer(is.na(d0$z_fathexp)); d0$z_fathexp[is.na(d0$z_fathexp)] <- 0
# job controls
d0$firmsize_f <- factor(ifelse(is.na(d0$firmsize), 0, d0$firmsize))
d0$sector_f <- factor(ifelse(is.na(d0$sector), 0, d0$sector))
d0$contract_f <- factor(ifelse(is.na(d0$contract), 0, d0$contract))
d0$isic1_f <- factor(ifelse(is.na(d0$isic1), "NA", d0$isic1))
d0$tenure <- pmin(pmax(2023 - d0$tenure_y, 0), 50); d0$tenure[is.na(d0$tenure)] <- median(d0$tenure, na.rm=TRUE)
for (v in c("taskdisc","influence","learnwork","complexps")) d0[[v]][is.na(d0[[v]])] <- mean(d0[[v]], na.rm=TRUE)
d0$repet01[is.na(d0$repet01)] <- 0
# individual complementarity index (within-country z of task discretion, influence, complex PS, learning)
zc <- function(x) ave(x, d0$country, FUN = function(v) (v - mean(v)) / sd(v))
d0$comp_ind <- (zc(d0$taskdisc) + zc(d0$influence) + zc(d0$complexps) + zc(d0$learnwork)) / 4
d0$exposed <- as.integer(d0$ai_exposure > quantile(d0$ai_exposure, 0.5))
d0$ict_top <- as.integer(d0$ict_work > ave(d0$ict_work, d0$country, FUN = function(x) quantile(x, 2/3)))
prep_main <- function(pv = 1, min_cell = 20) {
  d <- d0; d$skill_lit <- d[[paste0("lit", pv)]]
  d <- d[d$exposed == 1 & !is.na(d$ict_home) & !is.na(d$skill_lit), ]
  d$d <- d$ict_top
  tab <- table(d$occ2, d$d); ok <- rownames(tab)[apply(tab, 1, min) >= min_cell]
  droplevels(d[d$occ2 %in% ok, ])
}
CTRL <- "skill_lit + schooling + age + age2 + female + country + occ2"
JOB  <- "+ firmsize_f + sector_f + contract_f + isic1_f + tenure + taskdisc + influence + learnwork + complexps + repet01"
