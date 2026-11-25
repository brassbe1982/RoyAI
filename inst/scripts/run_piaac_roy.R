# First applied RoyAI run: PIAAC Cycle 2 (2023) x O*NET AI exposure.
# Usage: Rscript run_piaac_roy.R <puf_dir> <exposure_isco08_2d.csv> <out_dir>
#   puf_dir: folder with the OECD PUF CSVs prg<cnt>p2.csv; fetch them with
#            download_piaac_puf(dest = puf_dir) from https://webfs.oecd.org/piaac/
# Sample: piaac_analysis_countries() (PUF keeps earnings and 2-digit ISCO).
# Exclusion restriction (assumption to defend): ICT use at home (ICTHOMEC2)
# shifts the propensity to hold AI-ready work but not wages given skills.
library(RoyAI)
args <- commandArgs(trailingOnly = TRUE)
puf_dir <- args[1]; exposure_csv <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

keep <- c("SEQID", "CNTRYID", "SPFWT0", "EARNHRBONUSPPPC2", "ISCO2C", "AGE_R", "AGEG10LFS",
          "GENDER_R", "YRSQUALC2", "ICTWORKC2", "ICTHOMEC2",
          paste0("PVLIT", 1:10), paste0("PVNUM", 1:10), paste0("PVAPS", 1:10))
puf <- read_piaac_puf(puf_dir, countries = piaac_analysis_countries()$iso3, columns = keep)
ex <- read.csv(exposure_csv, colClasses = c(isco08_2d = "character"))
log("PUF rows:", nrow(puf))

rhs <- "skill_lit + schooling + age + I(age^2 / 100) + female + country"
prep <- function(pv, exposure_col = "aioe") {
  d <- piaac_onet_prep(puf, ex, pv = pv, exposure_col = exposure_col, extra = c(ict_home = "ICTHOMEC2"))
  d <- d[stats::complete.cases(d[, c("log_wage", "skill_lit", "age", "schooling", "ict_work", "ict_home")]), ]
  d$country <- factor(d$country); d$occ1 <- factor(substr(d$isco08_2d, 1, 1))
  d
}
fitA <- function(d, extra_rhs = "") roy(as.formula(paste("log_wage ~", rhs, extra_rhs)),
                                        as.formula(paste("d ~", rhs, extra_rhs, "+ ict_home")), data = d)
fitB <- function(d, ...) roy(as.formula(paste("log_wage ~", rhs)),
                             as.formula(paste("s ~", rhs, "+ ict_home")), data = d, ...)

# Rubin's rules over plausible values for the Model A effects
rubin <- function(effs) {
  M <- length(effs); est <- sapply(effs, `[[`, "estimate"); se <- sapply(effs, `[[`, "se")
  qbar <- rowMeans(est); W <- rowMeans(se^2); B <- apply(est, 1, stats::var)
  data.frame(effs[[1]][, 1:2], estimate = qbar, se = sqrt(W + (1 + 1 / M) * B))
}

res <- list()
d1 <- prep(1)
res$sample <- list(n = nrow(d1), by_country = table(d1$country, d1$s), notes = attr(d1, "notes"))
log("sample:", nrow(d1))

log("Model A, PV1-PV10")
cache <- file.path(out_dir, "fitsA_pv.rds")
fitsA <- if (file.exists(cache)) readRDS(cache) else lapply(1:10, function(pv) { log(" PV", pv); fitA(prep(pv)) })
saveRDS(fitsA, cache)
res$A_pv1 <- fitsA[[1]]
res$A_effects_rubin <- rubin(lapply(fitsA, roy_effects))
res$A_par_rubin <- data.frame(parameter = names(fitsA[[1]]$par),
                              estimate = rowMeans(sapply(fitsA, function(f) f$par)))
res$A_diag <- roy_diagnostics(fitsA[[1]])
res$A_lr_independence <- roy_lrtest(roy(as.formula(paste("log_wage ~", rhs)),
                                        as.formula(paste("d ~", rhs, "+ ict_home")), data = d1,
                                        dependence = "independence"), fitsA[[1]])

log("Model A robustness: GPT exposure; 1-digit occupation effects")
res$A_gpt <- fitA(prep(1, "gpt_beta")); res$A_gpt_effects <- roy_effects(res$A_gpt)
# occupation effects need every 1-digit group in both regimes: keep groups
# with at least 30 workers in each (adoption is zero by construction in
# low-exposure groups), i.e. the common-support subsample
tab <- table(d1$occ1, d1$d); cs <- rownames(tab)[apply(tab, 1, min) >= 30]
d_cs <- droplevels(d1[d1$occ1 %in% cs, ])
res$A_occ1_groups <- cs; res$A_occ1_n <- nrow(d_cs)
res$A_occ1 <- fitA(d_cs, "+ occ1"); res$A_occ1_effects <- roy_effects(res$A_occ1)

log("Model B (DMF), bootstrap B = 49")
set.seed(2026)
res$B <- fitB(d1, se = "bootstrap", B = 49)
res$B_effects <- roy_effects(res$B); res$B_diag <- roy_diagnostics(res$B)

for (k in c("A_pv1", "A_gpt", "A_occ1", "B")) res[[k]]$spec$data <- NULL
saveRDS(res, file.path(out_dir, "piaac_roy_results.rds"))
sink(file.path(out_dir, "piaac_roy_results.txt"))
cat("RoyAI first application: PIAAC Cycle 2 x O*NET AI exposure\n\n"); print(res$sample)
cat("\n=== Model A (PV1 fit) ===\n"); print(summary(fitsA[[1]]))
cat("\n=== Model A effects, Rubin's rules over 10 PVs ===\n"); print(res$A_effects_rubin)
print(res$A_diag); print(res$A_lr_independence)
cat("\n=== Robustness: GPT exposure ===\n"); print(res$A_gpt_effects)
cat("\n=== Robustness: + 1-digit occupation effects (groups", paste(cs, collapse = ","), "; n =", nrow(d_cs), ") ===\n"); print(res$A_occ1_effects)
print(round(res$A_occ1$par[c("rho0", "rho1")], 3))
cat("\n=== Model B (DMF two-step, bootstrap SEs) ===\n"); print(summary(res$B)); print(res$B_effects); print(res$B_diag)
sink()
log("done")
