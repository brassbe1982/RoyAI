# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# WBES firm-level extension of RoyAI: data preparation
# Input: ComprehData6.csv (user's harmonised WBES file, 2006-2026). Uses survey years 2019+.
# D = firm has its own website (the only digital-adoption item asked in every WBES round;
#     a precondition for AI use, NOT AI use itself).
suppressPackageStartupMessages(library(data.table))
wbes_prep <- function(path = file.path(ROYAI_ROOT, "ComprehData6.csv"), from = 2019, min_cell = 5) {
  d <- fread(path)[year >= from]
  d[, cy := paste(CNTname, year)]
  n0 <- nrow(d)
  # 1) flag values that look imputed: equal to the country-year modal value held by >= 5 firms
  flag_mode <- function(x) { t <- table(x); m <- as.numeric(names(t)[which.max(t)]); (x == m) & (max(t) >= 5) }
  d[, `:=`(imp_sales = flag_mode(Sales), imp_lc = flag_mode(laborCost), imp_emp = flag_mode(size_num)), by = cy]
  d <- d[!imp_sales & !imp_lc & Sales > 0 & laborCost > 0 & size_num >= 1 & size_num <= 50000 & !is.na(region) & region != "NA"]
  n1 <- nrow(d)
  # 2) outcomes, expressed relative to the country-year mean (removes currency and price level)
  d[, lp  := log(Sales / size_num)]
  d[, lw  := log(laborCost / size_num)]          # average labour cost per worker (firm wage)
  d[, lsh := log(laborCost / Sales)]             # labour share of sales
  for (v in c("lp", "lw", "lsh")) {
    d[, (v) := get(v) - mean(get(v)), by = cy]
    q <- quantile(d[[v]], c(.01, .99)); d <- d[get(v) >= q[1] & get(v) <= q[2]]
  }
  n2 <- nrow(d)
  # 3) covariates
  d[, `:=`(lsize = log(size_num), lage = log1p(pmin(nyearsOper, 150)), mexp = pmin(MangYrExpSect, 60) / 10,
           female = as.integer(femOwner == 1), manuf = as.integer(sector_MS == "Manufacturing"),
           exporter = as.integer(MainProdMkt == "Intl"), national = as.integer(MainProdMkt == "Natnal"),
           listed = as.integer(grepl("Pub Traded", legalStat)), soleprop = as.integer(grepl("Sole", legalStat)),
           foreign_tech = ForgLicenTechUse, training = formTraining, cert = iCert, audit = extAudit,
           credit = OverDraftFacility, elec_obst = ElectricObstOP)]
  d[, region := factor(region)]
  d[, D := as.integer(web == 1)]
  # 4) instrument: leave-one-out website share among other firms of the same country-year x sampling stratum
  d[, `:=`(ns = .N, sw = sum(D)), by = .(cy, strata)]
  d[, z := ifelse(ns - 1 >= min_cell, (sw - D) / (ns - 1), NA_real_)]
  d[, cyweb := mean(D), by = cy]                     # Mundlak term: z varies within country-year
  d[, cylsize := mean(lsize), by = cy]
  d[, clus := factor(cy)]
  attr(d, "flow") <- c(start = n0, after_imputed_drop = n1, after_trim = n2)
  d
}
CTRL_W <- "lsize + I(lsize^2) + lage + mexp + female + manuf + exporter + national + listed + soleprop + region + cyweb + cylsize"
EXT_W  <- "+ foreign_tech + training + cert + audit + credit + elec_obst"
