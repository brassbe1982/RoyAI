# Build the SOC-level AI exposure table used by piaac_onet_prep().
# Sources (downloaded 2026-11-25):
#   AIOE (Felten, Raj, Seamans 2021): https://raw.githubusercontent.com/AIOE-Data/AIOE/main/AIOE_DataAppendix.xlsx
#   GPT exposure (Eloundou et al. 2024): https://raw.githubusercontent.com/openai/GPTs-are-GPTs/main/data/occ_level.csv
#   O*NET 31.0: https://www.onetcenter.org/dl_files/database/db_31_0_csv/{work_context,job_zones}.csv
# Needs readxl. Then, with the BLS ISCO-08 <-> SOC 2010 crosswalk (manual download from
# https://www.bls.gov/soc/ : ISCO_SOC_Crosswalk.xls, columns renamed to soc / isco08),
#   onet_to_isco2(ex, crosswalk, c("aioe", "gpt_beta", "complementarity"))
library(RoyAI)
dir <- tempfile("onet"); dir.create(dir)
get <- function(url, f) { p <- file.path(dir, f); utils::download.file(url, p, mode = "wb", quiet = TRUE); p }
aioe <- as.data.frame(readxl::read_excel(
  get("https://raw.githubusercontent.com/AIOE-Data/AIOE/main/AIOE_DataAppendix.xlsx", "aioe.xlsx"),
  sheet = "Appendix A"))
gpt <- read.csv(get("https://raw.githubusercontent.com/openai/GPTs-are-GPTs/main/data/occ_level.csv", "gpt.csv"),
                check.names = FALSE)
base <- "https://www.onetcenter.org/dl_files/database/db_31_0_csv/"
wc <- read.csv(get(paste0(base, "work_context.csv"), "wc.csv"), check.names = FALSE)
jz <- read.csv(get(paste0(base, "job_zones.csv"), "jz.csv"), check.names = FALSE)
ex <- onet_exposure_soc(aioe, gpt, onet_complementarity(wc, jz))
write.csv(ex, "onet_exposure_soc.csv", row.names = FALSE)
