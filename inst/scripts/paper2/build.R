# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
suppressPackageStartupMessages(library(RoyAI))
extra_vars <- c(ict_home="ICTHOMEC2", newict="H2_Q18a02", newmach="H2_Q18a01", newmeth="H2_Q18a03",
  nochange="H2_Q18a07", repet="H2_Q17", firmsize="D2_Q07a", sector="D2_Q03", contract="D2_Q10_C",
  taskdisc="TASKDISCC2_T1", influence="INFLUENCEC2_T1", learnwork="LEARNATWORKC2", complexps="H2_Q06b",
  specsw="F2_Q05e", progr="F2_Q05f", isic2="ISIC2C", isic1n="ISIC1C", fath2="ISCO2_FATH", moth2="ISCO2_MOTH",
  books="J2_Q06", pared="PAREDC2", numwork="NUMWORKC2", writwork="WRITWORKC2", tenure_y="D2_Q05a2", hours="D2_Q11")
keep <- c("SEQID","CNTRYID","SPFWT0","EARNHRBONUSPPPC2","ISCO2C","AGE_R","AGEG10LFS","GENDER_R","YRSQUALC2",
  "ICTWORKC2", unname(extra_vars), "ISIC1C", paste0("PVLIT",1:10), paste0("PVNUM",1:10))
puf <- read_piaac_puf("puf", countries = piaac_analysis_countries()$iso3, columns = keep)
saveRDS(puf, "puf.rds")
ex <- read.csv(file.path(ROYAI_ROOT, "RoyAI-data/exposure_isco08_2d.csv"), colClasses=c(isco08_2d="character"))
d <- piaac_onet_prep(puf, ex, pv = 1, thresholds = list(ict = 2/3), extra = as.list(extra_vars))
# ISIC1C is letters: take raw
raw <- puf[match(paste(d$country, d$id), paste(puf$CNTRYID, puf$SEQID)), ]
d$isic1 <- raw$ISIC1C; d$isic2 <- raw$ISIC2C
for (k in 1:10) { d[[paste0("lit",k)]] <- NA }
lit <- sapply(1:10, function(k) as.numeric(raw[[paste0("PVLIT",k)]]))
mu <- mean(as.numeric(puf$PVLIT1), na.rm=TRUE); s <- sd(as.numeric(puf$PVLIT1), na.rm=TRUE)
for (k in 1:10) d[[paste0("lit",k)]] <- as.numeric(scale(lit[,k]))
saveRDS(d, "d_all.rds")
cat(nrow(d), "\n"); print(summary(d[, c("newict","repet","firmsize","sector","taskdisc","influence","learnwork","complexps","fath2","tenure_y")]))
print(table(d$isic1, useNA="ifany"))
