# ---------------------------------------------------------------------------
# Data preparation: PIAAC Cycle 2 public-use file linked to O*NET AI exposure
#
# Variable names, missing codes and suppressions were checked on 2026-11-25
# against the OECD files "piaac-cy2-international-codebook.xlsx" and
# "piaac-cy2-missing_variables-puf.xlsx". The regime definition remains a
# documented proxy: the Cycle 2 background questionnaire has no item on AI use.
# ---------------------------------------------------------------------------

#' PIAAC variable map
#'
#' Names of the PIAAC Cycle 2 public-use-file (PUF) variables used by
#' [piaac_onet_prep()], confirmed against the OECD Cycle 2 international
#' codebook. Each element is a vector of candidates tried in order (the Cycle 2
#' name first, then the Cycle 1 name), so the same code also reads Cycle 1
#' files. The names used are reported in the `variables` attribute of the
#' prepared data. Override any element, e.g. `piaac_varmap(earnings =
#' "EARNHRPPPC2")` to use earnings without bonuses.
#'
#' Notes from the codebook: `ISCO2C` is the 2-digit ISCO-08 code of the
#' current job; `EARNHRBONUSPPPC2` is hourly earnings including bonuses for
#' wage and salary earners in PPP-corrected US dollars; `ICTWORKC2` is the
#' index of ICT use at work; `AGE_R` is suppressed in 15 PUFs, so
#' `AGEG10LFS` (10-year bands) is used as a fallback; plausible values run
#' from 1 to 10; replicate weights are `SPFWT1`-`SPFWT80` with the method in
#' `VEMETHOD`.
#'
#' @param ... Elements to override (character vectors).
#' @return A named list of character vectors.
#' @export
piaac_varmap <- function(...) {
  def <- list(
    id = c("SEQID"),
    country = c("CNTRYID", "CNTRYID_E"),
    weight = c("SPFWT0"),
    earnings = c("EARNHRBONUSPPPC2", "EARNHRBONUSPPP"),
    occupation = c("ISCO2C"),
    age = c("AGE_R"),
    age_group = c("AGEG10LFS"),
    gender = c("GENDER_R"),
    schooling = c("YRSQUALC2", "YRSQUAL"),
    ict_work = c("ICTWORKC2", "ICTWORK"),
    pv_literacy = c("PVLIT"),
    pv_numeracy = c("PVNUM"),
    pv_aps = c("PVAPS")
  )
  utils::modifyList(def, list(...))
}

#' Countries whose Cycle 2 PUF has every variable the RoyAI application needs
#'
#' From the OECD file `piaac-cy2-missing_variables-puf.xlsx`: the countries
#' whose public-use file keeps `ISCO2C`, `EARNHRBONUSPPPC2`, `YRSQUALC2`,
#' `GENDER_R`, `ICTWORKC2` and the plausible values. Earnings are suppressed
#' in 17 PUFs and 2-digit occupation in 12, which leaves ten. The released
#' CSV files (checked 2026-11-25) suppress more than that table says: Croatia
#' and Lithuania code every `ISCO2C` as 999 ("not stated"), and Spain and the
#' Slovak Republic have no valid hourly earnings. That leaves these six.
#'
#' @return A data frame with ISO3 code, numeric `CNTRYID` and country name.
#' @export
piaac_analysis_countries <- function() {
  data.frame(iso3 = c("CHL", "CZE", "FRA", "ISR", "LVA", "POL"),
             cntryid = c(152, 203, 250, 376, 428, 616),
             name = c("Chile", "Czechia", "France", "Israel", "Latvia", "Poland"))
}

# PIAAC Cycle 2 missing-value codes (numeric exports); the CSV PUFs instead
# carry ".n", ".v", ... strings, which become NA on numeric conversion.
.piaac_missing <- list(
  earnings = 999999999996, occupation = 995, age = 99, age_group = 93, gender = 6,
  schooling = 96, ict_work = 9996)

.piaac_num <- function(x, key = NULL) {
  v <- suppressWarnings(as.numeric(as.character(x)))
  if (!is.null(key) && !is.null(.piaac_missing[[key]])) v[v >= .piaac_missing[[key]]] <- NA
  v
}

#' Read PIAAC Cycle 2 public-use CSV files
#'
#' Reads and stacks the country files `prg<cnt>p2.csv` (semicolon separated,
#' as distributed by the OECD) with base R. All columns are read as character
#' so that the PUF missing codes (".n", ".v", ...) survive; [piaac_onet_prep()]
#' converts them to `NA`. Columns absent from a country's file (suppressed
#' variables) are filled with `NA`.
#'
#' @param folder Folder holding the PUF CSV files.
#' @param countries Optional ISO3 codes to read (matched against file names),
#'   e.g. `piaac_analysis_countries()$iso3`.
#' @param columns Optional subset of columns to keep.
#' @return A data frame.
#' @export
read_piaac_puf <- function(folder, countries = NULL, columns = NULL) {
  files <- list.files(folder, pattern = "^prg.*p2\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (!is.null(countries))
    files <- files[tolower(substr(basename(files), 4, 6)) %in% tolower(countries)]
  if (!length(files)) stop("no prg*p2.csv files found in ", folder)
  parts <- lapply(files, function(f) {
    d <- utils::read.csv(f, sep = ";", colClasses = "character", check.names = FALSE)
    if (!is.null(columns)) d <- d[, intersect(columns, names(d)), drop = FALSE]
    d
  })
  all_cols <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(d) { for (m in setdiff(all_cols, names(d))) d[[m]] <- NA_character_; d[, all_cols] })
  do.call(rbind, parts)
}

# resolve candidate names against the columns of a data frame
.resolve_var <- function(cands, cols, what, required = TRUE) {
  hit <- cands[toupper(cands) %in% toupper(cols)]
  if (!length(hit)) {
    if (required) stop("PIAAC variable for '", what, "' not found; tried ",
                       paste(cands, collapse = ", "), ". Set it with piaac_varmap(", what, " = ...)")
    return(NA_character_)
  }
  cols[match(toupper(hit[1]), toupper(cols))]
}

#' Aggregate an O*NET-SOC exposure measure to ISCO-08 2-digit occupations
#'
#' Occupation-level AI-exposure measures are published on the O*NET-SOC /
#' SOC classification (for example the AI Occupational Exposure index of
#' Felten, Raj and Seamans, or the GPT-exposure scores of Eloundou et al.),
#' while PIAAC reports ISCO-08. This helper maps SOC codes to ISCO-08 with a
#' crosswalk (e.g. the BLS ISCO-08 to SOC 2010 crosswalk) and averages to the
#' 2-digit ISCO level, optionally weighting by employment.
#'
#' @param exposure Data frame with a `soc` column and one or more numeric
#'   measure columns.
#' @param crosswalk Data frame with columns `soc` and `isco08`; defaults to
#'   [isco_soc_crosswalk()]. ISCO codes of 3 or more digits;
#'   many-to-many links are allowed.
#' @param measures Names of the measure columns to aggregate.
#' @param weight Optional name of an employment-weight column in `exposure`.
#' @return Data frame with `isco08_2d` and the averaged measures.
#' @export
onet_to_isco2 <- function(exposure, crosswalk = isco_soc_crosswalk(), measures, weight = NULL) {
  stopifnot(all(c("soc") %in% names(exposure)), all(c("soc", "isco08") %in% names(crosswalk)))
  norm_soc <- function(x) gsub("[^0-9]", "", substr(as.character(x), 1, 7))
  exposure$soc_key <- norm_soc(exposure$soc); crosswalk$soc_key <- norm_soc(crosswalk$soc)
  m <- merge(exposure, crosswalk[, c("soc_key", "isco08")], by = "soc_key")
  if (!nrow(m)) stop("no SOC codes matched the crosswalk")
  m$isco08_2d <- substr(gsub("[^0-9]", "", as.character(m$isco08)), 1, 2)
  w <- if (is.null(weight)) rep(1, nrow(m)) else m[[weight]]
  out <- do.call(rbind, lapply(split(seq_len(nrow(m)), m$isco08_2d), function(i) {
    data.frame(isco08_2d = m$isco08_2d[i[1]],
               t(vapply(measures, function(v) stats::weighted.mean(m[[v]][i], w[i], na.rm = TRUE), 1)),
               check.names = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Prepare a PIAAC Cycle 2 + O*NET analysis file for RoyAI
#'
#' Merges an occupation-level AI-exposure table onto PIAAC respondents,
#' builds the outcome (log hourly earnings), human-capital covariates, a
#' skill measure from one plausible value, and the selection variables for
#' Model A (`d`) and Model B (`s`).
#'
#' **Regime definition (a documented limitation).** PIAAC Cycle 2
#' (fieldwork 2022-23) is treated as having no direct question on
#' generative-AI use, so adoption is proxied and the regimes measure
#' AI-ready work rather than observed AI use. Unless `adoption` is supplied:
#' \itemize{
#'   \item Model A: `d = 1` if the occupation's AI exposure is above the
#'     `thresholds$exposure` quantile **and** the worker's ICT use at work is
#'     above the within-country `thresholds$ict` quantile.
#'   \item Model B, following the exposure-complementarity split of
#'     Pizzinelli et al. (2023, IMF): `s = 0` (human-only) if not adopting;
#'     among adopters `s = 1` (AI-augmented) if occupational complementarity
#'     is above the `thresholds$complementarity` quantile, else `s = 2`
#'     (AI-dominant).
#' }
#' Supply `adoption = function(df) <integer 0/1/2 vector>` to use your own
#' definition (for example a survey item on AI tool use, if available).
#'
#' @param piaac Data frame of the PIAAC public-use file (one or more
#'   countries).
#' @param exposure Data frame with `isco08_2d`, the exposure column named by
#'   `exposure_col`, and (for Model B) `complementarity`, e.g. from
#'   [onet_to_isco2()].
#' @param varmap Variable map, see [piaac_varmap()].
#' @param pv Which plausible value to use (1-10). Full analyses should repeat
#'   over all plausible values and combine with Rubin's rules.
#' @param thresholds List with quantiles `exposure`, `ict` and
#'   `complementarity`.
#' @param adoption Optional function returning the regime for each row.
#' @param exposure_col Column of `exposure` used as the AI-exposure measure.
#'   The main specification uses the Felten-Raj-Seamans AI Occupational
#'   Exposure index (`"aioe"`); rerun with the Eloundou et al. GPT-exposure
#'   score (e.g. `"gpt_beta"`) as a robustness check.
#' @param extra Optional named character vector of further PUF variables to
#'   carry along as numeric columns, e.g. `c(ict_home = "ICTHOMEC2")`. PUF
#'   CSV missing codes (".v", ".n", ...) become `NA`.
#' @param trim Share of hourly earnings trimmed at each tail within country
#'   (default 1 percent), against misreporting and extreme values.
#' @param min_wage_ppp Drop hourly earnings below this value (top/bottom
#'   coding and misreporting).
#' @return A data frame with `id`, `country`, `weight`, `log_wage`, `skill_lit`,
#'   `skill_num`, `skill_aps`, `age`, `age2`, `female`, `schooling`, `ict_work`,
#'   `isco08_2d`, `ai_exposure`, `complementarity`, `d`, `s`; attribute `notes`
#'   lists the assumptions used.
#' @export
piaac_onet_prep <- function(piaac, exposure, varmap = piaac_varmap(), pv = 1,
                            thresholds = list(exposure = 0.5, ict = 2 / 3, complementarity = 0.5),
                            adoption = NULL, min_wage_ppp = 1, exposure_col = "aioe",
                            extra = NULL, trim = 0.01) {
  th <- utils::modifyList(list(exposure = 0.5, ict = 2 / 3, complementarity = 0.5), thresholds)
  cols <- names(piaac)
  keys <- c("id", "country", "weight", "earnings", "occupation", "gender", "schooling", "ict_work")
  vn <- stats::setNames(vapply(keys, function(k) .resolve_var(varmap[[k]], cols, k), ""), keys)
  age_v <- .resolve_var(varmap$age, cols, "age", required = FALSE)
  ageg_v <- .resolve_var(varmap$age_group, cols, "age_group", required = FALSE)
  if (is.na(age_v) && is.na(ageg_v)) stop("neither an age nor an age-group variable was found")
  pvn <- vapply(c("pv_literacy", "pv_numeracy", "pv_aps"), function(k)
    .resolve_var(paste0(varmap[[k]], pv), cols, k, required = FALSE), "")
  num <- .piaac_num

  age <- if (!is.na(age_v)) num(piaac[[age_v]], "age") else rep(NA_real_, nrow(piaac))
  if (!is.na(ageg_v)) {
    # AGEG10LFS: 1 = 24 or less, 2 = 25-34, ..., 5 = 55 plus -> band midpoints
    mid <- c(20, 29.5, 39.5, 49.5, 60)[num(piaac[[ageg_v]], "age_group")]
    age[is.na(age)] <- mid[is.na(age)]
  }
  occ <- num(piaac[[vn["occupation"]]], "occupation")
  out <- data.frame(
    id = piaac[[vn["id"]]], country = as.character(piaac[[vn["country"]]]),
    weight = num(piaac[[vn["weight"]]]),
    wage = num(piaac[[vn["earnings"]]], "earnings"),
    age = age,
    female = as.integer(num(piaac[[vn["gender"]]], "gender") == 2),
    schooling = num(piaac[[vn["schooling"]]], "schooling"),
    ict_work = num(piaac[[vn["ict_work"]]], "ict_work"),
    isco08_2d = ifelse(is.na(occ), NA_character_, sprintf("%02d", as.integer(occ))),
    stringsAsFactors = FALSE)
  for (k in 1:3) {
    nm <- c("skill_lit", "skill_num", "skill_aps")[k]
    out[[nm]] <- if (!is.na(pvn[k])) as.numeric(scale(num(piaac[[pvn[k]]]))) else NA_real_
  }
  vn <- c(vn, age = age_v, age_group = ageg_v)
  # PIAAC routes respondents who report no ICT use around the ICT items, so
  # the derived ICT indices are "valid skip" (".v" / 9996) for them. They are
  # recoded to one unit below the lowest observed value (no ICT use) instead
  # of being dropped.
  skip_low <- function(raw, v) {
    rs <- trimws(as.character(raw))
    vs <- rs %in% c(".v", "9996")
    if (any(vs)) v[vs] <- min(v, na.rm = TRUE) - 1
    v
  }
  out$ict_work <- skip_low(piaac[[vn["ict_work"]]], out$ict_work)
  for (k in names(extra)) {
    out[[k]] <- if (extra[[k]] %in% cols) num(piaac[[extra[[k]]]]) else NA_real_
    if (grepl("^ICT", extra[[k]]) && extra[[k]] %in% cols) out[[k]] <- skip_low(piaac[[extra[[k]]]], out[[k]])
  }
  out$age2 <- out$age^2 / 100
  n0 <- nrow(out)
  out <- out[is.finite(out$wage) & out$wage >= min_wage_ppp, , drop = FALSE]
  if (trim > 0) {
    lo <- stats::ave(out$wage, out$country, FUN = function(x) stats::quantile(x, trim))
    hi <- stats::ave(out$wage, out$country, FUN = function(x) stats::quantile(x, 1 - trim))
    out <- out[out$wage >= lo & out$wage <= hi, , drop = FALSE]
  }
  out$log_wage <- log(out$wage)

  ex <- exposure
  if (!exposure_col %in% names(ex)) stop("exposure column '", exposure_col, "' not found")
  ex$exposure <- ex[[exposure_col]]
  ex$isco08_2d <- sprintf("%02d", as.integer(ex$isco08_2d))
  out <- out[!is.na(out$isco08_2d), , drop = FALSE]
  if (!"complementarity" %in% names(ex)) ex$complementarity <- NA_real_
  out <- merge(out, ex[, c("isco08_2d", "exposure", "complementarity")], by = "isco08_2d", all.x = FALSE)
  names(out)[names(out) == "exposure"] <- "ai_exposure"

  if (is.null(adoption)) {
    hi_exp <- out$ai_exposure > stats::quantile(out$ai_exposure, th$exposure, na.rm = TRUE)
    ict_q <- stats::ave(out$ict_work, out$country, FUN = function(x)
      stats::quantile(x, th$ict, na.rm = TRUE))
    adopt <- hi_exp & out$ict_work > ict_q
    hi_comp <- out$complementarity > stats::quantile(out$complementarity, th$complementarity, na.rm = TRUE)
    out$d <- as.integer(adopt)
    out$s <- ifelse(!adopt, 0L, ifelse(hi_comp, 1L, 2L))
  } else {
    reg <- as.integer(adoption(out))
    out$s <- reg; out$d <- as.integer(reg > 0)
  }
  out$wage <- NULL
  out <- out[, c("id", "country", "weight", "log_wage", "skill_lit", "skill_num", "skill_aps",
                 "age", "age2", "female", "schooling", "ict_work", "isco08_2d", "ai_exposure",
                 "complementarity", "d", "s", names(extra))]
  attr(out, "variables") <- c(vn, pvn)
  attr(out, "notes") <- c(
    sprintf("PIAAC variables used: %s", paste(sprintf("%s=%s", names(c(vn, pvn)), c(vn, pvn)), collapse = ", ")),
    sprintf("AI exposure measure: %s", exposure_col),
    sprintf("plausible value %d used; repeat over PV1-PV10 and combine with Rubin's rules", pv),
    sprintf("%d of %d respondents kept (valid earnings, matched occupation, earnings trimmed at %.1f%% per tail within country)", nrow(out), n0, 100 * trim),
    if (is.null(adoption)) c(
      "LIMITATION: adoption is a PROXY (high-exposure occupation x high within-country ICT use at work); it measures AI-ready work, not observed AI use",
      "Model B split by occupational complementarity (Pizzinelli et al. 2023 style)") else
      "adoption supplied by user function",
    "ICT indices: valid skips (no ICT use) recoded to one unit below the observed minimum",
    "PIAAC Cycle 2 names confirmed against the OECD international codebook (2026-11-25)")
  out
}

#' Synthetic PIAAC-shaped data for testing the preparation pipeline
#'
#' Generates fake respondents with the column names of [piaac_varmap()] and a
#' matching fake exposure table. The numbers are synthetic and carry no
#' information about any country.
#'
#' @param n Number of respondents.
#' @param seed Random seed.
#' @return A list with `piaac` and `exposure` data frames.
#' @export
piaac_example <- function(n = 3000, seed = 1) {
  set.seed(seed)
  vm <- lapply(piaac_varmap(), `[`, 1); vm$age_group <- NULL
  occ <- sample(c(11:14, 21:26, 31:35, 41:44, 51:54, 61:63, 71:75, 81:83, 91:96), n, replace = TRUE)
  lit <- stats::rnorm(n, 270, 45)
  p <- data.frame(seq_len(n), sample(c("AAA", "BBB", "CCC"), n, TRUE), stats::runif(n, 50, 300),
                  exp(1.2 + 0.004 * (lit - 270) + 0.05 * (occ < 40) + stats::rnorm(n, 2.5, 0.5)),
                  occ, sample(16:65, n, TRUE), sample(1:2, n, TRUE), sample(6:20, n, TRUE),
                  stats::rnorm(n, 2, 0.7) + 0.4 * (occ < 40))
  names(p) <- c(vm$id, vm$country, vm$weight, vm$earnings, vm$occupation, vm$age, vm$gender,
                vm$schooling, vm$ict_work)
  p[[paste0(vm$pv_literacy, 1)]] <- lit
  p[[paste0(vm$pv_numeracy, 1)]] <- lit + stats::rnorm(n, 0, 20)
  p[[paste0(vm$pv_aps, 1)]] <- lit + stats::rnorm(n, 0, 25)
  codes <- sort(unique(occ))
  ex <- data.frame(isco08_2d = codes,
                   aioe = stats::plogis(2 - codes / 20 + stats::rnorm(length(codes), 0, 0.3)),
                   complementarity = stats::runif(length(codes)))
  ex$gpt_beta <- stats::plogis(stats::qlogis(ex$aioe) + stats::rnorm(length(codes), 0, 0.5))
  list(piaac = p, exposure = ex)
}

#' Occupational complementarity to AI from O*NET work context
#'
#' Builds an occupation-level complementarity score in the spirit of
#' Pizzinelli, Panton, Mendes Tavares, Cazzaniga and Li (2023, IMF Working
#' Paper 23/216): occupations whose tasks involve communication, responsibility
#' for others, physical presence, critical decisions, unstructured work and
#' high skill requirements are treated as shielded by, and complementary to,
#' AI. Each O*NET work-context item (scale "CX", 1-5) and the Job Zone (1-5)
#' are rescaled to \[0, 1\]; items are averaged within six components and the
#' score is the mean of the components. This is an **approximation** of the
#' published index: the item list follows the paper's components as closely as
#' O*NET 31.0 allows, and routine items enter reversed.
#'
#' @param work_context Data frame read from O*NET `work_context.csv`.
#' @param job_zones Data frame read from O*NET `job_zones.csv`.
#' @return Data frame with `soc` (6-digit, "XX-XXXX"), the six components and
#'   `complementarity`.
#' @export
onet_complementarity <- function(work_context, job_zones) {
  wc <- work_context[work_context[["Scale ID"]] == "CX", ]
  comp <- list(
    communication = c("Face-to-Face Discussions with Individuals and Within Teams", "Public Speaking"),
    responsibility = c("Health and Safety of Other Workers", "Work Outcomes and Results of Other Workers"),
    physical = c("Outdoors, Exposed to All Weather Conditions", "Physical Proximity"),
    criticality = c("Consequence of Error", "Freedom to Make Decisions", "Frequency of Decision Making"),
    unstructured = c("Determine Tasks, Priorities and Goals", "-Degree of Automation",
                     "-Importance of Repeating Same Tasks"))
  wc$soc <- substr(wc[["O*NET-SOC Code"]], 1, 7)
  wc$val <- (as.numeric(wc[["Data Value"]]) - 1) / 4
  agg <- stats::aggregate(val ~ soc + `Element Name`, data = wc, FUN = mean)
  names(agg)[2] <- "item"
  wide <- stats::reshape(agg, idvar = "soc", timevar = "item", direction = "wide")
  names(wide) <- sub("^val\\.", "", names(wide))
  out <- data.frame(soc = wide$soc)
  for (k in names(comp)) {
    items <- comp[[k]]; rev <- startsWith(items, "-"); items <- sub("^-", "", items)
    miss <- setdiff(items, names(wide))
    if (length(miss)) stop("O*NET work-context items not found: ", paste(miss, collapse = "; "))
    M <- as.matrix(wide[, items, drop = FALSE]); M[, rev] <- 1 - M[, rev]
    out[[k]] <- rowMeans(M, na.rm = TRUE)
  }
  jz <- job_zones
  jz$soc <- substr(jz[["O*NET-SOC Code"]], 1, 7)
  jz <- stats::aggregate(`Job Zone` ~ soc, data = jz, FUN = mean)
  out$skills <- (jz[["Job Zone"]][match(out$soc, jz$soc)] - 1) / 4
  out$complementarity <- rowMeans(out[, c(names(comp), "skills")], na.rm = TRUE)
  out
}

#' Assemble SOC-level AI exposure measures
#'
#' Merges the Felten-Raj-Seamans AIOE (6-digit SOC 2010), the Eloundou et
#' al. GPT-exposure ratings (O*NET-SOC 2019, averaged to 6 digits) and the
#' complementarity score from [onet_complementarity()] into one table ready
#' for [onet_to_isco2()].
#'
#' @param aioe Data frame with columns `SOC Code` and `AIOE` (Appendix A of
#'   `AIOE_DataAppendix.xlsx`).
#' @param gpt Data frame read from `data/occ_level.csv` of the GPTs-are-GPTs
#'   repository (`O*NET-SOC Code`, `dv_rating_beta`, `human_rating_beta`).
#' @param complementarity Output of [onet_complementarity()].
#' @param soc_xwalk Optional SOC 2010-2018 crosswalk from
#'   [read_soc_2010_2018()]. When supplied, the GPT and complementarity
#'   measures (SOC 2018 based) are converted to SOC 2010 to match the AIOE and
#'   the ISCO crosswalk; otherwise codes are matched as they are.
#' @return Data frame with `soc`, `aioe`, `gpt_beta`, `gpt_beta_human`,
#'   `complementarity`.
#' @export
onet_exposure_soc <- function(aioe, gpt, complementarity, soc_xwalk = NULL) {
  a <- data.frame(soc = substr(aioe[["SOC Code"]], 1, 7), aioe = as.numeric(aioe[["AIOE"]]))
  g <- data.frame(soc = substr(gpt[["O*NET-SOC Code"]], 1, 7),
                  gpt_beta = as.numeric(gpt[["dv_rating_beta"]]),
                  gpt_beta_human = as.numeric(gpt[["human_rating_beta"]]))
  g <- stats::aggregate(cbind(gpt_beta, gpt_beta_human) ~ soc, data = g, FUN = mean)
  cm <- complementarity[, c("soc", "complementarity")]
  if (!is.null(soc_xwalk)) {
    g <- soc18_to_soc10(g, soc_xwalk, c("gpt_beta", "gpt_beta_human"))
    cm <- soc18_to_soc10(cm, soc_xwalk, "complementarity")
  }
  out <- merge(a, g, by = "soc", all = TRUE)
  out <- merge(out, cm, by = "soc", all = TRUE)
  out[order(out$soc), ]
}

#' ISCO-08 to SOC 2010 crosswalk shipped with RoyAI
#'
#' Reads `inst/extdata/isco08_soc10.csv`: a one-to-many crosswalk from 3-digit
#' ISCO-08 groups to detailed SOC 2010 occupations, taken from the
#' `iscoCrosswalks` package (MIT licence; IBS, after Hardy, Keister and
#' Lewandowski 2018, built on the BLS ISCO-08/SOC 2010 crosswalk). See
#' `isco08_soc10_SOURCE.txt` in the same folder. The official BLS file can be
#' used instead by passing it to [onet_to_isco2()].
#'
#' @return Data frame with `isco08`, `soc` (SOC 2010, "XX-XXXX"), and labels.
#' @export
isco_soc_crosswalk <- function() {
  x <- utils::read.csv(system.file("extdata", "isco08_soc10.csv", package = "RoyAI"),
                       colClasses = "character")
  x$soc <- paste0(substr(x$soc10, 1, 2), "-", substr(x$soc10, 3, 6))
  x
}

#' Convert SOC 2018 occupation measures to SOC 2010
#'
#' O*NET 31.0 and the GPT-exposure ratings use O*NET-SOC 2019 codes (built on
#' SOC 2018), while the AIOE and the ISCO crosswalk use SOC 2010. This
#' averages SOC 2018 measures onto SOC 2010 codes with the BLS
#' `soc_2010_to_2018_crosswalk.xlsx`.
#'
#' @param x Data frame with `soc` (SOC 2018, "XX-XXXX") and measure columns.
#' @param xwalk Data frame with columns `soc2010` and `soc2018`.
#' @param measures Measure columns to convert.
#' @return Data frame with `soc` (SOC 2010) and the averaged measures.
#' @export
soc18_to_soc10 <- function(x, xwalk, measures) {
  m <- merge(xwalk[, c("soc2010", "soc2018")], x[, c("soc", measures)], by.x = "soc2018", by.y = "soc")
  out <- stats::aggregate(m[, measures, drop = FALSE], by = list(soc = m$soc2010),
                          FUN = function(v) mean(v, na.rm = TRUE))
  out[order(out$soc), ]
}

#' Read the BLS SOC 2010 to 2018 crosswalk
#'
#' @param path Path to `soc_2010_to_2018_crosswalk.xlsx` from BLS.
#' @return Data frame with `soc2010` and `soc2018`.
#' @export
read_soc_2010_2018 <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) stop("install the 'readxl' package")
  raw <- as.data.frame(readxl::read_excel(path, col_names = FALSE, .name_repair = "minimal"))
  hdr <- which(raw[[1]] == "2010 SOC Code")[1]
  if (is.na(hdr)) stop("header row '2010 SOC Code' not found")
  d <- raw[(hdr + 1):nrow(raw), c(1, 3)]
  names(d) <- c("soc2010", "soc2018")
  d[stats::complete.cases(d) & grepl("^[0-9]{2}-[0-9]{4}$", d$soc2010), ]
}

#' Download PIAAC Cycle 2 public-use CSV files from the OECD
#'
#' Fetches `prg<cnt>p2.csv` (or the zipped CSV where the OECD distributes a
#' zip) from the OECD PIAAC public-use file server
#' (\url{https://webfs.oecd.org/piaac/index.html}). By downloading you accept
#' the OECD terms of use for the PIAAC data.
#'
#' @param countries ISO3 codes; defaults to [piaac_analysis_countries()].
#' @param dest Destination folder.
#' @param base Base URL of the Cycle 2 CSV folder.
#' @return Paths of the CSV files, invisibly.
#' @export
download_piaac_puf <- function(countries = piaac_analysis_countries()$iso3, dest = ".",
                               base = "https://webfs.oecd.org/piaac/cy2-puf-data/CSV/") {
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  old <- options(timeout = max(1800, getOption("timeout"))); on.exit(options(old))
  out <- character(0)
  for (cc in tolower(countries)) {
    f <- file.path(dest, sprintf("prg%sp2.csv", cc))
    if (file.exists(f)) { out <- c(out, f); next }
    ok <- tryCatch(utils::download.file(paste0(base, basename(f)), f, mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) {
      if (file.exists(f)) file.remove(f)
      z <- tempfile(fileext = ".zip")
      utils::download.file(paste0(base, sprintf("prg%sp2_csv.zip", cc)), z, mode = "wb", quiet = TRUE)
      utils::unzip(z, files = basename(f), exdir = dest)
    }
    out <- c(out, f)
  }
  invisible(out)
}
