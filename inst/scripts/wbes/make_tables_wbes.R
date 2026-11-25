# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# Tables and figure for the WBES firm-level paper
P <- file.path(ROYAI_ROOT, "RoyAI-wbes-extension/manuscript_source")
r <- readRDS(file.path(ROYAI_ROOT, "RoyAI-data/results/wbes/wbes_results.rds"))
f3 <- function(x, k = 3) ifelse(is.na(x), "", formatC(x, format = "f", digits = k))
se <- function(x, k = 3) paste0("(", f3(x, k), ")")
w <- function(file, lines) writeLines(lines, file.path(P, "tables", file))
tb <- function(cap, lab, spec, head, body, note) c(paste0("\\begin{table}[t]\\centering\\caption{", cap, "}\\label{", lab, "}"),
  paste0("\\begin{adjustbox}{max width=\\textwidth}\\begin{tabular}{", spec, "}\\toprule ", head, "\\\\\\midrule"), body, "\\bottomrule\\end{tabular}\\end{adjustbox}",
  paste0("\\begin{minipage}{0.97\\textwidth}\\footnotesize\\vspace{4pt} ", note, "\\end{minipage}\\end{table}"))
d <- as.data.frame(r$desc)
lab <- c(lp = "Log sales per worker (relative)", lw = "Log labour cost per worker (relative)", lsh = "Log labour share of sales (relative)", lsize = "Log employment",
  lage = "Log firm age", mexp = "Manager experience (decades)", female = "Female owner", manuf = "Manufacturing", exporter = "Main market international",
  listed = "Publicly listed", foreign_tech = "Uses foreign-licensed technology", training = "Formal training", cert = "Quality certification", credit = "Overdraft facility",
  z = "Instrument: peer website share")
w("tab_desc.tex", tb("Descriptive statistics by digital regime (firms with a valid instrument)", "tab:desc", "lccc",
  "& No website & Website & All", c(sapply(names(lab), function(k) sprintf("%s & %s & %s & %s\\\\", lab[k], f3(d[d$D == 0, k]), f3(d[d$D == 1, k]), f3(d[d$D == 9, k]))),
  sprintf("Firms & \\multicolumn{3}{c}{%s firms, %s countries, %s country-years; share with website %s}\\\\", format(r$n["with_iv"], big.mark = ","), r$n["countries"], r$n["cy"], f3(r$share_D, 2))),
  "World Bank Enterprise Surveys, standardized survey rounds 2019--2026. Outcomes are deviations from the country-year mean. Values that equal the country-year modal value (held by at least five firms) for sales or labour cost are treated as imputed and dropped; outcomes are trimmed at 1 percent per tail. Instrument: share of other firms in the same country-year and sampling stratum that have a website, for strata with at least five other firms."))
o <- r$ols; iv <- r$iv; a <- r$aipw
yl <- c(lp = "Labour productivity", lw = "Average wage", lsh = "Labour share")
row <- function(y) { ob <- o[o$y == y & o$controls == "base", ]; oe <- o[o$y == y & o$controls == "ext", ]
  ib <- iv[iv$y == y & iv$controls == "base", ]; ie <- iv[iv$y == y & iv$controls == "ext", ]; aa <- a[[y]]
  c(sprintf("%s & %s & %s %s & %s %s & %s & %s %s & %s %s & %s %s\\\\", yl[y], f3(ob$raw), f3(ob$est), se(ob$se), f3(oe$est), se(oe$se), f3(oe$oster),
    f3(aa["est","ATT"]), se(aa["se","ATT"]), f3(aa["est","ATU"]), se(aa["se","ATU"]), f3(ie$est), se(ie$se))) }
w("tab_main.tex", tb("The productivity and wage premium of digital readiness", "tab:main", "lcccccccc",
 "& Raw & OLS base & OLS extended & Oster & \\multicolumn{2}{c}{Doubly robust} & 2SLS extended\\\\\\cmidrule(lr){6-7} & gap & & & bound & ATT & ATU & ",
 c(sapply(names(yl), row), sprintf("First stage (extended): coefficient & \\multicolumn{8}{l}{%s %s, $F=%s$}\\\\", f3(iv$fs[iv$y == "lp" & iv$controls == "ext"]), se(iv$fs_se[iv$y == "lp" & iv$controls == "ext"]), f3(iv$F[iv$y == "lp" & iv$controls == "ext"], 0))),
 "Dependent variables in logs, relative to the country-year mean. Base controls: log employment and its square, log firm age, manager experience, female owner, manufacturing, main market (national, international), listed company, sole proprietorship, region, and the country-year means of website adoption and log employment (Mundlak terms). Extended controls add foreign-licensed technology, formal training, quality certification, external audit, overdraft facility and the electricity obstacle. Oster bound: $\\delta=1$, $R_{\\max}=1.3R^2$. Doubly robust: augmented inverse probability weighting on all 123,175 firms, 99 cluster-bootstrap draws. Other columns: 72,155 firms with a valid instrument. Standard errors clustered by country-year."))
wi <- r$within; cl <- r$conley
dz <- function(y) { x <- cl[cl$y == y, ]; i <- which(diff(sign(x$est)) != 0)[1]; x$delta[i] + (0 - x$est[i]) * (x$delta[i+1] - x$delta[i]) / (x$est[i+1] - x$est[i]) }
pi_ <- iv$fs[iv$y == "lp" & iv$controls == "ext"]
lr <- function(y) { e <- iv$est[iv$y == y & iv$controls == "ext"]; d0 <- wi$est[wi$y == y & wi$D == 0]; d1 <- wi$est[wi$y == y & wi$D == 1]
  sprintf("%s & %s %s & %s %s & %s & %s & %s\\\\", yl[y], f3(d0), se(wi$se[wi$y == y & wi$D == 0]), f3(d1), se(wi$se[wi$y == y & wi$D == 1]), f3(dz(y), 2), f3(e - d1 / pi_), f3(e - d0 / pi_)) }
w("tab_excl.tex", tb("How much direct effect of peer adoption would overturn the instrumental-variable estimate?", "tab:excl", "lccccc",
 "& \\multicolumn{2}{c}{Peer share and outcome, within regime} & Direct effect that & \\multicolumn{2}{c}{2SLS net of direct effect}\\\\\\cmidrule(lr){2-3}\\cmidrule(lr){5-6} & No website & Website & sets 2SLS to zero & $\\delta$ = website & $\\delta$ = no website",
 sapply(names(yl), lr),
 "Columns 1--2: coefficient of the peer website share in regressions of the outcome on the instrument and extended controls, separately by regime; standard errors clustered by country-year. Because the regime is endogenous, these are informative about, not estimates of, the direct effect $\\delta$. Column 3: value of $\\delta$ at which the plausibly-exogenous 2SLS estimate of \\citet{conley2012} equals zero. Columns 4--5: 2SLS minus $\\delta/\\pi$, with $\\pi$ the first-stage coefficient."))
h <- r$het; h$level <- as.character(h$level)
w("tab_het.tex", tb("Heterogeneity of the labour-productivity premium by firm size and region", "tab:het", "lcccccc",
 "Group & Firms & Share with website & OLS extended & First-stage $F$ & 2SLS extended",
 sprintf("%s & %s & %s & %s %s & %s & %s %s\\\\", c(small = "Small (<20)", medium = "Medium (20--99)", large = "Large (100+)", AFR = "Sub-Saharan Africa", EAP = "East Asia and Pacific", ECA = "Europe and Central Asia", LAC = "Latin America and Caribbean", MNA = "Middle East and North Africa", SAR = "South Asia")[h$level],
   format(h$n, big.mark = ","), f3(h$shareD, 2), f3(h$est), se(h$se), f3(h$F, 0), f3(h$iv), se(h$iv_se)),
 "Outcome: log sales per worker relative to the country-year mean. Extended controls as in Table~\\ref{tab:main} (region effects dropped within regions). Standard errors clustered by country-year."))
# copula table (if available)
fs <- list.files(file.path(ROYAI_ROOT, "RoyAI-data/results/wbes/copula"), pattern = "\\.rds$", full.names = TRUE)
if (length(fs)) { cc <- do.call(rbind, lapply(fs, function(f) { x <- readRDS(f); if (!is.null(x$error)) return(NULL)
  tt <- x$coef; data.frame(f = x$family, y = x$y, dl = x$delta, aic = x$aic, ate = x$te$estimate[1], ate_se = x$te$se[1], att = x$te$estimate[2], att_se = x$te$se[2], atu = x$te$estimate[3], atu_se = x$te$se[3],
    tau0 = unname(tt["tau0"]), tau1 = unname(tt["tau1"]), irr = any(grepl("regular", unlist(x$multistart))), min = x$minutes) }))
  cc <- cc[order(cc$y, cc$dl, cc$aic), ]; saveRDS(cc, file.path(ROYAI_ROOT, "RoyAI-data/results/wbes/copula_summary.rds")); print(cc, digits = 3)
  w("tab_copula.tex", tb("Copula Roy model of the website decision: labour productivity and average wage", "tab:copula", "llcrccccc",
   "Outcome & Copula & $\\delta$ & AIC & ATE & ATT & ATU & $\\tau_0$ & $\\tau_1$",
   sprintf("%s & %s & %s & %s & %s %s & %s %s & %s %s & %s & %s\\\\", c(lp = "Productivity", lw = "Wage")[cc$y], tools::toTitleCase(cc$f), f3(cc$dl, 2), format(round(cc$aic), big.mark = ","),
     f3(cc$ate), se(cc$ate_se), f3(cc$att), se(cc$att_se), f3(cc$atu), se(cc$atu_se), f3(cc$tau0, 2), f3(cc$tau1, 2)),
   "Full-information maximum likelihood of the binary copula switching model (\\texttt{RoyAI} R package, multistart), random subsample of 20,000 firms with a valid instrument. Selection equation: extended controls and the peer website share; outcome equations: extended controls. $\\delta$: direct effect of the peer share subtracted from the outcome. $\\tau_s$: Kendall's tau between the outcome error in regime $s$ and the selection error ($\\tau_1>0$: firms with higher productivity in the website regime are more likely to adopt). AIC comparable within outcome and $\\delta$.")) }
cat("ok\n")
