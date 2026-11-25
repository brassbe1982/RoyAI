# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# Tables and figures for Paper 2 from paper2_results.rds
P <- file.path(ROYAI_ROOT, "RoyAI-paper2-submission/manuscript_source")
r <- readRDS(file.path(ROYAI_ROOT, "RoyAI-data/results/paper2/paper2_results.rds"))
f3 <- function(x, k = 3) ifelse(is.na(x), "", formatC(x, format = "f", digits = k))
se <- function(x, k = 3) paste0("(", f3(x, k), ")")
w <- function(file, lines) writeLines(lines, file.path(P, "tables", file))
tb <- function(cap, lab, spec, head, body, note, size = "") c(paste0("\\begin{table}[t]\\centering", size, "\\caption{", cap, "}\\label{", lab, "}"),
  paste0("\\begin{adjustbox}{max width=\\textwidth}\\begin{tabular}{", spec, "}\\toprule ", head, "\\\\\\midrule"), body, "\\bottomrule\\end{tabular}\\end{adjustbox}",
  paste0("\\begin{minipage}{0.97\\textwidth}\\footnotesize\\vspace{4pt} ", note, "\\end{minipage}\\end{table}"))
# T1 descriptives
D <- r$desc; lab <- c(log_wage = "Log hourly wage (PPP USD)", skill_lit = "Literacy (std.)", skill_num = "Numeracy (std.)", schooling = "Years of schooling",
  age = "Age", female = "Female", ict_home = "ICT use at home (index)", comp_ind = "Complementarity index", taskdisc = "Task discretion", influence = "Influence", tenure = "Tenure (years)")
w("tab_desc.tex", tb("Descriptive statistics by work arrangement (AI-exposed occupations)", "tab:desc", "lcccc",
  "& Human-only & AI-augmented & AI-dominant & All",
  c(sapply(names(lab), function(k) paste0(lab[k], " & ", paste(f3(D[k, ], 2), collapse = " & "), "\\\\")),
    paste0("Workers & ", paste(format(as.integer(r$n_arr), big.mark = ","), collapse = " & "), " & ", format(sum(r$n_arr), big.mark = ","), "\\\\")),
  "PIAAC Cycle 2, Chile, Czechia, France, Israel, Latvia, Poland; wage and salary earners in the 16 AI-exposed two-digit occupations with common support; literacy plausible value 1. Human-only: ICT use at work below the country's top tercile. Adopters are AI-augmented if their complementarity index (mean of within-country z-scores of task discretion, influence, complex problem solving and learning at work) is at or above the country median among adopters, AI-dominant otherwise."))
# T2 instruments
iv <- r$partA$iv; pl <- r$partA$placebo
nm <- c(ict_home = "ICT use at home", ind_newict = "Industry share reporting new ICT (leave-one-out)", fath_exp = "Father's occupation AI exposure")
plk <- c(ict_home = "ict_home", ind_newict = "z_indict", fath_exp = "z_fathexp")
body <- unlist(lapply(names(nm), function(k) lapply(c("base","job"), function(cc) { x <- iv[iv$iv == k & iv$controls == cc, ]
  p <- pl[[cc]][plk[k], ]
  sprintf("%s & %s & %s & %s %s & %s %s & %s %s\\\\", if (cc == "base") nm[k] else "", if (cc == "base") "Base" else "+ Job", f3(x$F, 1),
    f3(x$late.d), se(x$se), f3(p[1]), se(p[2]), if (k == "ict_home") f3(r$partA$placebo_uncond[[if (cc=="base") 1 else 2]][1]) else "", if (k == "ict_home") se(r$partA$placebo_uncond[[if (cc=="base") 1 else 2]][2]) else "") })))
w("tab_iv.tex", tb("Candidate instruments for AI-ready work: strength and placebo tests", "tab:iv", "llcccc",
  "Instrument & Controls & First-stage $F$ & 2SLS & \\multicolumn{2}{c}{Placebo: effect on wages in non-exposed occupations}\\\\ & & & & given own ICT at work & unconditional", body,
  "First-stage $F$ and standard errors clustered by country $\\times$ two-digit industry. Base controls: literacy, schooling, age, age squared, female, country and two-digit occupation effects. Job controls add firm size, sector, contract type, one-digit industry, tenure, task discretion, influence, learning at work, complex problem solving and repetitive tasks. Placebo: reduced-form coefficient of the instrument in non-exposed occupations, where AI-ready work is zero by construction, so any effect is a direct effect on wages. Sample with a valid industry instrument."))
# T3 bounds
a <- r$aipw; v <- r$vkr; o <- r$ols; os <- r$partA$oster
rows <- c(
 sprintf("OLS, 10 plausible values & %s %s & %s %s\\\\", f3(o["base","est"]), se(o["base","se"]), f3(o["job","est"]), se(o["job","se"])),
 sprintf("Oster bound ($\\delta=1$, $R_{\\max}=1.3R^2$) & %s & %s\\\\", f3(os["base", 2]), f3(os["job", 2])),
 sprintf("Doubly robust ATE & %s %s & %s %s\\\\", f3(a$base["est","ATE"]), se(a$base["se","ATE"]), f3(a$job["est","ATE"]), se(a$job["se","ATE"])),
 sprintf("Doubly robust ATT & %s %s & %s %s\\\\", f3(a$base["est","ATT"]), se(a$base["se","ATT"]), f3(a$job["est","ATT"]), se(a$job["se","ATT"])),
 sprintf("Doubly robust ATU & %s %s & %s %s\\\\", f3(a$base["est","ATU"]), se(a$base["se","ATU"]), f3(a$job["est","ATU"]), se(a$job["se","ATU"])),
 sprintf("ATT $-$ ATU (sorting on observed gains) & %s %s & %s %s\\\\", f3(a$base["est","ATT_ATU"]), se(a$base["se","ATT_ATU"]), f3(a$job["est","ATT_ATU"]), se(a$job["se","ATT_ATU"])),
 "\\midrule",
 sprintf("2SLS, ICT at home, exact exclusion & %s %s & %s %s\\\\", f3(iv$late.d[1]), se(iv$se[1]), f3(iv$late.d[2]), se(iv$se[2])),
 sprintf("2SLS, direct effect from conditional placebo & %s [%s, %s] & %s [%s, %s]\\\\", f3(v$late[1]), f3(v$lo[1], 2), f3(v$hi[1], 2), f3(v$late[3]), f3(v$lo[3], 2), f3(v$hi[3], 2)),
 sprintf("2SLS, direct effect from unconditional placebo & %s [%s, %s] & %s [%s, %s]\\\\", f3(v$late[2]), f3(v$lo[2], 2), f3(v$hi[2], 2), f3(v$late[4]), f3(v$lo[4], 2), f3(v$hi[4], 2)))
w("tab_bounds.tex", tb("The premium for AI-ready work: selection on observables and calibrated instrumental variables", "tab:bounds", "lcc",
 "& Base controls & + Job controls", rows,
 "Log hourly wages; exposed occupations. Standard errors (in parentheses) clustered by country $\\times$ two-digit industry; doubly robust estimates use augmented inverse probability weighting \\citep{robins1994} with a logit propensity score and regime-specific linear outcome models, standard errors from 199 cluster-bootstrap draws. Calibrated 2SLS: outcome $\\ln w-\\hat\\delta\\,\\text{ICT}_{\\text{home}}$ with $\\hat\\delta$ estimated in non-exposed occupations and its sampling error propagated \\citep{vkr2018}; 95 percent intervals in brackets. Plausible value 1 unless stated."))
# T4 heterogeneity
h <- r$het
hr <- function(lbl, m, k) sprintf("%s & %s %s\\\\", lbl, f3(m[k, "est"]), se(m[k, "se"]))
w("tab_het.tex", tb("Heterogeneity of the premium by human capital and task complementarity (job controls)", "tab:het", "lc",
 "Interaction of AI-ready work with & Coefficient", c(hr("Literacy: middle tercile (vs bottom)", h$lit_t, "d:lit_tmid"), hr("Literacy: top tercile (vs bottom)", h$lit_t, "d:lit_thigh"),
   hr("Numeracy (std.)", h$skill_num, "d:skill_num"), hr("Tertiary education", h$tertiary, "d:tertiary"), hr("Complementarity index", h$comp, "d:comp_ind"), hr("Female", h$female, "d:female")),
 "Each row is a separate regression of log hourly wages on AI-ready work, its interaction with the stated characteristic, the characteristic and all base and job controls. Standard errors clustered by country $\\times$ two-digit industry. Plausible value 1."))
# T5 arrangements
m <- r$mn; mo <- r$mn_ols
rr <- function(lbl, kb, ko = NULL) sprintf("%s & %s %s & %s %s & %s %s & %s %s\\\\", lbl,
  if (is.null(ko)) "" else f3(mo$base[ko, "est"]), if (is.null(ko)) "" else se(mo$base[ko, "se"]), f3(m$base["est", kb]), se(m$base["se", kb]),
  if (is.null(ko)) "" else f3(mo$job[ko, "est"]), if (is.null(ko)) "" else se(mo$job[ko, "se"]), f3(m$job["est", kb]), se(m$job["se", kb]))
w("tab_arr.tex", tb("Three work arrangements: AI-augmented and AI-dominant work relative to human-only work", "tab:arr", "lcccc",
 "& \\multicolumn{2}{c}{Base controls} & \\multicolumn{2}{c}{+ Job controls}\\\\\\cmidrule(lr){2-3}\\cmidrule(lr){4-5} & OLS & Doubly robust & OLS & Doubly robust",
 c(rr("Augmented vs human-only (ATE)", "ATE_aug_hum", "arraugmented"), rr("Dominant vs human-only (ATE)", "ATE_dom_hum", "arrdominant"),
   rr("Augmented vs dominant (ATE)", "ATE_aug_dom"), rr("Augmented, effect on the treated", "ATT_aug"), rr("Dominant, effect on the treated", "ATT_dom")),
 "Doubly robust estimates use a multinomial-logit generalized propensity score \\citep{imbens2000,cattaneo2010} and arrangement-specific outcome regressions; 199 cluster-bootstrap draws. OLS: coefficients on arrangement indicators. The job controls include task discretion, influence, learning at work and complex problem solving, which also define the complementarity index; the job-control columns therefore compare workers with similar task content and identify the arrangement effect from the remaining variation."))
# Figure: MTE
gr <- seq(0.05, 0.95, by = 0.05); cols <- c(d0 = "#1f4e79", d017 = "#c0504d", d023 = "#7f7f7f")
pdf(file.path(P, "figures", "fig_mte_p2.pdf"), 6.2, 3.8); par(mar = c(4, 4, 1, 1))
plot(NA, xlim = c(0.05, 0.95), ylim = c(-0.8, 1.2), xlab = expression("Unobserved resistance to AI-ready work, " * u[D]), ylab = "MTE (log points)")
for (k in names(cols)) { x <- r$mte[[k]]; e <- x["est", paste0("u", gr)]; polygon(c(gr, rev(gr)), c(x["lo", paste0("u", gr)], rev(x["hi", paste0("u", gr)])), col = adjustcolor(cols[k], .10), border = NA); lines(gr, e, col = cols[k], lwd = 2) }
abline(h = 0, col = "grey70"); abline(v = r$mte$d0["est", c("q01","q99")], lty = 3, col = "grey50")
legend("topleft", c(expression(delta == 0), expression(delta == 0.017), expression(delta == 0.023)), col = cols, lwd = 2, bty = "n", cex = .8); dev.off()
# MTE summary table
ms <- function(k) { x <- r$mte[[k]]; c(f3(x["est","u0.2"]), f3(x["est","u0.5"]), f3(x["est","u0.8"]), sprintf("%s [%s, %s]", f3(x["est","slope_avg"]), f3(x["lo","slope_avg"], 2), f3(x["hi","slope_avg"], 2))) }
w("tab_mte.tex", tb("Marginal treatment effects of AI-ready work with a calibrated exclusion restriction", "tab:mte", "lcccc",
 "Direct effect of home ICT, $\\delta$ & MTE($u=0.2$) & MTE($u=0.5$) & MTE($u=0.8$) & MTE(0.95) $-$ MTE(0.05) [95\\% CI]",
 sapply(c(d0 = "0 (exact exclusion)", d017 = "0.017 (conditional placebo)", d023 = "0.023 (unconditional placebo)"), function(l) NULL) |> (\(z) mapply(function(k, l) paste0(l, " & ", paste(ms(k), collapse = " & "), "\\\\"), c("d0","d017","d023"), c("0 (exact exclusion)", "0.017 (conditional placebo)", "0.023 (unconditional placebo)")))(),
 "Local instrumental variables \\citep{heckmanvytlacil2005}: $E[\\ln w-\\delta\\,\\text{ICT}_{\\text{home}}\\mid X,P]$ is modelled as linear in all base and job controls, a cubic in the probit propensity score $P$ (instrument: ICT use at home) and interactions of $P$ with literacy, schooling, sex and the complementarity index; the MTE is the derivative in $P$ at the mean of $X$. Common support of $P$: 0.14 to 0.90. Percentile intervals from 199 cluster-bootstrap draws."))
cat("ok\n")
# T7 copula Roy cells (from the credible-result step)
cc <- do.call(rbind, lapply(r$copula_cells, function(c) data.frame(f = c$family, ctl = c$controls, dl = c$delta, aic = c$aic,
  ate = c$te$estimate[1], att = c$te$estimate[2], atu = c$te$estimate[3])))
cc <- cc[order(cc$f, cc$ctl, cc$dl), ]
w("tab_copula.tex", tb("Copula Roy (Model A) estimates when the exclusion restriction is relaxed", "tab:copula", "llcrccc",
 "Copula & Controls & $\\delta$ & AIC & ATE & ATT & ATU",
 sprintf("%s & %s & %s & %s & %s & %s & %s\\\\", tools::toTitleCase(cc$f), ifelse(cc$ctl == "base", "Base", "+ Job"), f3(cc$dl, 3), format(round(cc$aic), big.mark = ","), f3(cc$ate, 2), f3(cc$att, 2), f3(cc$atu, 2)),
 "Full-information maximum likelihood of the binary copula switching model with the \\texttt{RoyAI} R package (multistart). Outcome $\\ln w-\\delta\\,\\text{ICT}_{\\text{home}}$; ICT use at home is excluded from the wage equations. Plausible value 1; $n=8{,}366$. All fits have a positive definite Hessian but a large condition number (about $10^7$) from the occupation effects in the selection equation."))
cat("ok2\n")
