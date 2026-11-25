# Project root holding RoyAI-data/ and the paper folders; set ROYAI_ROOT or run from that root.
ROYAI_ROOT <- Sys.getenv("ROYAI_ROOT", ".")
# Tables and figures for Paper 3 from paper3_results.rds
P <- file.path(ROYAI_ROOT, "RoyAI-paper3-submission/manuscript_source")
r <- readRDS(file.path(ROYAI_ROOT, "RoyAI-data/results/paper3/paper3_results.rds"))
f3 <- function(x, k = 3) formatC(x, format = "f", digits = k)
w <- function(file, lines) writeLines(lines, file.path(P, "tables", file))
# T1 observed inequality
g <- r$obs_by_group; lab <- c(var = "Variance of log wage", gini = "Gini (wage levels)", p90p10 = "P90/P10 (log points)", p50p10 = "P50/P10", p90p50 = "P90/P50", mean = "Mean log wage")
c1 <- r$cells[[1]]
w("tab_obs.tex", c("\\begin{table}[t]\\centering\\caption{Observed hourly-wage inequality, six countries}\\label{tab:obs}",
 "\\begin{adjustbox}{max width=\\textwidth}\\begin{tabular}{lcccc}\\toprule & All & Non-exposed & \\multicolumn{2}{c}{AI-exposed occupations}\\\\\\cmidrule(lr){4-5} & & occupations & $D=0$ & $D=1$\\\\\\midrule",
 sapply(names(lab), function(k) paste0(lab[k], " & ", f3(r$obs_all[k]), " & ", paste(f3(g[k, ]), collapse = " & "), "\\\\")),
 { gc <- readRDS(file.path(ROYAI_ROOT, "RoyAI-data/results/paper3/group_counts.rds")); sprintf("Workers & %s & %s & %s & %s\\\\", format(sum(gc), big.mark = ","), format(gc[["0.0"]], big.mark = ","), format(gc[["1.0"]], big.mark = ","), format(gc[["1.1"]], big.mark = ",")) },
 "\\bottomrule\\end{tabular}\\end{adjustbox}", "\\begin{minipage}{0.92\\textwidth}\\footnotesize\\vspace{4pt} PIAAC Cycle 2 public-use files for Chile, Czechia, France, Israel, Latvia and Poland; wage and salary earners; hourly earnings including bonuses in PPP US dollars, trimmed at 1 percent per tail within country. Exposed occupations: AIOE above the pooled median, restricted to the 16 two-digit occupations with at least 20 workers in each regime. $D=1$: ICT use at work in the country's top tercile (AI-ready work). Literacy plausible value 1.\\end{minipage}\\end{table}"))
# T2 QTT
q <- r$qtt_ci
w("tab_qtt.tex", c("\\begin{table}[t]\\centering\\caption{Unconditional quantile effects of AI-ready work on the treated (exposed occupations)}\\label{tab:qtt}",
 "\\begin{adjustbox}{max width=\\textwidth}\\begin{tabular}{lccccccccc}\\toprule Quantile & 0.1 & 0.2 & 0.3 & 0.4 & 0.5 & 0.6 & 0.7 & 0.8 & 0.9\\\\\\midrule",
 paste0("Raw quantile gap & ", paste(f3(r$raw_qd), collapse = " & "), "\\\\"),
 paste0("QTT (reweighted) & ", paste(f3(q["est", ]), collapse = " & "), "\\\\"),
 paste0("95\\% CI, lower & ", paste(f3(q["lo", ]), collapse = " & "), "\\\\"),
 paste0("95\\% CI, upper & ", paste(f3(q["hi", ]), collapse = " & "), "\\\\"),
 "\\bottomrule\\end{tabular}\\end{adjustbox}", "\\begin{minipage}{0.95\\textwidth}\\footnotesize\\vspace{4pt} QTT: quantile of log wages of adopters minus the same quantile of non-adopters reweighted by $p(X)/(1-p(X))$ \\citep{firpo2007}, with the propensity score from a logit on literacy, schooling, age, age squared, sex, country, two-digit occupation and job controls (firm size, sector, contract, one-digit industry, tenure, task discretion, influence, learning at work, complex problem solving, repetitive tasks). Percentile intervals from 99 bootstrap draws clustered by country $\\times$ two-digit industry. Identification requires selection on observables.\\end{minipage}\\end{table}"))
# T3 Shapley decomposition
b <- r$boot; pv <- r$pv
row <- function(lbl, j) sprintf("%s & %s & (%s, %s) & [%s, %s]\\\\", lbl, f3(c1[[if (j <= 3) "shapley" else "dvar"]][if (j <= 3) j else 1], 4),
  f3(b["lo", j], 4), f3(b["hi", j], 4), f3(min(pv[, j]), 4), f3(max(pv[, j]), 4))
w("tab_shapley.tex", c("\\begin{table}[t]\\centering\\caption{Contribution of AI-ready work to the variance of log wages: Shapley decomposition}\\label{tab:shapley}",
 "\\begin{adjustbox}{max width=\\textwidth}\\begin{tabular}{lccc}\\toprule Component & Estimate & 95\\% CI & Range over 10 PVs\\\\\\midrule",
 row("Reallocation (mean gain of adopters)", 1), row("Skill gradient of gains", 2), row("Task complementarity and residual heterogeneity", 3),
 "\\midrule", row("Total: $\\mathrm{Var}(\\ln Y)-\\mathrm{Var}(\\ln Y^{0})$", 4),
 sprintf("Memo: $\\mathrm{Var}(\\ln Y)$ observed & %s & & \\\\", f3(c1$tab["observed","var"], 4)),
 sprintf("Memo: ATT used & %s & (%s, %s) & \\\\", f3(c1$att), f3(b["lo",5]), f3(b["hi",5])),
 sprintf("Memo: change in Gini & %s & (%s, %s) & \\\\", f3(c1$tab["observed","gini"] - c1$tab["no_ai","gini"], 4), f3(b["lo",9], 4), f3(b["hi",9], 4)),
 "\\bottomrule\\end{tabular}\\end{adjustbox}", "\\begin{minipage}{0.95\\textwidth}\\footnotesize\\vspace{4pt} Population: all 16,751 wage earners in the six countries; only workers in exposed occupations can be adopters. $\\ln Y^0$ is the counterfactual without AI-ready work. Individual gains are $\\hat\\Delta(X)=\\hat m_1(X)-\\hat m_0(X)$ from regime-specific regressions with job controls, under rank invariance of residuals. The three components are Shapley values over the six orderings in which the mean gain, the gain's projection on literacy, numeracy and schooling, and the remainder are switched on. Bootstrap: 99 draws clustered by country $\\times$ two-digit industry.\\end{minipage}\\end{table}"))
# T4 sensitivity grid
cells <- r$cells; key <- sapply(cells, function(x) paste(x$cell$resid, x$cell$mg, x$cell$sigma))
gl <- function(mg, rs, sg) cells[[which(key == paste(rs, mg, sg))]]
mlab <- c(soo = "Selection on observables", lo = "Lower credible bound ($-0.01$)", hi = "Upper credible bound ($0.06$)", copula = "Copula Roy ATT ($0.75$)")
lines <- c("\\begin{table}[t]\\centering\\small\\caption{Distributional counterfactuals under alternative mean gains, residual dependence and equilibrium responses}\\label{tab:grid}",
 "\\begin{adjustbox}{max width=\\textwidth}\\begin{tabular}{llccccc}\\toprule Mean gain of adopters & Residuals & \\multicolumn{2}{c}{Observed vs no AI-ready work} & \\multicolumn{3}{c}{$\\Delta$Var, adoption +20 pp}\\\\\\cmidrule(lr){3-4}\\cmidrule(lr){5-7} & & $\\Delta$Var & $\\Delta$Gini & $\\sigma=\\infty$ & $\\sigma=4$ & $\\sigma=2$\\\\\\midrule")
for (mg in names(mlab)) for (rs in c("rank","indep")) { a <- gl(mg, rs, Inf); b4 <- gl(mg, rs, 4); b2 <- gl(mg, rs, 2)
  dv <- function(x) x$tab["expand_ge","var"] - x$tab["observed","var"]
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s\\\\", if (rs == "rank") mlab[mg] else "", if (rs == "rank") "Rank invariant" else "Independent",
    f3(a$tab["observed","var"] - a$tab["no_ai","var"], 4), f3(a$tab["observed","gini"] - a$tab["no_ai","gini"], 4),
    f3(dv(a), 4), f3(dv(b4), 4), f3(dv(b2), 4))) }
lines <- c(lines, "\\bottomrule\\end{tabular}\\end{adjustbox}", "\\begin{minipage}{0.97\\textwidth}\\footnotesize\\vspace{4pt} Each row re-centres the individual gains $\\hat\\Delta(X)$ so that the mean gain of adopters equals the stated value, keeping their dependence on $X$. Rank invariant: each adopter's residual is the same in both regimes (comonotone potential outcomes). Independent: residuals are drawn independently from the regime-specific residual distributions. Expansion: the 20 percent of exposed workers with the highest propensity among non-adopters switch to AI-ready work. $\\sigma$: elasticity of substitution between adopter and non-adopter labour in a CES aggregate; $\\sigma=\\infty$ is partial equilibrium.\\end{minipage}\\end{table}")
w("tab_grid.tex", lines)
# Figure 1: QTT
pdf(file.path(P, "figures", "fig_qtt.pdf"), 6, 3.6); par(mar = c(4, 4, 1, 1))
tt <- seq(.1, .9, .1); plot(tt, q["est", ], ylim = range(c(q, r$raw_qd)), type = "b", pch = 19, xlab = "Quantile", ylab = "Log points", col = "#1f4e79")
polygon(c(tt, rev(tt)), c(q["lo", ], rev(q["hi", ])), col = adjustcolor("#1f4e79", .15), border = NA)
lines(tt, r$raw_qd, type = "b", pch = 1, lty = 2, col = "grey40"); abline(h = 0, col = "grey70")
legend("topright", c("QTT, reweighted (95% CI)", "Raw quantile gap"), pch = c(19, 1), lty = c(1, 2), col = c("#1f4e79", "grey40"), bty = "n", cex = .8); dev.off()
cat("ok\n")
