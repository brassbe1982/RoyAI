# Build LaTeX tables and PDF figures for the Labour Economics manuscript from
# the saved results of run_piaac_labour.R and run_piaac_sensitivity.R.
# Usage: Rscript make_tables_figures.R <out_dir> <tables_dir> <figures_dir>
suppressPackageStartupMessages(library(RoyAI))
a <- commandArgs(trailingOnly = TRUE); od <- a[1]; td <- a[2]; fd <- a[3]
dir.create(td, FALSE, TRUE); dir.create(fd, FALSE, TRUE)
res <- readRDS(file.path(od, "piaac_labour_results.rds"))
sel <- readRDS(file.path(od, "sel_pv1.rds"))
rob <- readRDS(file.path(od, "rob.rds")); names(rob) <- sapply(rob, `[[`, "spec")
sens <- readRDS(file.path(od, "sensitivity_delta.rds"))
f2 <- function(x, d = 3) formatC(x, format = "f", digits = d)
se <- function(x) paste0("(", f2(x), ")")
wr <- function(lines, f) writeLines(lines, file.path(td, f))
fam_lab <- c(gaussian = "Gaussian", t = "Student-$t$", clayton = "Clayton", gumbel = "Gumbel", frank = "Frank")

# Table 1: descriptives
ds <- res$desc
lab <- c(log_wage = "Log hourly earnings (PPP USD)", skill_lit = "Literacy (PV1, standardised)",
         skill_num = "Numeracy (PV1, standardised)", schooling = "Years of schooling", age = "Age",
         female = "Female", ict_work = "ICT use at work (index)", ict_home = "ICT use at home (index)",
         ai_exposure = "AI occupational exposure (AIOE)")
n0 <- sum(res$sample$by_country[, "0"]); n1 <- sum(res$sample$by_country[, "1"])
wr(c(sprintf("%s & %s & %s & %s & %s\\\\", lab[ds$variable], f2(ds$mean_d0, 2), f2(ds$mean_d1, 2),
             f2(ds$mean_d1 - ds$mean_d0, 2), f2(ds$sd_all, 2)),
     "\\midrule", sprintf("Workers & %s & %s & & %s\\\\", format(n0, big.mark = ","), format(n1, big.mark = ","),
                          format(n0 + n1, big.mark = ","))), "tab_desc.tex")

# Table 2: family comparison, PV1
tab <- sel$table
rows <- sapply(tab$family, function(f) {
  te <- res$sel_te[[f]]; dp <- res$sel_dep[[f]]
  sprintf("%s & %s & %s & %s & %s & %s (%s) & %s (%s) & %s (%s)\\\\", fam_lab[f],
          format(round(tab$AIC[tab$family == f], 1), nsmall = 1, big.mark = ","),
          f2(tab$AIC[tab$family == f] - min(tab$AIC), 1),
          f2(dp[["tau0"]], 2), f2(dp[["tau1"]], 2),
          f2(te$estimate[1]), f2(te$se[1]), f2(te$estimate[2]), f2(te$se[2]), f2(te$estimate[3]), f2(te$se[3]))
})
wr(rows, "tab_families.tex")

# Table 3: main estimates (Rubin over 10 PVs) with OLS benchmark
pv <- res$pv
ols <- res$ols
line <- function(nm, k) sprintf("%s & %s & %s & %s\\\\", nm, f2(pv$frank$te$estimate[k]),
                                f2(pv$gaussian$te$estimate[k]), if (k == 1) f2(ols[["Estimate"]]) else "")
lineS <- function(k) sprintf(" & %s & %s & %s\\\\", se(pv$frank$te$se[k]), se(pv$gaussian$te$se[k]),
                             if (k == 1) se(ols[["Std. Error"]]) else "")
wr(c(line("ATE", 1), lineS(1), line("ATT", 2), lineS(2), line("ATU", 3), lineS(3), "\\midrule",
     sprintf("Kendall's $\\tau_0$ & %s & %s & \\\\", f2(pv$frank$dep[["tau0"]], 2), f2(pv$gaussian$dep[["tau0"]], 2)),
     sprintf("Kendall's $\\tau_1$ & %s & %s & \\\\", f2(pv$frank$dep[["tau1"]], 2), f2(pv$gaussian$dep[["tau1"]], 2)),
     sprintf("LATE, ICT at home p25$\\to$p75 & %s %s & & \\\\", f2(res$late$estimate[1]), se(res$late$se[1])),
     sprintf("Complier share & %s & & \\\\", f2(res$late$estimate[2]))), "tab_main.tex")

# Table 4: exclusion diagnostics and sensitivity to a direct effect delta
fs <- res$first_stage$coef; ex <- res$excl_check
diag <- c(sprintf("ICT at home in probit for $D$ & %s & %s\\\\", f2(fs[["Estimate"]]), se(fs[["Std. Error"]])),
          sprintf("ICT at home in earnings, $D=0$ (OLS) & %s & %s\\\\", f2(ex[1, 1]), se(ex[2, 1])),
          sprintf("ICT at home in earnings, $D=1$ (OLS) & %s & %s\\\\", f2(ex[1, 2]), se(ex[2, 2])))
wr(diag, "tab_excl.tex")
sens <- sens[order(sens$family, sens$delta), ]
srows <- sapply(seq_len(nrow(sens)), function(i) with(sens[i, ],
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s\\\\", fam_lab[family], f2(delta, 2),
          format(round(AIC, 1), nsmall = 1, big.mark = ","), f2(tau0, 2), f2(tau1, 2),
          f2(ATE), f2(ATT), f2(ATU), f2(OLS))))
fams <- unique(sens$family)
out <- c(); for (f in fams) { out <- c(out, srows[sens$family == f]); if (f != tail(fams, 1)) out <- c(out, "\\midrule") }
wr(out, "tab_sens.tex")

# Table 5: robustness
rl <- c(main = "Main specification (Frank)", no_occ_fe = "No occupation effects", numeracy = "Numeracy instead of literacy",
        gpt_exposure = "GPT exposure defines the sample", ict_top_half = "Adoption = top half of ICT use",
        t_copula = "Student-$t$ copula", no_exclusion = "No exclusion restriction")
fr <- sel$fits$frank
main <- list(n = fr$n, te = res$sel_te$frank, dep = res$sel_dep$frank)
rr <- c(list(main = main), rob)
wr(sapply(names(rl), function(k) { x <- rr[[k]]
  sprintf("%s & %s & %s & %s & %s & %s & %s\\\\", rl[k], format(x$n, big.mark = ","), f2(x$te$estimate[1]),
          f2(x$te$estimate[2]), f2(x$te$estimate[3]), f2(x$dep[["tau0"]], 2), f2(x$dep[["tau1"]], 2)) }), "tab_rob.tex")

# Table 6: heterogeneity
h <- res$het; gl <- c(women = "Women", men = "Men", age_u40 = "Under 40", age_40p = "40 and over",
  tertiary = "15+ years of schooling", non_tertiary = "Under 15 years", lit_top_half = "Literacy, top half",
  lit_bottom_half = "Literacy, bottom half", country_152 = "Chile", country_203 = "Czechia", country_250 = "France",
  country_376 = "Israel", country_428 = "Latvia", country_616 = "Poland")
wr(sapply(names(gl), function(g) { x <- h[h$group == g, ]
  sprintf("%s & %s & %s %s & %s %s & %s %s\\\\", gl[g], format(x$n[1], big.mark = ","), f2(x$estimate[1]), se(x$se[1]),
          f2(x$estimate[2]), se(x$se[2]), f2(x$estimate[3]), se(x$se[3])) }), "tab_het.tex")

# Table 7: three regimes
j3 <- res$J3
wr(unlist(lapply(j3, function(x) {
  r <- function(nm, t) sprintf("%s & %s & %s & %s & %s & %s\\\\", nm, f2(t$estimate[1]), f2(t$estimate[2]),
                               f2(t$estimate[3]), f2(t$estimate[4]), "")
  c(sprintf("\\multicolumn{6}{l}{\\emph{%s copula; AIC %s; $\\tau=(%s)$}}\\\\", fam_lab[x$family],
            format(round(x$AIC, 1), nsmall = 1, big.mark = ","), paste(f2(x$dep, 2), collapse = ", ")),
    r("Augmented vs human-only", x$te10), r("Dominant vs human-only", x$te20), r("Dominant vs augmented", x$te21))
})), "tab_j3.tex")

# Figure 2: MTE; Figure 3: propensity overlap; Figure 4: sensitivity
pdf(file.path(fd, "fig_mte.pdf"), width = 6.5, height = 4)
m <- res$mte; par(mar = c(4.2, 4.2, 1, 1))
plot(m$u, m$mte, type = "n", ylim = range(c(m$lower, m$upper, 0)), xlab = "Unobserved resistance to adoption, u",
     ylab = "MTE (log points)")
polygon(c(m$u, rev(m$u)), c(m$lower, rev(m$upper)), col = "grey85", border = NA)
lines(m$u, m$mte, lwd = 2); abline(h = res$ols[["Estimate"]], lty = 2)
legend("topright", c("MTE, Frank copula (PV1)", "95% pointwise band", "OLS"), lwd = c(2, 8, 1), lty = c(1, 1, 2),
       col = c("black", "grey85", "black"), bty = "n", cex = 0.85)
dev.off()
pdf(file.path(fd, "fig_overlap.pdf"), width = 6.5, height = 4)
ps <- res$pscore; par(mar = c(4.2, 4.2, 1, 1)); br <- seq(0, 1, by = 0.025)
h0 <- hist(ps$p[ps$d == 0], breaks = br, plot = FALSE); h1 <- hist(ps$p[ps$d == 1], breaks = br, plot = FALSE)
plot(h0$mids, h0$density, type = "s", ylim = c(0, max(h0$density, h1$density)), xlab = "Estimated propensity score",
     ylab = "Density", lwd = 2, lty = 2)
lines(h1$mids, h1$density, type = "s", lwd = 2)
legend("topright", c("Non-adopters (D = 0)", "Adopters (D = 1)"), lty = c(2, 1), lwd = 2, bty = "n", cex = 0.85)
dev.off()
pdf(file.path(fd, "fig_sens.pdf"), width = 6.5, height = 4)
par(mar = c(4.2, 4.2, 1, 1)); cols <- c(frank = "black", gaussian = "grey45", clayton = "grey70")
yl <- range(c(sens$ATE, sens$OLS, 0))
plot(NA, xlim = range(sens$delta), ylim = yl, xlab = expression("Assumed direct effect of ICT at home on log earnings, " * delta),
     ylab = "ATE (log points)")
for (f in fams) { s <- sens[sens$family == f, ]; lines(s$delta, s$ATE, type = "b", pch = 19, col = cols[f], lwd = 2) }
s <- sens[sens$family == fams[1], ]; lines(s$delta, s$OLS, lty = 2, lwd = 2)
legend("bottomleft", c(paste(fam_lab[fams], "copula"), "OLS"), col = c(cols[fams], "black"), lty = c(rep(1, length(fams)), 2),
       pch = c(rep(19, length(fams)), NA), lwd = 2, bty = "n", cex = 0.85)
dev.off()

# numbers quoted in the text
key <- list(n = res$sample$n, n_occ = length(res$sample$occ2), share_d = res$sample$share_d,
            ols = res$ols[["Estimate"]], raw = res$desc$mean_d1[1] - res$desc$mean_d0[1],
            overlap_out = res$overlap$share_out, late = res$late$estimate[1], compl = res$late$estimate[2],
            J3_n = res$J3_n, J3_shares = as.numeric(res$J3_shares))
dput(key, file.path(td, "key_numbers.R"))
cat("tables and figures written\n")

# typeset minus signs; tidy the three-regime table
for (f in list.files(td, "^tab_.*\\.tex$", full.names = TRUE)) {
  s <- readLines(f)
  s <- gsub("(^|[^[:alnum:]$(])-([0-9])", "\\1$-$\\2", s, perl = TRUE)
  s <- gsub("\\(-([0-9])", "($-$\\1", s)
  if (basename(f) == "tab_j3.tex") {
    hd <- grepl("^\\\\multicolumn", s)
    s[hd] <- gsub("multicolumn\\{6\\}", "multicolumn{5}", gsub("\\$-\\$", "-", s[hd]))
    s[!hd] <- sub(" & \\\\\\\\$", "\\\\\\\\", s[!hd])
  }
  writeLines(s, f)
}
