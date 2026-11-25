# =============================================================================
# Summarise the Model C Monte Carlo: tables (CSV + LaTeX) and figures for Paper 1
# Usage: Rscript mc_tables.R      (reads results/raw/*.rds, writes results/)
# Works on partial results: cells that have not finished are skipped.
# =============================================================================
MC_LOAD_ONLY <- TRUE
MC_HERE <- local({
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
source(file.path(MC_HERE, "mc_copula.R"))
res_dir <- file.path(MC_HERE, "results")
tab_dir <- file.path(res_dir, "tables"); dir.create(tab_dir, FALSE, TRUE)
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, FALSE, TRUE)

cells <- list.files(file.path(res_dir, "raw"), pattern = "\\.rds$", full.names = TRUE)
load_cell <- function(f) readRDS(f)
parse_cell <- function(name) {
  p <- strsplit(sub("\\.rds$", "", basename(name)), "_")[[1]]
  list(design = p[1], family = p[2], tau = as.numeric(sub("tau", "", p[3])),
       n = as.numeric(sub("n", "", p[4])))
}
FAMLAB <- c(gaussian = "Gaussian", t = "Student-t", clayton = "Clayton",
            gumbel = "Gumbel", frank = "Frank")

summ <- function(est, se, truth) {
  ok <- is.finite(est)
  e <- est[ok]; s <- se[ok]
  c(true = truth, mean = mean(e), bias = mean(e) - truth,
    sd = stats::sd(e), rmse = sqrt(mean((e - truth)^2)),
    mean_se = if (all(is.na(s))) NA else mean(s, na.rm = TRUE),
    coverage = if (all(is.na(s))) NA else mean(abs(e - truth) <= 1.96 * s, na.rm = TRUE),
    R = sum(ok))
}

fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
write_tex <- function(df, file, caption, label, digits = 3) {
  num <- vapply(df, is.numeric, TRUE)
  body <- apply(df, 1, function(r) paste(r, collapse = " & "))
  df[num] <- lapply(df[num], fmt, d = digits)
  body <- apply(df, 1, function(r) paste(r, collapse = " & "))
  tex <- c("\\begin{table}[htbp]\\centering\\small",
           sprintf("\\caption{%s}\\label{%s}", caption, label),
           sprintf("\\begin{tabular}{%s}", paste(rep("l", ncol(df)), collapse = "")),
           "\\hline", paste(gsub("_", "\\\\_", names(df)), collapse = " & "), "\\\\ \\hline",
           paste0(body, " \\\\"), "\\hline\\end{tabular}\\end{table}")
  writeLines(tex, file)
}

main_files <- cells[grepl("^main_", basename(cells))]
rows_par <- list(); rows_te <- list(); mte_store <- list()
for (f in main_files) {
  info <- parse_cell(f); obj <- load_cell(f)
  reps <- Filter(Negate(is.null), obj$reps)
  tr <- get_truth(info$family, info$tau)
  pn <- names(reps[[1]]$fiml[[info$family]]$par)
  truth_par <- stats::setNames(tr$par, pn)
  # structural parameters, FIML (correct family) and two-step
  for (p in pn) {
    fe <- sapply(reps, function(r) { z <- r$fiml[[info$family]]; if (is.null(z)) NA else z$par[[p]] })
    fs <- sapply(reps, function(r) { z <- r$fiml[[info$family]]; if (is.null(z)) NA else z$se[[p]] })
    rows_par[[length(rows_par) + 1]] <- data.frame(info, estimator = "FIML", parameter = p,
                                                   t(summ(fe, fs, truth_par[[p]])))
    if (p %in% names(reps[[1]]$twostep$par)) {
      te2 <- sapply(reps, function(r) r$twostep$par[[p]])
      rows_par[[length(rows_par) + 1]] <- data.frame(info, estimator = "Two-step", parameter = p,
                                                     t(summ(te2, NA * te2, truth_par[[p]])))
    }
    if (grepl("^beta", p)) {
      k <- match(p, pn[grepl("^beta", pn)])
      ob <- sapply(reps, function(r) r$ols$beta[[k]])
      rows_par[[length(rows_par) + 1]] <- data.frame(info, estimator = "OLS", parameter = p,
                                                     t(summ(ob, NA * ob, truth_par[[p]])))
    }
  }
  for (q in c("ATE", "ATT", "ATU", "LATE")) {
    est <- list(
      OLS = cbind(sapply(reps, function(r) r$ols$te[[q]]), sapply(reps, function(r) r$ols$te_se[[match(q, c("ATE","ATT","ATU","LATE"))]])),
      `Two-step` = cbind(sapply(reps, function(r) r$twostep$te[[q]]), NA),
      FIML = cbind(sapply(reps, function(r) { z <- r$fiml[[info$family]]; if (is.null(z)) NA else z$te[[q]] }),
                   sapply(reps, function(r) { z <- r$fiml[[info$family]]; if (is.null(z)) NA else z$te_se[[match(q, c("ATE","ATT","ATU","LATE"))]] })))
    for (e in names(est))
      rows_te[[length(rows_te) + 1]] <- data.frame(info, estimator = e, parameter = q,
                                                   t(summ(est[[e]][, 1], est[[e]][, 2], tr$te[[q]])))
  }
  m <- do.call(rbind, lapply(reps, function(r) r$fiml[[info$family]]$mte))
  if (!is.null(m)) mte_store[[basename(f)]] <- list(info = info, truth = tr$mte, est = m)
}

if (length(rows_par)) {
  par_tab <- do.call(rbind, rows_par); te_tab <- do.call(rbind, rows_te)
  write.csv(par_tab, file.path(tab_dir, "mc_parameters_all.csv"), row.names = FALSE)
  write.csv(te_tab, file.path(tab_dir, "mc_treatment_effects_all.csv"), row.names = FALSE)
  # Table 1: FIML recovery of key parameters
  key <- c("beta1:x1", "beta0:x1", "tau1", "tau0")
  t1 <- subset(par_tab, estimator == "FIML" & parameter %in% key)
  t1 <- t1[order(match(t1$family, names(FAMLAB)), t1$tau, t1$n, match(t1$parameter, key)), ]
  t1o <- data.frame(Copula = FAMLAB[t1$family], tau = t1$tau, n = t1$n, Parameter = t1$parameter,
                    True = t1$true, Bias = t1$bias, RMSE = t1$rmse, Coverage = t1$coverage, R = t1$R)
  write.csv(t1o, file.path(tab_dir, "table1_fiml_parameters.csv"), row.names = FALSE)
  write_tex(t1o, file.path(tab_dir, "table1_fiml_parameters.tex"),
            "Copula FIML recovery of outcome slopes and Kendall's tau (bias, RMSE, 95\\% Wald coverage).",
            "tab:mc-params")
  # Table 2: treatment effects by estimator
  t2 <- te_tab[order(match(te_tab$family, names(FAMLAB)), te_tab$tau, te_tab$n,
                     match(te_tab$parameter, c("ATE", "ATT", "ATU", "LATE")),
                     match(te_tab$estimator, c("OLS", "Two-step", "FIML"))), ]
  t2o <- data.frame(Copula = FAMLAB[t2$family], tau = t2$tau, n = t2$n, Effect = t2$parameter,
                    Estimator = t2$estimator, True = t2$true, Bias = t2$bias, RMSE = t2$rmse,
                    Coverage = t2$coverage)
  write.csv(t2o, file.path(tab_dir, "table2_treatment_effects.csv"), row.names = FALSE)
  write_tex(subset(t2o, n == 2000 & tau > 0), file.path(tab_dir, "table2_treatment_effects_n2000.tex"),
            "Treatment parameters by estimator, $n = 2000$. Two-step coverage omitted (no generated-regressor correction).",
            "tab:mc-te")
}

# ---- Table 3 and 4: misspecification and copula selection -------------------
mis_files <- cells[grepl("^misspec_", basename(cells))]
rows3 <- list(); rows4 <- list(); fig_mis <- list()
for (f in mis_files) {
  info <- parse_cell(f); obj <- load_cell(f)
  reps <- Filter(Negate(is.null), obj$reps)
  tr <- get_truth(info$family, info$tau)
  for (fitfam in unique(c(info$family, "gaussian"))) {
    for (q in c("ATE", "ATT", "ATU", "LATE")) {
      e <- sapply(reps, function(r) { z <- r$fiml[[fitfam]]; if (is.null(z$te)) NA else z$te[[q]] })
      s <- sapply(reps, function(r) { z <- r$fiml[[fitfam]]; if (is.null(z$te)) NA else z$te_se[[match(q, c("ATE","ATT","ATU","LATE"))]] })
      rows3[[length(rows3) + 1]] <- data.frame(DGP = FAMLAB[info$family], tau = info$tau,
                                               Fitted = FAMLAB[fitfam], Effect = q,
                                               t(summ(e, s, tr$te[[q]])))
    }
    e <- sapply(reps, function(r) { z <- r$fiml[[fitfam]]; if (is.null(z)) NA else z$par[["tau1"]] })
    rows3[[length(rows3) + 1]] <- data.frame(DGP = FAMLAB[info$family], tau = info$tau,
                                             Fitted = FAMLAB[fitfam], Effect = "tau1",
                                             t(summ(e, NA * e, info$tau)))
  }
  ic <- lapply(reps, function(r) {
    ll <- sapply(FAMS, function(g) if (is.null(r$fiml[[g]])) NA else r$fiml[[g]]$loglik)
    k  <- sapply(FAMS, function(g) if (is.null(r$fiml[[g]])) NA else r$fiml[[g]]$npar)
    c(aic = names(which.min(-2 * ll + 2 * k)), bic = names(which.min(-2 * ll + log(2000) * k)))
  })
  ic <- do.call(rbind, ic)
  for (crit in c("aic", "bic")) {
    tabc <- table(factor(ic[, crit], levels = FAMS)) / nrow(ic)
    rows4[[length(rows4) + 1]] <- data.frame(DGP = FAMLAB[info$family], tau = info$tau,
                                             Criterion = toupper(crit),
                                             t(setNames(as.numeric(tabc), FAMLAB[FAMS])),
                                             Hit = as.numeric(tabc[info$family]), R = nrow(ic),
                                             check.names = FALSE)
  }
  if (info$family == "gumbel" && info$tau == 0.6) {
    fig_mis <- list(truth = tr$mte,
                    gumbel = do.call(rbind, lapply(reps, function(r) r$fiml$gumbel$mte)),
                    gaussian = do.call(rbind, lapply(reps, function(r) r$fiml$gaussian$mte)))
  }
}
if (length(rows3)) {
  t3 <- do.call(rbind, rows3)
  t3o <- data.frame(t3[, c("DGP", "tau", "Fitted", "Effect")], True = t3$true, Bias = t3$bias,
                    RMSE = t3$rmse, Coverage = t3$coverage, R = t3$R)
  write.csv(t3o, file.path(tab_dir, "table3_misspecification.csv"), row.names = FALSE)
  write_tex(subset(t3o, tau == 0.6 & DGP != "Gaussian"), file.path(tab_dir, "table3_misspecification.tex"),
            "Gaussian copula fitted to non-Gaussian data versus the correct copula ($\\tau = 0.6$, $n = 2000$).",
            "tab:mc-misspec")
  t4 <- do.call(rbind, rows4)
  write.csv(t4, file.path(tab_dir, "table4_copula_selection.csv"), row.names = FALSE)
  write_tex(t4, file.path(tab_dir, "table4_copula_selection.tex"),
            "Share of replications in which AIC/BIC selects each copula ($n = 2000$). Hit = correct family.",
            "tab:mc-select", digits = 2)
}

# ---- weak instrument -----------------------------------------------------
wk_files <- cells[grepl("^weakiv_", basename(cells))]
rows_w <- list()
for (f in wk_files) {
  info <- parse_cell(f); reps <- Filter(Negate(is.null), load_cell(f)$reps)
  tr <- get_truth(info$family, info$tau, gamma = DGP$gamma_weak)
  for (q in c("ATE", "ATT", "ATU", "LATE")) for (e in c("OLS", "Two-step", "FIML")) {
    est <- switch(e, OLS = sapply(reps, function(r) r$ols$te[[q]]),
                  `Two-step` = sapply(reps, function(r) r$twostep$te[[q]]),
                  FIML = sapply(reps, function(r) r$fiml[[info$family]]$te[[q]]))
    se <- switch(e, OLS = sapply(reps, function(r) r$ols$te_se[[match(q, c("ATE","ATT","ATU","LATE"))]]), `Two-step` = NA * est,
                 FIML = sapply(reps, function(r) r$fiml[[info$family]]$te_se[[match(q, c("ATE","ATT","ATU","LATE"))]]))
    rows_w[[length(rows_w) + 1]] <- data.frame(Copula = FAMLAB[info$family], Effect = q,
                                               Estimator = e, t(summ(est, se, tr$te[[q]])))
  }
}
if (length(rows_w)) {
  tw <- do.call(rbind, rows_w)
  write.csv(tw, file.path(tab_dir, "table_weak_instrument.csv"), row.names = FALSE)
  write_tex(tw[, c("Copula", "Effect", "Estimator", "true", "bias", "rmse", "coverage")],
            file.path(tab_dir, "table_weak_instrument.tex"),
            "Weak excluded instruments ($\\gamma_{z1} = 0.15$, $\\gamma_{z2} = 0.10$), $\\tau = 0.6$, $n = 2000$.",
            "tab:mc-weakiv")
}

# ---- Table 5: J = 3 ordered model ----------------------------------------------
j3_files <- cells[grepl("^j3_", basename(cells))]
rows5 <- list()
for (f in j3_files) {
  info <- parse_cell(f); reps <- Filter(Negate(is.null), load_cell(f)$reps)
  tr <- get_truth(info$family, info$tau, J = 3)
  pn <- names(reps[[1]]$par); tp <- stats::setNames(tr$par, pn)
  for (p in intersect(c("beta1:x1", "beta2:x1", "kappa2", "tau0", "tau1", "tau2"), pn)) {
    rows5[[length(rows5) + 1]] <- data.frame(Copula = FAMLAB[info$family], Quantity = p,
      t(summ(sapply(reps, function(r) r$par[[p]]), sapply(reps, function(r) r$se[[p]]), tp[[p]])))
  }
  for (k in seq_along(tr$te)) {
    q <- names(tr$te)[k]
    rows5[[length(rows5) + 1]] <- data.frame(Copula = FAMLAB[info$family], Quantity = q,
      t(summ(sapply(reps, function(r) r$te[[k]]), sapply(reps, function(r) r$te_se[[k]]), tr$te[[q]])))
  }
}
if (length(rows5)) {
  t5 <- do.call(rbind, rows5)
  write.csv(t5, file.path(tab_dir, "table5_three_regimes.csv"), row.names = FALSE)
  write_tex(t5[, c("Copula", "Quantity", "true", "bias", "rmse", "coverage", "R")],
            file.path(tab_dir, "table5_three_regimes.tex"),
            "Three ordered regimes (human-only, AI-augmented, AI-dominant), $\\tau = 0.6$, $n = 3000$.",
            "tab:mc-j3")
}

# ---- Figure 1: MTE curves ----------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  curves <- list()
  add <- function(label, panel, m, truth) {
    if (is.null(m)) return()
    curves[[length(curves) + 1]] <<- data.frame(
      u = UGRID, panel = panel, series = label, mean = colMeans(m),
      lo = apply(m, 2, stats::quantile, 0.05), hi = apply(m, 2, stats::quantile, 0.95))
  }
  truth_rows <- list()
  for (fam in c("gaussian", "gumbel")) {
    key <- sprintf("main_%s_tau0.6_n2000.rds", fam)
    if (!is.null(mte_store[[key]])) {
      add(paste(FAMLAB[fam], "FIML (correct)"), paste(FAMLAB[fam], "DGP"), mte_store[[key]]$est, NULL)
      truth_rows[[fam]] <- data.frame(u = UGRID, panel = paste(FAMLAB[fam], "DGP"),
                                      series = "True MTE", mean = mte_store[[key]]$truth, lo = NA, hi = NA)
    }
  }
  if (length(fig_mis)) add("Gaussian FIML (misspecified)", "Gumbel DGP", fig_mis$gaussian, NULL)
  if (length(curves)) {
    dd <- do.call(rbind, c(curves, truth_rows))
    p <- ggplot(dd, aes(u, mean, colour = series)) +
      geom_ribbon(data = dd[!is.na(dd$lo), ], aes(ymin = lo, ymax = hi, fill = series),
                  alpha = 0.15, colour = NA) +
      geom_line(aes(linetype = series == "True MTE"), linewidth = 0.8) +
      scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid"), guide = "none") +
      facet_wrap(~panel) + guides(fill = "none") +
      labs(x = expression(u[S] ~ "(unobserved resistance to AI adoption)"),
           y = "MTE at x1 = 0", colour = NULL, fill = NULL,
           caption = "tau = 0.6, n = 2000; bands = 5-95% of estimates across replications") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(fig_dir, "figure1_mte.pdf"), p, width = 8, height = 4)
    ggsave(file.path(fig_dir, "figure1_mte.png"), p, width = 8, height = 4, dpi = 200)
    write.csv(dd, file.path(fig_dir, "figure1_mte_data.csv"), row.names = FALSE)
  }
}
cat("Tables written to", tab_dir, "\n")
