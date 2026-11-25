# Monte Carlo checks for Model D and Model E

Run from this folder (needs RoyAI >= 0.2.0 installed):

```
Rscript mc_modelD.R 200 4   # replications, cores
Rscript mc_modelE.R 100 4
```

Outputs go to `results/`: `mc_modelD_summary.csv`, `mc_modelE_summary.csv` (bias, RMSE,
95% coverage, share of regular maxima, mean seconds per fit), the raw draws
(`*_raw.rds`, excluded from the package tarball) and the run logs.

## Model D (`roy_continuous()`)

DGP `roy_continuous_sim()`: mu(x, a) = 1 + 0.5 x1 + 0.8 a - 0.4 a^2,
log sigma(a) = log 0.6 + 0.5 a, selection Z = (1, x1, z1, z2), gamma = (0.6, 0.3, 0.8, 0.6)
(about 27% non-users), Lambda(t) = 1 - exp(-t), Kendall tau = 0.4.
Cells: Gaussian and Frank copulas x n in {1000, 3000}, 200 replications each.
Estimators: copula ML with the true family (`ml`), Gaussian ML on Frank data
(`ml_gauss`, misspecified copula), two-step control function (`cf`, Gaussian and
homoskedastic) and naive OLS (`ols`). Estimands: dose coefficients, alpha_a, tau,
the dose-response E[Y(a)] at a = 0, 0.5, 0.95, and at a = 0.2, 0.5, 0.8 the effects
TT (0.1 more intensity for workers at a), ATE (0.1 more for everyone) and TUT
(non-users moved to a). True effects: `roy_continuous_truth()` over 200,000
covariate draws.

## Model E (`roy_dynamic()`)

DGP `roy_dynamic_sim(N = 1000, T = 10)`: mu0 = 1, sigma0 = 0.3, kappa_1 ~ N(0, 0.25^2),
sigma_e = 0.3, cost = 0.1 + 0.1 x, lambda = 0.1, beta = 0.95 (fixed), psi in {0, 0.1}.
100 replications per psi. Estimators: CCP (Hotz-Miller with the human-only path) and
NFXP (full-information ML). Static comparison on the period-10 cross-section: the naive
earnings gap and the Gaussian Roy model (`roy()`, x as the excluded cost shifter)
against the true ATT and ATE of period 10 (computed from the simulated kappa).
CCP stage-2 standard errors are conditional on stage 1, so CCP coverage for the
cost and lambda parameters is not a valid test of the estimator.
