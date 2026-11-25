# Model C Monte Carlo (Paper 1, Section 7)

Run: `Rscript mc_copula.R <main_n2000|main_n500|misspec|weakiv|j3|main_large> [cores]` or `./run_all.sh`.
Each cell is checkpointed to `results/raw/<cell>.rds`; reruns skip finished cells.
Summarise: `Rscript mc_tables.R` -> `results/tables/` (CSV + LaTeX) and `results/figures/`.

DGP: X=(1,x1), Z=(1,x1,z1,z2), x1,z1~N(0,1), z2~Bern(0.5); beta0=(0.5,0.3), beta1=(1,0.5),
gamma=(0.2,0.3,0.8,0.6), sigma=(0.8,1), (tau0,tau1)=(-tau,+tau), t df=5. LATE: z2 0->1.
Weak IV: gamma=(0.2,0.3,0.15,0.10). J=3 ordered: beta2=(1.3,0.8), sigma2=1.2, kappa2=1.2,
taus=(-0.6,0.3,0.6), n=3000. True effects: analytic, averaged over 1e6 covariate draws.

| Output | Design | R |
|---|---|---|
| table1_fiml_parameters | 5 copulas x tau {0,.3,.6} x n {500,2000,10000}, FIML | 500 (t: 250); n=10000: 100 (t: 50) |
| table2_treatment_effects | same cells; OLS vs normal two-step vs copula FIML | same |
| table3_misspecification | Gaussian FIML on every DGP vs correct copula, n=2000 | 100 |
| table4_copula_selection | AIC/BIC choice among all 5 copulas, n=2000 | 100 |
| table_weak_instrument | Gaussian/Clayton/Gumbel, tau=.6, n=2000 | 300 |
| table5_three_regimes | J=3 ordered, tau=.6, n=3000 | 200 (t: 100) |
| figure1_mte | true vs estimated MTE, Gaussian and Gumbel DGP, tau=.6, n=2000 | 500 / 100 |

Two-step coverage is not reported (no generated-regressor correction). 0 failed fits in 63 cells.
