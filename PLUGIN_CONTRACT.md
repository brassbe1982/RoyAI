# RoyAI plug-in contract (core <-> roy_copula / roy_mte)

Conventions (all modules):
- Selection: `D = 1(Z'gamma + V > 0)`, `V ~ N(0,1)` (probit margin).
- Outcomes: `Y_s = X'beta_s + U_s`, `U_s ~ N(0, sigma_s^2)` margins.
- Regimes are coded 0-based: Model A 0 = no AI, 1 = AI adopter; Model B 0 = human-only, 1 = AI-augmented, 2 = AI-dominant.

## 1. Dependence families (Model A likelihood)

```r
fam <- roy_dependence(name, h, npar, linkinv, link, start,
                      tau = NULL, cond_mean = NULL, par_names = ...)
roy_register_dependence(fam)
```
- `h(v, u, theta)` = P(V <= v | U = u), uniforms `u = F_s(U_s)`, `v = Phi(V)`, vectorised in `u`, `v`.
- `linkinv(eta)` maps unconstrained reals to `theta`; `link` is its inverse.
- Optional `tau(theta)` (Kendall's tau), `cond_mean(c, theta, select)` = E[e_s | D = select] for standardized e_s and index c (otherwise integrated by Gauss-Hermite).
- Simulation (`roy_sim_A(dependence = name)`) inverts `h` in its first argument, which assumes an exchangeable copula; rotated families should ship their own simulator.
- Natural parameter names in a fit: `paste0(par_names, 0)` for (U0, V) and `paste0(par_names, 1)` for (U1, V).

Registration at load time: if the package namespace contains a function `roy_copula_register()`, `.onLoad()` (R/zzz.R) calls it.

## 2. Estimation methods

```r
roy_register_method(model = "A" | "B", method = "name", fitter = function(spec, ...) ...)
```
`roy(..., method = "name")` then calls `fitter(spec, dependence =, se =, B =, control =, ...)`.
The fitter must return an object inheriting from `"roy_fit"` with at least:
`spec`, `model`, `method`, `type`, `family`, `par` (named natural-scale vector), `vcov`,
`coef` (list: `selection`, `outcome` = list per regime, plus `sigma`, `dependence`, `cov_sd` for Model A),
`relist` (function(par) -> coef list), `loglik`, `npar`, `converged`, `se_type`.

## 3. What roy_mte can rely on

- `roy_propensity(fit, newdata)`: n x J matrix of regime probabilities.
- `roy_potential(fit, regime, given = NULL | k, newdata, coef_list)`: E[Y_s | X] or E[Y_s | X, Z, S = k].
- `fit$coef$selection`, `fit$coef$outcome`, `fit$coef$sigma`, `fit$coef$dependence`, `fit$family$h`.
- `roy_effects(fit)`: ATE / ATT / ATU (Model A) and pairwise effects (Model B), delta-method or bootstrap SEs. MTE, LATE, PRTE belong to roy_mte.
