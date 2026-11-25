# RoyAI

Generalized Roy models of human-AI complementarity and endogenous selection, in R.

| Module | What it does |
|---|---|
| `roy_core` | `roy()`: Model A (binary AI adoption; FIML or Heckman two-step; generalized, extended and classical Roy) and Model B (human-only / AI-augmented / AI-dominant; MNL + Dubin-McFadden or Lee correction) |
| `roy_select` | probit and multinomial-logit selection stages, correction terms, propensities |
| `roy_sim` | data-generating processes for Models A, B and the classical Roy model; `roy_mc()` Monte Carlo runner |
| `roy_diagnostics` | convergence, selection-on-unobservables test, exclusion-restriction strength, overlap, boundary checks; `roy_lrtest()`, `roy_ic()` |
| effects | `roy_potential()`, `roy_effects()` (ATE/ATT/ATU and pairwise regime effects) |
| `roy_copula`, `roy_mte` | Model C: copula Roy model (Gaussian, t, Clayton, Gumbel, Frank), multistart with regularity checks; ATE/ATT/ATU/MTE/LATE |
| `roy_continuous` | Model D: continuous AI-use intensity with a mass of non-users; copula ML (`roy_continuous()`), dose-response and effects at a given intensity (`roy_continuous_effects()`), simulator |
| `roy_dynamic` | Model E: dynamic Roy selection with Kalman learning about complementarity and learning by doing; backward-induction solver, simulator, NFXP and CCP estimators (`roy_dynamic()`), dynamic effects |
| data | `piaac_onet_prep()`, `onet_to_isco2()`, `piaac_varmap()` (stub, see VERIFY notes) |
| plug-ins | `roy_register_dependence()`, `roy_register_method()`; see `PLUGIN_CONTRACT.md` |

The interactive Shiny platform built on this package is at [brassbe1982/RoyAIapp](https://github.com/brassbe1982/RoyAIapp).

## Install (RStudio / Positron)

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("brassbe1982/RoyAI")

# or from a local copy of the source folder
install.packages("RoyAI", repos = NULL, type = "source")
# or, during development
devtools::load_all("RoyAI"); devtools::test("RoyAI")
```

## Quick start

```r
library(RoyAI)
dat <- roy_sim_A(3000)
fit <- roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = dat)
summary(fit); roy_effects(fit); roy_diagnostics(fit)
```

Model D and Model E:

```r
d <- roy_continuous_sim(2000, tau = 0.4)
fitD <- roy_continuous(y ~ x1, a ~ x1 + z1 + z2, data = d, family = "frank")
roy_continuous_effects(fitD)

p <- roy_dynamic_sim(1000, T = 10, psi = 0.05)
fitE <- roy_dynamic(p, cost = ~ x, method = "nfxp")   # or "ccp"
roy_dynamic_effects(fitE)
```

Monte Carlo checks of both models: `inst/montecarlo-DE/` (README there).

See `vignette("RoyAI-intro")`.

## Reference publication

This package is part of the RoyAI project, which builds on the following peer-reviewed publication:

> Niankara, I., Aljallad, Y., Soussi, O., Maghrabi, K., Alkhyeli, M. and Alghool, S. (2025). *The Impact of Artificial Intelligence on Fresh Graduate Labor Demand: A Theoretical and Empirical Analysis with Industry and Firm Heterogeneity.* In 2025 Global Congress on Emerging Technologies (GCET), Lyon, France, 2-5 December 2025. IEEE. [doi:10.1109/GCET68529.2025.11450624](https://doi.org/10.1109/GCET68529.2025.11450624). Open version: [BRASS Digital Lab](https://brassbe1982-bdlweb.share.connect.posit.cloud/research/ai-graduate-labor-demand-agentic.html).

`citation("RoyAI")` prints this reference and the software citation in BibTeX.

## Author and funding

Ibrahim Niankara (author, maintainer; Ibrahim.niankara@aau.ac.ae). Funded by BRASS Digital Lab, an LLC registered in Burkina Faso (RCCM BF-OUA-01-2026-B13-02466; https://brassbe1982-bdlweb.share.connect.posit.cloud/).
