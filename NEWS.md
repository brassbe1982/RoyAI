# RoyAI 0.3.0

First release prepared for CRAN. No change to estimates, likelihoods or
defaults relative to 0.2.0.

* Every exported function has a runnable example and a documented return
  value. Examples that need the internet (`download_piaac_puf()`) or a file
  from the U.S. Bureau of Labor Statistics (`read_soc_2010_2018()`) are in
  `\dontrun{}`; `roy_copula_select()` and `roy_mc()` are in `\donttest{}`.
* Documentation is written in roxygen2 markdown, so code and links render
  in the help pages.
* `DESCRIPTION` cites the methods with DOIs and gains `URL` and `BugReports`
  (<https://github.com/brassbe1982/RoyAI>).
* The ten slowest test blocks (more than 3 seconds each) call
  `testthat::skip_on_cran()`; they run when `NOT_CRAN` is `"true"`.
* GitHub Actions checks the package on Windows, macOS and Linux for R
  release, devel and oldrel.
* Analysis scripts in `inst/scripts/` read paths relative to a project root
  (`Sys.getenv("ROYAI_ROOT", ".")`) instead of fixed folders.

# RoyAI 0.2.0

* Model D: `roy_continuous()`, a copula Roy model for continuous AI-use
  intensity with a mass of non-users, estimated by maximum likelihood with
  the multistart and regularity guard of Model C;
  `roy_continuous_effects()` gives the dose-response curve and the TT, ATE
  and TUT of more intensity; `roy_continuous_sim()` and
  `roy_continuous_truth()` simulate and compute true effects.
* Model E: `roy_dynamic()`, a dynamic Roy model with Kalman learning about
  complementarity with AI and learning by doing, solved by backward
  induction (`roy_dynamic_solve()`) and estimated by nested fixed point or
  conditional choice probabilities; `roy_dynamic_effects()` and
  `roy_dynamic_sim()`.
* Multistart maximum likelihood keeps the best *regular* maximum (positive
  definite Hessian, small score, dependence away from the boundary) and
  warns when it skips a higher irregular boundary maximum.
* Monte Carlo studies of Models D and E in `inst/montecarlo-DE/`.

# RoyAI 0.1.0

* Model A (binary Roy model; full ML and Heckman two-step; generalized,
  extended and classical variants), Model B (three production arrangements,
  multinomial-logit selection with Dubin-McFadden or Lee correction) and
  Model C (copula generalized Roy model with Gaussian, Student-t, Clayton,
  Gumbel and Frank families; ATE, ATT, ATU, MTE and LATE).
* Pluggable dependence families and estimation methods
  (`roy_register_dependence()`, `roy_register_method()`).
* PIAAC Cycle 2 and O*NET data tools.
