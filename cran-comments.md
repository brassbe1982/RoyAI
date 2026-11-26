## Submission

This is a first submission of RoyAI (version 0.3.0).

RoyAI estimates generalized Roy selection models of human-AI complementarity:
binary, multinomial and ordered copula Roy models, a copula model for
continuous treatment intensity with a mass of non-users, and a dynamic Roy
model with learning, with average, marginal and local average treatment
effects.

## Test environments

* Ubuntu 24.04, R 4.3.3: `R CMD check --as-cran` on the tarball built by
  `R CMD build`: Status OK.
* Windows 11 (maintainer's machine, R 4.6.2): `R CMD check --as-cran` of the
  development version: 0 errors, 0 warnings, 1 NOTE on "Lost braces" in three
  Rd files, all fixed in this version. [Replace with the run on 0.3.0.]
* GitHub Actions (`R CMD check --no-manual`, complete test suite with
  `NOT_CRAN = "true"`): Windows, macOS and Ubuntu with R release, devel and
  oldrel-1, <https://github.com/brassbe1982/RoyAI/actions>:
  all 9 jobs passed (0 errors, 0 warnings) on 2026-11-26.
* win-builder, R-devel and R-release: [TO BE RUN BY THE MAINTAINER:
  devtools::check_win_devel(), devtools::check_win_release(); results arrive by
  e-mail]
* macOS builder: [TO BE RUN BY THE MAINTAINER: devtools::check_mac_release()]

## R CMD check results

0 errors | 0 warnings | 1 note

* `checking CRAN incoming feasibility ... NOTE`
  `Maintainer: 'Ibrahim Niankara <Ibrahim.niankara@aau.ac.ae>'`
  `New submission`

  This is a first submission.

Words in DESCRIPTION that a spell checker may flag are proper names
(Bourguignon, Dubin, Fournier, Gurgand, Hotz, McFadden, Niankara, Roy) or
terms of the field (copula).

## Timing

The complete check takes about 2 minutes: examples 14 seconds (26 seconds
with `--run-donttest`), tests 23 seconds, vignette 35 seconds. The ten
slowest test blocks (about 100 seconds in total) call
`testthat::skip_on_cran()`; they run when `NOT_CRAN` is `"true"` and pass.

## Other

* No compiled code. Imports only stats, utils and numDeriv. Packages in
  Suggests are used conditionally.
* Examples, tests and the vignette use one core.
* Examples that need the internet (`download_piaac_puf()`) or a file
  distributed by the U.S. Bureau of Labor Statistics
  (`read_soc_2010_2018()`) are wrapped in `\dontrun{}`; every other exported
  function has an example that runs.
* Nothing is written outside `tempdir()`; the example of `read_piaac_puf()`
  writes a small file to `tempdir()` and deletes it.
* Source and issue tracker: <https://github.com/brassbe1982/RoyAI>.

## Downstream dependencies

None (new package).
