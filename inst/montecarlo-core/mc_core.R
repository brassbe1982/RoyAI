# Monte Carlo sanity study for the RoyAI core estimators.
# Run: Rscript inst/montecarlo-core/mc_core.R  (about 10-15 minutes)
library(RoyAI)
R <- 200
out <- list()
out$A_ml <- roy_mc(function() roy_sim_A(2000),
                   function(d) roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = d), R = R, seed = 101)
out$A_twostep <- roy_mc(function() roy_sim_A(2000),
                        function(d) roy(y ~ x1 + x2, d ~ x1 + x2 + z, data = d, method = "twostep"),
                        R = R, seed = 101)
out$classical <- roy_mc(function() roy_sim_classical(2000),
                        function(d) roy(y ~ x1 + x2, d ~ x1 + x2, data = d, type = "classical"),
                        R = R, seed = 102)
out$B_dmf <- roy_mc(function() roy_sim_B(4000),
                    function(d) roy(y ~ x1 + x2, s ~ x1 + x2 + z1 + z2, data = d), R = R, seed = 103)
for (k in names(out)) { cat("\n#####", k, "\n"); print(out[[k]]) }
saveRDS(lapply(out, `[[`, "summary"), "mc_core_results.rds")
