#!/usr/bin/env Rscript
# Run one independent R reproduction cell.
# Args: scenario sample_size n_sim seed cores out_json
suppressPackageStartupMessages({
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("Usage: scenario sample_size n_sim seed cores out_json")
scenario <- args[[1]]
ss <- as.integer(args[[2]])
n_sim <- as.integer(args[[3]])
seed <- as.integer(args[[4]])
num.core <- max(1L, as.integer(args[[5]]))
out_json <- args[[6]]

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
orig_dir <- normalizePath(file.path(dirname(this_file), "original_r"))
data_sim <- file.path(orig_dir, "data_sim.R")
test_sim <- file.path(orig_dir, "test_sim.R")
source(data_sim, local = FALSE)
source(test_sim, local = FALSE)

make_fu <- function(eff, eff_efu, et_trt, et_pbo, efu_rate) {
  list(eff = eff, eff.efu = eff_efu, et.rate.trt = et_trt, et.rate.pbo = et_pbo, efu.rate = efu_rate)
}
sc <- switch(
  scenario,
  Ex1 = list(fu1 = make_fu(4, 2, 0.2, 0.1, 0.5), fu2 = make_fu(6, 3, 0.2, 0.1, 0.5)),
  Ex2 = list(fu1 = make_fu(4, 2, 0.2, 0.1, 0), fu2 = make_fu(6, 3, 0.2, 0.1, 0)),
  Ex3 = list(fu1 = make_fu(4, 0, 0.2, 0.1, 0.5), fu2 = make_fu(6, 0, 0.2, 0.1, 0.5)),
  Null = list(fu1 = make_fu(0, 0, 0.2, 0.1, 0.5), fu2 = make_fu(0, 0, 0.2, 0.1, 0.5)),
  stop("Unknown scenario")
)

wilson_ci <- function(p, n, z = 1.96) {
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den
  c(max(0, center - half), min(1, center + half))
}

cl <- parallel::makeCluster(num.core)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

# Source original functions on every worker (avoids fragile clusterExport of closures).
parallel::clusterExport(cl, c("data_sim", "test_sim"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(mmrm)
    library(dplyr)
    library(MASS)
  })
  source(data_sim, local = FALSE)
  source(test_sim, local = FALSE)
  NULL
}))

fu1 <- sc$fu1
fu2 <- sc$fu2
parallel::clusterExport(cl, c("ss", "fu1", "fu2"), envir = environment())
parallel::clusterSetRNGStream(cl, seed)

p_vals <- as.numeric(unlist(parallel::parLapply(cl, seq_len(n_sim), function(i) {
  d <- MMRM.sim.1(i, ss, 0.5, 38, 10.5, 0.5, fu1, fu2, NULL)
  MMRM.test.1(d)
})))

ok <- is.finite(p_vals)
m <- sum(ok)
power <- if (m > 0) mean(p_vals[ok] < 0.05) else NA_real_
mcse <- if (m > 0) sqrt(power * (1 - power) / m) else NA_real_
ci <- if (m > 0) wilson_ci(power, m) else c(NA_real_, NA_real_)

jsonlite::write_json(list(
  scenario = scenario,
  sample_size = ss,
  n_sim = n_sim,
  seed = seed,
  power_pct = 100 * power,
  mcse_power_pct = 100 * mcse,
  ci95_power_lo_pct = 100 * ci[1],
  ci95_power_hi_pct = 100 * ci[2],
  R_version = R.version.string,
  mmrm_version = as.character(utils::packageVersion("mmrm")),
  valid_p = m
), out_json, auto_unbox = TRUE, pretty = TRUE)
