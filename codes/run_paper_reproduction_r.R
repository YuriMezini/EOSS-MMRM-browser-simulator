#!/usr/bin/env Rscript
# Independent R/mmrm reproduction using original data_sim.R + test_sim.R
# (PAP_AS03 settings). This is NOT paired with the JavaScript generator.
#
# Usage:
#   Rscript paper/codes/run_paper_reproduction_r.R [--n-sim 5000] [--smoke]
#
# Label everywhere:
#   "R/mmrm reproduction using the settings reported in PAP_AS03"

suppressPackageStartupMessages({
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
n_sim <- 5000L
seed <- 1234L
alpha <- 0.05
smoke <- FALSE
if ("--smoke" %in% args) {
  smoke <- TRUE
  n_sim <- 20L
}
if ("--n-sim" %in% args) {
  i <- match("--n-sim", args)
  n_sim <- as.integer(args[[i + 1L]])
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript so --file= is available")
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
root <- normalizePath(file.path(dirname(this_file), "../.."))

orig_dir <- file.path(root, "paper", "codes", "original_r")
out_dir <- file.path(root, "paper", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_id <- sprintf(
  "paper-repro-r_n%d_seed%d_%s",
  n_sim, seed,
  format(as.POSIXlt(Sys.time(), tz = "UTC"), "%Y%m%dT%H%M%SZ")
)
if (smoke) run_id <- paste0(run_id, "_smoke")

checkpoint <- file.path(out_dir, paste0("paper_reproduction_r_checkpoint_", run_id, ".csv"))
final_csv <- file.path(out_dir, if (smoke) "paper_reproduction_r_smoke.csv" else "paper_reproduction_r.csv")
manifest_path <- file.path(out_dir, if (smoke) "paper_reproduction_r_smoke_manifest.json" else "paper_reproduction_r_manifest.json")

detected <- parallel::detectCores()
if (is.na(detected) || detected < 1) detected <- 2L
num.core <- max(1L, as.integer(detected) - 1L)
assign("num.core", num.core, envir = .GlobalEnv)

# Source original generators/testers (use MMRM.sim.1 / MMRM.test.1 only;
# avoid MMRM.sim.n / MMRM.test.n which create a new cluster every call).
source(file.path(orig_dir, "data_sim.R"), local = FALSE)
source(file.path(orig_dir, "test_sim.R"), local = FALSE)

r_version <- R.version.string
mmrm_version <- as.character(utils::packageVersion("mmrm"))

sample_sizes <- c(40L, 60L, 80L, 100L, 120L, 140L, 160L, 180L, 200L, 220L, 240L)

make_fu <- function(eff, eff_efu, et_trt, et_pbo, efu_rate) {
  list(
    eff = eff,
    eff.efu = eff_efu,
    et.rate.trt = et_trt,
    et.rate.pbo = et_pbo,
    efu.rate = efu_rate
  )
}

scenarios <- list(
  list(
    id = "Ex1",
    fu1 = make_fu(4, 2, 0.2, 0.1, 0.5),
    fu2 = make_fu(6, 3, 0.2, 0.1, 0.5)
  ),
  list(
    id = "Ex2",
    fu1 = make_fu(4, 2, 0.2, 0.1, 0),
    fu2 = make_fu(6, 3, 0.2, 0.1, 0)
  ),
  list(
    id = "Ex3",
    fu1 = make_fu(4, 0, 0.2, 0.1, 0.5),
    fu2 = make_fu(6, 0, 0.2, 0.1, 0.5)
  ),
  list(
    id = "Null",
    fu1 = make_fu(0, 0, 0.2, 0.1, 0.5),
    fu2 = make_fu(0, 0, 0.2, 0.1, 0.5)
  )
)

header <- c(
  "scenario", "sample_size", "n_sim", "seed", "power_pct", "mcse_power_pct",
  "ci95_power_lo_pct", "ci95_power_hi_pct", "R_version", "mmrm_version", "run_id",
  "label"
)

read_done <- function(path) {
  if (!file.exists(path)) return(character())
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  paste(d$scenario, d$sample_size, sep = "|")
}

done_keys <- read_done(checkpoint)
if (!file.exists(checkpoint)) {
  cat(paste(header, collapse = ","), "\n", file = checkpoint, sep = "")
}

wilson_ci <- function(p, n, z = 1.96) {
  if (n <= 0) return(c(NA_real_, NA_real_))
  den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den
  c(max(0, center - half), min(1, center + half))
}

cat(sprintf("Independent R reproduction run_id=%s n_sim=%d cores=%d\n", run_id, n_sim, num.core))
cat("Label: R/mmrm reproduction using the settings reported in PAP_AS03\n")
cat("Generator: original data_sim.R / test_sim.R (Kenward-Roger / Kenward-Roger-Linear)\n")

cl <- parallel::makeCluster(num.core)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(mmrm)
    library(dplyr)
    library(MASS)
  })
  NULL
}))
parallel::clusterExport(cl, c("MMRM.test.1", "MMRM.sim.1"), envir = .GlobalEnv)

for (sc in scenarios) {
  for (ss in sample_sizes) {
    key <- paste(sc$id, ss, sep = "|")
    if (key %in% done_keys) {
      cat(sprintf("[skip] %s n=%d\n", sc$id, ss))
      next
    }
    t0 <- Sys.time()

    # Same seeding approach as original MMRM.sim.n: clusterSetRNGStream(seed).
    # Generate and test on the worker in one pass to avoid shipping trial data twice.
    fu1 <- sc$fu1
    fu2 <- sc$fu2
    parallel::clusterExport(cl, c("ss", "fu1", "fu2"), envir = environment())
    parallel::clusterSetRNGStream(cl, seed)
    p_vals <- unlist(parallel::parLapply(cl, seq_len(n_sim), function(i) {
      d <- MMRM.sim.1(
        nsim = i,
        ss = ss,
        trt.rate = 0.5,
        mu = 38,
        sd = 10.5,
        corr = 0.5,
        fu1 = fu1,
        fu2 = fu2,
        fu3 = NULL
      )
      MMRM.test.1(d)
    }))
    p_vals <- as.numeric(p_vals)
    ok <- is.finite(p_vals)
    m <- sum(ok)
    power <- if (m > 0) mean(p_vals[ok] < alpha) else NA_real_
    mcse <- if (m > 0) sqrt(power * (1 - power) / m) else NA_real_
    ci <- if (m > 0) wilson_ci(power, m) else c(NA_real_, NA_real_)

    row <- data.frame(
      scenario = sc$id,
      sample_size = ss,
      n_sim = n_sim,
      seed = seed,
      power_pct = 100 * power,
      mcse_power_pct = 100 * mcse,
      ci95_power_lo_pct = 100 * ci[1],
      ci95_power_hi_pct = 100 * ci[2],
      R_version = r_version,
      mmrm_version = mmrm_version,
      run_id = run_id,
      label = "R/mmrm reproduction using the settings reported in PAP_AS03",
      stringsAsFactors = FALSE
    )
    utils::write.table(
      row, file = checkpoint, sep = ",", row.names = FALSE, col.names = FALSE,
      append = TRUE, qmethod = "double"
    )
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf(
      "[%s n=%d] power=%.2f%% mcse=%.2f%% (%.1fs)\n",
      sc$id, ss, row$power_pct, row$mcse_power_pct, elapsed
    ))
    flush.console()
  }
}

final <- utils::read.csv(checkpoint, stringsAsFactors = FALSE)
final <- final[, header]
utils::write.csv(final, final_csv, row.names = FALSE)

manifest <- list(
  run_id = run_id,
  mode = if (smoke) "paper-reproduction-r-smoke" else "paper-reproduction-r",
  n_sim = n_sim,
  seed = seed,
  alpha = alpha,
  scenarios = vapply(scenarios, function(s) s$id, character(1)),
  sample_sizes = as.integer(sample_sizes),
  expected_cells = length(scenarios) * length(sample_sizes),
  completed_cells = nrow(final),
  r_method = "Kenward-Roger / Kenward-Roger-Linear (original test_sim.R)",
  generator = "original data_sim.R (not JavaScript)",
  label = "R/mmrm reproduction using the settings reported in PAP_AS03",
  R_version = r_version,
  mmrm_version = mmrm_version,
  output_file = basename(final_csv),
  status = if (nrow(final) == length(scenarios) * length(sample_sizes)) "COMPLETE" else "INCOMPLETE",
  notes = c(
    "Independent of the paired JS-R analysis.",
    "Do not label these values as exact paper / Figure 4 source data.",
    "Do not combine Monte Carlo errors with paired-run errors as if paired."
  )
)
jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
cat("Wrote", final_csv, "\n")
cat("Manifest", manifest_path, "\n")
if (manifest$status != "COMPLETE") quit(status = 2)
