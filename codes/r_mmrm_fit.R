#!/usr/bin/env Rscript
# Fit original EOSS_MMRM contrast with R mmrm on identical observed datasets.
# Input JSON (stdin or --file=): { alpha, methods: ["kenward-roger",...], trials: [{id, nf, rows:[{id,trt,y0,y1,y2,...}]}] }
# Output JSON (stdout): { session, results: [...] }

suppressPackageStartupMessages({
  ok <- requireNamespace("mmrm", quietly = TRUE) &&
        requireNamespace("jsonlite", quietly = TRUE) &&
        requireNamespace("dplyr", quietly = TRUE)
  if (!ok) {
    cat(jsonlite::toJSON(list(error = "Need packages: mmrm, jsonlite, dplyr"), auto_unbox = TRUE))
    quit(status = 1)
  }
})

args <- commandArgs(trailingOnly = TRUE)
infile <- NULL
for (a in args) {
  if (startsWith(a, "--file=")) infile <- sub("^--file=", "", a)
}

raw <- if (!is.null(infile)) paste(readLines(infile, warn = FALSE), collapse = "\n") else paste(readLines("stdin", warn = FALSE), collapse = "\n")
cfg <- jsonlite::fromJSON(raw, simplifyVector = FALSE)

alpha <- if (!is.null(cfg$alpha)) as.numeric(cfg$alpha) else 0.05
methods <- if (!is.null(cfg$methods)) unlist(cfg$methods) else c("kenward-roger")
trials <- cfg$trials

fit_one <- function(rows, nf, method) {
  t0 <- proc.time()[["elapsed"]]
  df <- tryCatch({
    d <- dplyr::bind_rows(lapply(rows, function(r) {
      # jsonlite null -> NULL; coerce to NA for numeric outcomes
      rr <- lapply(r, function(x) if (is.null(x)) NA else x)
      as.data.frame(rr, stringsAsFactors = FALSE)
    }))
    d$id <- as.integer(d$id)
    d$trt <- as.integer(d$trt)
    d$y0 <- as.numeric(d$y0)
    for (v in seq_len(nf)) {
      col <- paste0("y", v)
      d[[col]] <- as.numeric(d[[col]])
    }
    d
  }, error = function(e) NULL)

  if (is.null(df)) {
    return(list(
      ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
      failureReason = "data_parse_failed", converged = FALSE,
      runtime_seconds = proc.time()[["elapsed"]] - t0
    ))
  }

  varying <- paste0("y", seq_len(nf))
  longData <- tryCatch({
    ld <- reshape(df, varying = varying, direction = "long", sep = "", idvar = "id")
    ld <- ld[order(ld$id, ld$time), , drop = FALSE]
    ld$time <- factor(ld$time)
    ld$trt <- factor(ld$trt)
    ld$id <- factor(ld$id)
    ld
  }, error = function(e) NULL)

  if (is.null(longData)) {
    return(list(
      ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
      failureReason = "reshape_failed", converged = FALSE,
      runtime_seconds = proc.time()[["elapsed"]] - t0
    ))
  }

  ctrl <- if (identical(method, "kenward-roger")) {
    mmrm::mmrm_control(method = "Kenward-Roger", vcov = "Kenward-Roger-Linear")
  } else if (identical(method, "satterthwaite")) {
    mmrm::mmrm_control(method = "Satterthwaite")
  } else if (identical(method, "residual")) {
    # Closest available non-KR option in mmrm >= 0.3.x (no separate "Asymptotic" method name)
    mmrm::mmrm_control(method = "Residual")
  } else {
    mmrm::mmrm_control(method = "Kenward-Roger", vcov = "Kenward-Roger-Linear")
  }

  fit <- tryCatch(
    suppressWarnings(mmrm::mmrm(
      formula = y ~ y0 + trt + time + trt * time + y0 * time + us(time | id),
      data = longData,
      control = ctrl
    )),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
      failureReason = paste0("mmrm_error:", conditionMessage(fit)),
      converged = FALSE,
      runtime_seconds = proc.time()[["elapsed"]] - t0
    ))
  }

  beta <- tryCatch(mmrm::component(fit, "beta_est"), error = function(e) NULL)
  if (is.null(beta)) {
    return(list(
      ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
      failureReason = "beta_extract_failed", converged = FALSE,
      runtime_seconds = proc.time()[["elapsed"]] - t0
    ))
  }

  contrast <- numeric(length(beta))
  if (nf == 2) {
    # Original test_sim.R: contrast[c(3,5)] <- 1
    if (length(contrast) < 5) {
      return(list(ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
                  failureReason = "contrast_index_oob", converged = FALSE,
                  runtime_seconds = proc.time()[["elapsed"]] - t0))
    }
    contrast[c(3, 5)] <- 1
  } else {
    if (length(contrast) < 7) {
      return(list(ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
                  failureReason = "contrast_index_oob", converged = FALSE,
                  runtime_seconds = proc.time()[["elapsed"]] - t0))
    }
    contrast[c(3, 7)] <- 1
  }

  d1 <- tryCatch(mmrm::df_1d(fit, contrast), error = function(e) e)
  if (inherits(d1, "error")) {
    return(list(
      ok = FALSE, estimate = NA, standardError = NA, pValue = NA,
      failureReason = paste0("df_1d_error:", conditionMessage(d1)),
      converged = FALSE,
      runtime_seconds = proc.time()[["elapsed"]] - t0
    ))
  }

  est <- as.numeric(if (!is.null(d1$est)) d1$est else d1$estimate)
  se <- as.numeric(d1$se)
  pval <- as.numeric(d1$p_val)
  ok <- isTRUE(is.finite(est) && is.finite(se) && is.finite(pval) && se > 0)

  list(
    ok = ok,
    estimate = if (isTRUE(ok)) est else NA,
    standardError = if (isTRUE(ok)) se else NA,
    pValue = if (isTRUE(ok)) pval else NA,
    failureReason = if (isTRUE(ok)) NULL else "nonfinite_result",
    converged = isTRUE(ok),
    method = method,
    runtime_seconds = proc.time()[["elapsed"]] - t0
  )
}

results <- list()
for (tr in trials) {
  tid <- tr$id
  nf <- as.integer(tr$nf)
  rows <- tr$rows
  for (method in methods) {
    one <- fit_one(rows, nf, method)
    one$id <- tid
    one$method <- method
    results[[length(results) + 1]] <- one
  }
}

session <- list(
  r_version = R.version.string,
  mmrm_version = as.character(utils::packageVersion("mmrm")),
  platform = R.version$platform,
  methods = methods
)

out <- list(session = session, results = results)
# digits=NA keeps full numeric precision; default digits=4 rounded R estimates
# and made paired JS-minus-R differences uninterpretable at the 1e-5 scale.
cat(jsonlite::toJSON(out, auto_unbox = TRUE, null = "null", na = "null", digits = NA))
