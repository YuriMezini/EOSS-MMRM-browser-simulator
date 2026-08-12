/**
 * Paired JS vs R summary metrics.
 * A = bias vs prespecified full-data estimand
 * B = bias vs observed-data target
 * C = paired JS−R computational difference
 * Never label B or C as A.
 */

'use strict';

function isFin(x) {
  return typeof x === 'number' && Number.isFinite(x);
}

function mean(xs) {
  if (!xs.length) return null;
  var s = 0;
  for (var i = 0; i < xs.length; i++) s += xs[i];
  return s / xs.length;
}

function sd(xs) {
  if (xs.length < 2) return null;
  var m = mean(xs);
  var s = 0;
  for (var i = 0; i < xs.length; i++) {
    var d = xs[i] - m;
    s += d * d;
  }
  return Math.sqrt(s / (xs.length - 1));
}

function quantile(xs, q) {
  if (!xs.length) return null;
  var a = xs.slice().sort(function (u, v) { return u - v; });
  var pos = (a.length - 1) * q;
  var base = Math.floor(pos);
  var rest = pos - base;
  if (a[base + 1] !== undefined) return a[base] + rest * (a[base + 1] - a[base]);
  return a[base];
}

function corr(xs, ys) {
  var n = Math.min(xs.length, ys.length);
  if (n < 2) return null;
  var mx = mean(xs), my = mean(ys);
  var num = 0, dx = 0, dy = 0;
  for (var i = 0; i < n; i++) {
    var a = xs[i] - mx, b = ys[i] - my;
    num += a * b;
    dx += a * a;
    dy += b * b;
  }
  if (dx <= 0 || dy <= 0) return null;
  return num / Math.sqrt(dx * dy);
}

/**
 * @param {object[]} rows paired trial rows
 */
function summarizePaired(rows, meta) {
  var diffs = [], absDiffs = [], seRatios = [], absP = [], D = [];
  var jsEst = [], rEst = [], jsOk = 0, rOk = 0, bothOk = 0;
  var jsFail = {}, rFail = {};
  var runtimeJs = [], runtimeR = [];
  var biasFullJs = [], biasObsJs = [], biasFullR = [], biasObsR = [];
  var jsRej = 0, rRej = 0;
  var coveredJs = 0, coveredR = 0;
  var alpha = meta.alpha;
  var trueFull = meta.true_full_data_target;

  for (var i = 0; i < rows.length; i++) {
    var row = rows[i];
    if (row.js_ok) jsOk++;
    else {
      var jf = row.js_failure_reason || 'unknown';
      jsFail[jf] = (jsFail[jf] || 0) + 1;
    }
    if (row.r_ok) rOk++;
    else {
      var rf = row.r_failure_reason || 'unknown';
      rFail[rf] = (rFail[rf] || 0) + 1;
    }
    if (isFin(row.runtime_js)) runtimeJs.push(row.runtime_js);
    if (isFin(row.runtime_r)) runtimeR.push(row.runtime_r);

    if (row.js_ok && row.r_ok && isFin(row.js_estimate) && isFin(row.r_estimate)) {
      bothOk++;
      var d = row.js_estimate - row.r_estimate;
      diffs.push(d);
      absDiffs.push(Math.abs(d));
      jsEst.push(row.js_estimate);
      rEst.push(row.r_estimate);
      if (isFin(row.js_se) && isFin(row.r_se) && row.r_se > 0) {
        seRatios.push(row.js_se / row.r_se);
      }
      if (isFin(row.js_pvalue) && isFin(row.r_pvalue)) {
        absP.push(Math.abs(row.js_pvalue - row.r_pvalue));
      }
      var dj = row.js_reject ? 1 : 0;
      var dr = row.r_reject ? 1 : 0;
      D.push(dj - dr);
      if (row.js_reject) jsRej++;
      if (row.r_reject) rRej++;
    }

    if (row.js_ok && isFin(row.js_estimate) && isFin(trueFull)) {
      biasFullJs.push(row.js_estimate - trueFull);
      if (isFin(row.js_se)) {
        var lo = row.js_estimate - 1.96 * row.js_se;
        var hi = row.js_estimate + 1.96 * row.js_se;
        if (trueFull >= lo && trueFull <= hi) coveredJs++;
      }
    }
    if (row.r_ok && isFin(row.r_estimate) && isFin(trueFull)) {
      biasFullR.push(row.r_estimate - trueFull);
      if (isFin(row.r_se)) {
        lo = row.r_estimate - 1.96 * row.r_se;
        hi = row.r_estimate + 1.96 * row.r_se;
        if (trueFull >= lo && trueFull <= hi) coveredR++;
      }
    }
    if (row.js_ok && isFin(row.js_estimate) && isFin(row.observed_data_target)) {
      biasObsJs.push(row.js_estimate - row.observed_data_target);
    }
    if (row.r_ok && isFin(row.r_estimate) && isFin(row.observed_data_target)) {
      biasObsR.push(row.r_estimate - row.observed_data_target);
    }
  }

  var meanDiff = mean(diffs);
  var sdDiff = sd(diffs);
  var mcseDiff = sdDiff != null && diffs.length ? sdDiff / Math.sqrt(diffs.length) : null;
  var meanD = mean(D);
  var sdD = sd(D);
  var mcseD = sdD != null && D.length ? sdD / Math.sqrt(D.length) : null;

  function biasBlock(arr, label) {
    var m = mean(arr);
    var s = sd(arr);
    var n = arr.length;
    var mcse = s != null && n ? s / Math.sqrt(n) : null;
    return {
      label: label,
      n: n,
      mean_bias: m,
      mcse: mcse,
      ci95_lo: m != null && mcse != null ? m - 1.96 * mcse : null,
      ci95_hi: m != null && mcse != null ? m + 1.96 * mcse : null,
      rmse: n ? Math.sqrt(mean(arr.map(function (x) { return x * x; }))) : null
    };
  }

  return {
    scenario: meta.scenario,
    sample_size: meta.sample_size,
    n_sim_requested: meta.n_sim,
    seed: meta.seed,
    alpha: alpha,
    r_method: meta.r_method || 'kenward-roger',
    true_full_data_target: trueFull,
    analytic_observed_mixture_target: meta.analytic_observed_mixture_target,
    // C: paired computational difference
    paired_n: bothOk,
    mean_estimate_difference_js_minus_r: meanDiff,
    mcse_estimate_difference: mcseDiff,
    ci95_estimate_difference_lo: meanDiff != null && mcseDiff != null ? meanDiff - 1.96 * mcseDiff : null,
    ci95_estimate_difference_hi: meanDiff != null && mcseDiff != null ? meanDiff + 1.96 * mcseDiff : null,
    mean_abs_difference: mean(absDiffs),
    paired_rmse: absDiffs.length ? Math.sqrt(mean(absDiffs.map(function (x) { return x * x; }))) : null,
    median_abs_difference: quantile(absDiffs, 0.5),
    max_abs_difference: absDiffs.length ? Math.max.apply(null, absDiffs) : null,
    estimate_correlation: corr(jsEst, rEst),
    mean_js_over_r_se_ratio: mean(seRatios),
    mean_abs_pvalue_difference: mean(absP),
    rejection_disagreement_rate: D.length ? mean(D.map(function (x) { return Math.abs(x); })) : null,
    mean_paired_reject_diff_js_minus_r: meanD,
    mcse_paired_reject_diff: mcseD,
    js_rejection_rate: bothOk ? jsRej / bothOk : null,
    r_rejection_rate: bothOk ? rRej / bothOk : null,
    stored_reference_power_pct: meta.stored_reference_power_pct,
    // A / B blocks
    A_js_bias_vs_full: biasBlock(biasFullJs, 'A_js_vs_full'),
    A_r_bias_vs_full: biasBlock(biasFullR, 'A_r_vs_full'),
    B_js_bias_vs_observed: biasBlock(biasObsJs, 'B_js_vs_observed'),
    B_r_bias_vs_observed: biasBlock(biasObsR, 'B_r_vs_observed'),
    js_wald_coverage_vs_full: biasFullJs.length ? coveredJs / biasFullJs.length : null,
    r_wald_coverage_vs_full: biasFullR.length ? coveredR / biasFullR.length : null,
    js_valid_fits: jsOk,
    r_valid_fits: rOk,
    js_failed_fits: rows.length - jsOk,
    r_failed_fits: rows.length - rOk,
    js_failure_reasons: jsFail,
    r_failure_reasons: rFail,
    median_runtime_js: quantile(runtimeJs, 0.5),
    median_runtime_r: quantile(runtimeR, 0.5),
    iqr_runtime_js: runtimeJs.length ? quantile(runtimeJs, 0.75) - quantile(runtimeJs, 0.25) : null,
    iqr_runtime_r: runtimeR.length ? quantile(runtimeR, 0.75) - quantile(runtimeR, 0.25) : null,
    median_runtime_ratio_r_over_js: (function () {
      var mjs = quantile(runtimeJs, 0.5), mr = quantile(runtimeR, 0.5);
      if (mjs && mr && mjs > 0) return mr / mjs;
      return null;
    })(),
    timestamp: new Date().toISOString()
  };
}

module.exports = {
  summarizePaired: summarizePaired,
  mean: mean,
  sd: sd,
  corr: corr,
  quantile: quantile
};
