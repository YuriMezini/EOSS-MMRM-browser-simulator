/**
 * Model-bias summary statistics.
 * Formulas follow MODEL_BIAS_VALIDATION_GUIDE.md.
 */

'use strict';

var RESULT_COLUMNS = [
  'scenario',
  'target_definition',
  'sample_size',
  'n_sim_requested',
  'seed',
  'alpha',
  'true_effect',
  'mean_estimate',
  'bias',
  'relative_bias_pct',
  'empirical_sd',
  'mean_model_se',
  'se_sd_ratio',
  'rmse',
  'coverage_95_pct',
  'rejection_rate_pct',
  'rejection_label',
  'valid_fits',
  'failed_fits',
  'failure_rate_pct',
  'mcse_bias',
  'mcse_coverage_pct',
  'mcse_rejection_pct',
  'stored_reference_power_pct',
  'power_difference_pp',
  'runtime_seconds',
  'timestamp'
];

function isFiniteNumber(x) {
  return typeof x === 'number' && Number.isFinite(x);
}

/**
 * @param {object} opts
 * @param {object} opts.scenario - scenario definition
 * @param {number} opts.sampleSize
 * @param {number} opts.nSim
 * @param {number} opts.seed
 * @param {number} opts.alpha
 * @param {Array<{ok:boolean,estimate:number|null,standardError:number|null,pValue:number|null}>} opts.replicates
 * @param {number|null} opts.storedReferencePowerPct
 * @param {number} opts.runtimeSeconds
 */
function summarizeBias(opts) {
  var scenario = opts.scenario;
  var replicates = opts.replicates;
  var trueEffect = scenario.true_effect;
  var nSim = opts.nSim;
  var alpha = opts.alpha;

  var valid = [];
  var failed = 0;
  for (var i = 0; i < replicates.length; i++) {
    var r = replicates[i];
    if (r && r.ok && isFiniteNumber(r.estimate) && isFiniteNumber(r.standardError) && isFiniteNumber(r.pValue)) {
      valid.push(r);
    } else {
      failed++;
    }
  }

  var nValid = valid.length;
  var meanEstimate = null;
  var bias = null;
  var relativeBiasPct = null;
  var empiricalSd = null;
  var meanModelSe = null;
  var seSdRatio = null;
  var rmse = null;
  var coverage = null;
  var rejectionRate = null;
  var mcseBias = null;
  var mcseCoveragePct = null;
  var mcseRejectionPct = null;

  if (nValid > 0) {
    var sumEst = 0, sumSe = 0, sumSqErr = 0, covered = 0, rejected = 0;
    for (i = 0; i < nValid; i++) {
      sumEst += valid[i].estimate;
      sumSe += valid[i].standardError;
      var err = valid[i].estimate - trueEffect;
      sumSqErr += err * err;
      var lo = valid[i].estimate - 1.96 * valid[i].standardError;
      var hi = valid[i].estimate + 1.96 * valid[i].standardError;
      if (trueEffect >= lo && trueEffect <= hi) covered++;
      if (valid[i].pValue < alpha) rejected++;
    }
    meanEstimate = sumEst / nValid;
    bias = meanEstimate - trueEffect;
    relativeBiasPct = trueEffect === 0 ? null : 100 * bias / Math.abs(trueEffect);
    meanModelSe = sumSe / nValid;
    rmse = Math.sqrt(sumSqErr / nValid);
    coverage = covered / nValid;
    rejectionRate = rejected / nValid;

    var sumSqDev = 0;
    for (i = 0; i < nValid; i++) {
      var d = valid[i].estimate - meanEstimate;
      sumSqDev += d * d;
    }
    empiricalSd = Math.sqrt(sumSqDev / Math.max(nValid - 1, 1));
    seSdRatio = empiricalSd > 0 ? meanModelSe / empiricalSd : null;
    mcseBias = empiricalSd / Math.sqrt(nValid);
    mcseCoveragePct = 100 * Math.sqrt(coverage * (1 - coverage) / nValid);
    mcseRejectionPct = 100 * Math.sqrt(rejectionRate * (1 - rejectionRate) / nValid);
  }

  var storedReferencePowerPct = opts.storedReferencePowerPct == null ? null : opts.storedReferencePowerPct;
  var powerDifferencePp = null;
  if (storedReferencePowerPct != null && rejectionRate != null) {
    powerDifferencePp = 100 * rejectionRate - storedReferencePowerPct;
  }

  return {
    scenario: scenario.id,
    target_definition: scenario.target_definition,
    sample_size: opts.sampleSize,
    n_sim_requested: nSim,
    seed: opts.seed,
    alpha: alpha,
    true_effect: trueEffect,
    mean_estimate: meanEstimate,
    bias: bias,
    relative_bias_pct: relativeBiasPct,
    empirical_sd: empiricalSd,
    mean_model_se: meanModelSe,
    se_sd_ratio: seSdRatio,
    rmse: rmse,
    coverage_95_pct: coverage == null ? null : 100 * coverage,
    rejection_rate_pct: rejectionRate == null ? null : 100 * rejectionRate,
    rejection_label: scenario.rejection_label,
    valid_fits: nValid,
    failed_fits: failed,
    failure_rate_pct: nSim > 0 ? 100 * failed / nSim : null,
    mcse_bias: mcseBias,
    mcse_coverage_pct: mcseCoveragePct,
    mcse_rejection_pct: mcseRejectionPct,
    stored_reference_power_pct: storedReferencePowerPct,
    power_difference_pp: powerDifferencePp,
    runtime_seconds: opts.runtimeSeconds,
    timestamp: new Date().toISOString()
  };
}

function csvEscape(value) {
  if (value === null || value === undefined) return 'NA';
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) return 'NA';
    return String(value);
  }
  var s = String(value);
  if (/[",\n\r]/.test(s)) return '"' + s.replace(/"/g, '""') + '"';
  return s;
}

function rowToCsv(row) {
  return RESULT_COLUMNS.map(function (col) { return csvEscape(row[col]); }).join(',');
}

function csvHeader() {
  return RESULT_COLUMNS.join(',');
}

module.exports = {
  RESULT_COLUMNS: RESULT_COLUMNS,
  summarizeBias: summarizeBias,
  csvEscape: csvEscape,
  rowToCsv: rowToCsv,
  csvHeader: csvHeader,
  isFiniteNumber: isFiniteNumber
};
