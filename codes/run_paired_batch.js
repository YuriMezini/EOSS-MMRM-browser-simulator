#!/usr/bin/env node
/**
 * Paired JavaScript vs R mmrm evaluation runner.
 *
 * Example:
 *   node run_paired_batch.js --n-sim 50 --sample-sizes 40,100 --methods kenward-roger --seed 1234
 */

'use strict';

var fs = require('fs');
var path = require('path');
var os = require('os');
var { spawnSync } = require('child_process');
var crypto = require('crypto');

var engine = require('./src/eoss_engine');
var scenarios = require('./src/scenarios');
var { summarizePaired } = require('./src/paired_statistics');

var ROOT = __dirname;
var OUT_DIR = path.join(ROOT, 'paper', 'results');
var R_SCRIPT = path.join(ROOT, 'src', 'r_mmrm_fit.R');

function parseArgs(argv) {
  var args = {
    nSim: 50,
    seed: scenarios.COMMON.baseSeed,
    alpha: scenarios.COMMON.alpha,
    sampleSizes: [40, 100],
    scenarioIds: scenarios.SCENARIOS.map(function (s) { return s.id; }),
    methods: ['kenward-roger'],
    mode: 'paired-custom',
    rBatchSize: 25,
    outPrefix: 'paired'
  };
  for (var i = 2; i < argv.length; i++) {
    var a = argv[i];
    if (a === '--n-sim') args.nSim = parseInt(argv[++i], 10);
    else if (a === '--seed') args.seed = parseInt(argv[++i], 10);
    else if (a === '--alpha') args.alpha = parseFloat(argv[++i]);
    else if (a === '--sample-sizes') {
      args.sampleSizes = argv[++i].split(',').map(function (x) { return parseInt(x.trim(), 10); });
    } else if (a === '--scenarios') {
      args.scenarioIds = argv[++i].split(',').map(function (x) { return x.trim(); });
    } else if (a === '--methods') {
      args.methods = argv[++i].split(',').map(function (x) { return x.trim(); });
    } else if (a === '--mode') args.mode = argv[++i];
    else if (a === '--r-batch-size') args.rBatchSize = parseInt(argv[++i], 10);
    else if (a === '--out-prefix') args.outPrefix = argv[++i];
  }
  return args;
}

function ensureDir(d) {
  if (!fs.existsSync(d)) fs.mkdirSync(d, { recursive: true });
}

function relPath(p) {
  return path.relative(ROOT, p).split(path.sep).join('/');
}

function sha256File(p) {
  if (!fs.existsSync(p)) return null;
  var h = crypto.createHash('sha256');
  h.update(fs.readFileSync(p));
  return h.digest('hex');
}

function callRBatch(trials, alpha, methods) {
  var tmp = path.join(os.tmpdir(), 'eoss_mmrm_paired_' + process.pid + '_' + Date.now() + '.json');
  var payload = { alpha: alpha, methods: methods, trials: trials };
  fs.writeFileSync(tmp, JSON.stringify(payload));
  var res = spawnSync('Rscript', [R_SCRIPT, '--file=' + tmp], {
    encoding: 'utf8',
    maxBuffer: 64 * 1024 * 1024
  });
  try { fs.unlinkSync(tmp); } catch (e) { /* ignore */ }
  if (res.status !== 0) {
    throw new Error('Rscript failed: ' + (res.stderr || res.stdout || ('status ' + res.status)));
  }
  var out = JSON.parse(res.stdout);
  return out;
}

function runCell(scenario, sampleSize, args) {
  var params = scenarios.scenarioParams(scenario, sampleSize);
  var nf = params.fus.length;
  var pairedRows = [];
  var tCell0 = Date.now();
  var rSession = null;

  for (var start = 0; start < args.nSim; start += args.rBatchSize) {
    var end = Math.min(args.nSim, start + args.rBatchSize);
    var batchTrials = [];
    var batchMeta = [];

    for (var i = start; i < end; i++) {
      var seed = engine.trialSeed(args.seed, i);
      var t0 = process.hrtime.bigint();
      var structured = engine.simulateTrialStructured(params, seed);
      var observed = engine.toObservedWide(structured, nf);
      var legacy = engine.simTrial(params, seed);
      // Identical-dataset guard (legacy vs structured observed)
      for (var k = 0; k < observed.length; k++) {
        if (observed[k].y1 !== legacy[k].y1 || observed[k].y2 !== legacy[k].y2 ||
            observed[k].trt !== legacy[k].trt || observed[k].y0 !== legacy[k].y0) {
          throw new Error('LEGACY/STRUCTURED mismatch at replicate ' + i + ' subject ' + observed[k].id);
        }
      }
      var tJs0 = process.hrtime.bigint();
      var jsFit = engine.fitMMRM(observed, nf);
      var runtimeJs = Number(process.hrtime.bigint() - tJs0) / 1e9;

      var empFull = engine.empiricalFullDataTarget(structured, nf);
      var empObs = engine.empiricalObservedDataTarget(structured, nf);

      batchMeta.push({
        replicate: i,
        seed: seed,
        structured: structured,
        observed: observed,
        jsFit: jsFit,
        runtimeJs: runtimeJs,
        empFull: empFull,
        empObs: empObs,
        genNs: Number(process.hrtime.bigint() - t0) / 1e9
      });
      batchTrials.push({
        id: String(i),
        nf: nf,
        rows: observed
      });
    }

    var rOut = callRBatch(batchTrials, args.alpha, args.methods);
    if (!rSession) rSession = rOut.session;
    var byKey = {};
    (rOut.results || []).forEach(function (r) {
      byKey[r.id + '|' + r.method] = r;
    });

    var primaryMethod = args.methods[0];
    batchMeta.forEach(function (m) {
      var r = byKey[String(m.replicate) + '|' + primaryMethod] || {};
      var jsOk = !!m.jsFit.ok;
      var rOk = r.ok === true || r.ok === 'true';
      var jsEst = jsOk ? m.jsFit.estimate : null;
      var rEst = rOk && r.estimate != null ? Number(r.estimate) : null;
      var jsSe = jsOk ? m.jsFit.standardError : null;
      var rSe = rOk && r.standardError != null ? Number(r.standardError) : null;
      var jsP = jsOk ? m.jsFit.pValue : null;
      var rP = rOk && r.pValue != null ? Number(r.pValue) : null;
      var jsReject = jsOk && jsP < args.alpha;
      var rReject = rOk && rP < args.alpha;

      pairedRows.push({
        scenario: scenario.id,
        sample_size: sampleSize,
        replicate: m.replicate,
        seed: m.seed,
        true_full_data_target: scenario.true_full_data_target,
        observed_data_target: scenario.analytic_observed_mixture_target,
        empirical_full_data_md: m.empFull,
        empirical_observed_data_md: m.empObs,
        js_estimate: jsEst,
        r_estimate: rEst,
        estimate_difference: jsOk && rOk ? jsEst - rEst : null,
        js_se: jsSe,
        r_se: rSe,
        se_difference: jsOk && rOk ? jsSe - rSe : null,
        js_pvalue: jsP,
        r_pvalue: rP,
        js_reject: jsReject,
        r_reject: rReject,
        decision_disagreement: jsOk && rOk ? (jsReject !== rReject) : null,
        js_ok: jsOk,
        r_ok: rOk,
        js_convergence: jsOk ? !!m.jsFit.converged : false,
        r_convergence: rOk ? !!r.converged : false,
        js_failure_reason: jsOk ? null : (m.jsFit.failureReason || 'js_fail'),
        r_failure_reason: rOk ? null : (r.failureReason || 'r_fail'),
        js_iterations: m.jsFit.iterations != null ? m.jsFit.iterations : null,
        runtime_js: m.runtimeJs,
        runtime_r: r.runtime_seconds != null ? Number(r.runtime_seconds) : null,
        r_method: primaryMethod,
        dataset_fingerprint: crypto.createHash('sha256')
          .update(JSON.stringify(m.observed))
          .digest('hex')
          .slice(0, 16)
      });
    });
  }

  var summary = summarizePaired(pairedRows, {
    scenario: scenario.id,
    sample_size: sampleSize,
    n_sim: args.nSim,
    seed: args.seed,
    alpha: args.alpha,
    true_full_data_target: scenario.true_full_data_target,
    analytic_observed_mixture_target: scenario.analytic_observed_mixture_target,
    stored_reference_power_pct: scenarios.storedReferencePowerPct(scenario.id, sampleSize),
    r_method: args.methods[0]
  });
  summary.cell_runtime_seconds = (Date.now() - tCell0) / 1000;
  summary.r_session = rSession;

  return { rows: pairedRows, summary: summary, rSession: rSession };
}

function writeCsv(file, rows) {
  if (!rows.length) {
    fs.writeFileSync(file, '');
    return;
  }
  var cols = Object.keys(rows[0]);
  var lines = [cols.join(',')];
  rows.forEach(function (r) {
    lines.push(cols.map(function (c) {
      var v = r[c];
      if (v === null || v === undefined) return 'NA';
      if (typeof v === 'object') return '"' + JSON.stringify(v).replace(/"/g, '""') + '"';
      if (typeof v === 'string' && /[",\n]/.test(v)) return '"' + v.replace(/"/g, '""') + '"';
      return String(v);
    }).join(','));
  });
  fs.writeFileSync(file, lines.join('\n') + '\n');
}

function flattenSummary(s) {
  return {
    scenario: s.scenario,
    sample_size: s.sample_size,
    n_sim_requested: s.n_sim_requested,
    seed: s.seed,
    alpha: s.alpha,
    r_method: s.r_method,
    true_full_data_target: s.true_full_data_target,
    analytic_observed_mixture_target: s.analytic_observed_mixture_target,
    paired_n: s.paired_n,
    mean_estimate_difference_js_minus_r: s.mean_estimate_difference_js_minus_r,
    mcse_estimate_difference: s.mcse_estimate_difference,
    ci95_estimate_difference_lo: s.ci95_estimate_difference_lo,
    ci95_estimate_difference_hi: s.ci95_estimate_difference_hi,
    mean_abs_difference: s.mean_abs_difference,
    paired_rmse: s.paired_rmse,
    median_abs_difference: s.median_abs_difference,
    max_abs_difference: s.max_abs_difference,
    estimate_correlation: s.estimate_correlation,
    mean_js_over_r_se_ratio: s.mean_js_over_r_se_ratio,
    mean_abs_pvalue_difference: s.mean_abs_pvalue_difference,
    rejection_disagreement_rate: s.rejection_disagreement_rate,
    mean_paired_reject_diff_js_minus_r: s.mean_paired_reject_diff_js_minus_r,
    mcse_paired_reject_diff: s.mcse_paired_reject_diff,
    js_rejection_rate_pct: s.js_rejection_rate == null ? null : 100 * s.js_rejection_rate,
    r_rejection_rate_pct: s.r_rejection_rate == null ? null : 100 * s.r_rejection_rate,
    stored_reference_power_pct: s.stored_reference_power_pct,
    power_difference_js_vs_stored_reference_pp: s.js_rejection_rate == null || s.stored_reference_power_pct == null
      ? null
      : 100 * s.js_rejection_rate - s.stored_reference_power_pct,
    power_difference_r_vs_stored_reference_pp: s.r_rejection_rate == null || s.stored_reference_power_pct == null
      ? null
      : 100 * s.r_rejection_rate - s.stored_reference_power_pct,
    A_js_mean_bias_vs_full: s.A_js_bias_vs_full.mean_bias,
    A_js_mcse: s.A_js_bias_vs_full.mcse,
    A_r_mean_bias_vs_full: s.A_r_bias_vs_full.mean_bias,
    A_r_mcse: s.A_r_bias_vs_full.mcse,
    B_js_mean_bias_vs_observed: s.B_js_bias_vs_observed.mean_bias,
    B_js_mcse: s.B_js_bias_vs_observed.mcse,
    B_r_mean_bias_vs_observed: s.B_r_bias_vs_observed.mean_bias,
    B_r_mcse: s.B_r_bias_vs_observed.mcse,
    js_wald_coverage_vs_full_pct: s.js_wald_coverage_vs_full == null ? null : 100 * s.js_wald_coverage_vs_full,
    r_wald_coverage_vs_full_pct: s.r_wald_coverage_vs_full == null ? null : 100 * s.r_wald_coverage_vs_full,
    js_valid_fits: s.js_valid_fits,
    r_valid_fits: s.r_valid_fits,
    js_failed_fits: s.js_failed_fits,
    r_failed_fits: s.r_failed_fits,
    median_runtime_js: s.median_runtime_js,
    median_runtime_r: s.median_runtime_r,
    median_runtime_ratio_r_over_js: s.median_runtime_ratio_r_over_js,
    cell_runtime_seconds: s.cell_runtime_seconds,
    timestamp: s.timestamp
  };
}

function collectFailures(allRows) {
  var map = {};
  allRows.forEach(function (r) {
    function add(engineName, reason) {
      if (!reason) return;
      var key = r.scenario + '|' + r.sample_size + '|' + engineName + '|' + reason;
      if (!map[key]) {
        map[key] = {
          scenario: r.scenario,
          sample_size: r.sample_size,
          engine: engineName,
          failure_reason: reason,
          count: 0
        };
      }
      map[key].count++;
    }
    if (!r.js_ok) add('javascript', r.js_failure_reason);
    if (!r.r_ok) add('r_mmrm', r.r_failure_reason);
  });
  return Object.keys(map).map(function (k) { return map[k]; });
}

async function main() {
  var args = parseArgs(process.argv);
  ensureDir(OUT_DIR);

  var selected = scenarios.SCENARIOS.filter(function (s) {
    return args.scenarioIds.indexOf(s.id) >= 0;
  });
  if (!selected.length) throw new Error('No scenarios selected');

  console.log('Paired JS–R evaluation');
  console.log('mode=%s n_sim=%d sizes=%s methods=%s',
    args.mode, args.nSim, args.sampleSizes.join(','), args.methods.join(','));
  console.log('NOTE: compound-symmetry data generation only.');

  var allRows = [];
  var summaries = [];
  var rSession = null;
  var t0 = Date.now();
  var combos = selected.length * args.sampleSizes.length;
  var done = 0;

  for (var si = 0; si < selected.length; si++) {
    for (var zi = 0; zi < args.sampleSizes.length; zi++) {
      var sc = selected[si];
      var ss = args.sampleSizes[zi];
      var cell = runCell(sc, ss, args);
      if (!rSession) rSession = cell.rSession;
      allRows = allRows.concat(cell.rows);
      summaries.push(flattenSummary(cell.summary));
      done++;
      console.log(
        '[' + done + '/' + combos + '] ' + sc.id + ' n=' + ss +
        ' paired=' + cell.summary.paired_n +
        ' meanΔ=' + (cell.summary.mean_estimate_difference_js_minus_r == null
          ? 'NA'
          : cell.summary.mean_estimate_difference_js_minus_r.toFixed(6)) +
        ' disag=' + (cell.summary.rejection_disagreement_rate == null
          ? 'NA'
          : cell.summary.rejection_disagreement_rate.toFixed(3)) +
        ' (' + cell.summary.cell_runtime_seconds.toFixed(1) + 's)'
      );

      // incremental write
      writeCsv(path.join(OUT_DIR, args.outPrefix + '_results.csv'), allRows);
      writeCsv(path.join(OUT_DIR, args.outPrefix + '_summary.csv'), summaries);
    }
  }

  var flatBias = summaries.map(function (s) {
    return {
      scenario: s.scenario,
      sample_size: s.sample_size,
      n_sim: s.n_sim_requested,
      A_js_mean_bias_vs_full: s.A_js_mean_bias_vs_full,
      A_js_mcse: s.A_js_mcse,
      A_r_mean_bias_vs_full: s.A_r_mean_bias_vs_full,
      A_r_mcse: s.A_r_mcse,
      B_js_mean_bias_vs_observed: s.B_js_mean_bias_vs_observed,
      B_js_mcse: s.B_js_mcse,
      B_r_mean_bias_vs_observed: s.B_r_mean_bias_vs_observed,
      B_r_mcse: s.B_r_mcse,
      note: 'A=full-data estimand; B=observed-data target; not the same as C (paired JS-R)'
    };
  });
  var powerRows = summaries.map(function (s) {
    return {
      scenario: s.scenario,
      sample_size: s.sample_size,
      n_sim: s.n_sim_requested,
      js_rejection_rate_pct: s.js_rejection_rate_pct,
      r_rejection_rate_pct: s.r_rejection_rate_pct,
      stored_reference_power_pct: s.stored_reference_power_pct,
      power_difference_js_vs_stored_reference_pp: s.power_difference_js_vs_stored_reference_pp,
      power_difference_r_vs_stored_reference_pp: s.power_difference_r_vs_stored_reference_pp,
      rejection_disagreement_rate: s.rejection_disagreement_rate,
      label: 'power differences vs paper are not estimator bias'
    };
  });
  var runtimeRows = summaries.map(function (s) {
    return {
      scenario: s.scenario,
      sample_size: s.sample_size,
      median_runtime_js: s.median_runtime_js,
      median_runtime_r: s.median_runtime_r,
      median_runtime_ratio_r_over_js: s.median_runtime_ratio_r_over_js,
      cell_runtime_seconds: s.cell_runtime_seconds,
      note: 'Node JS fit time vs R mmrm fit time; not a browser benchmark'
    };
  });
  var failures = collectFailures(allRows);

  var prefix = args.outPrefix;
  var files = {
    results: path.join(OUT_DIR, prefix + '_results.csv'),
    summary: path.join(OUT_DIR, prefix + '_summary.csv'),
    bias: path.join(OUT_DIR, prefix + '_bias_summary.csv'),
    power: path.join(OUT_DIR, prefix + '_power_summary.csv'),
    runtime: path.join(OUT_DIR, prefix + '_runtime_summary.csv'),
    failures: path.join(OUT_DIR, prefix + '_failures.csv')
  };
  // Also write final_* aliases when mode is paired-final / paired-pilot
  if (args.mode.indexOf('final') >= 0 || prefix === 'final_paired') {
    files = {
      results: path.join(OUT_DIR, 'final_paired_results.csv'),
      summary: path.join(OUT_DIR, 'final_paired_summary.csv'),
      bias: path.join(OUT_DIR, 'final_bias_summary.csv'),
      power: path.join(OUT_DIR, 'final_power_summary.csv'),
      runtime: path.join(OUT_DIR, 'final_runtime_summary.csv'),
      failures: path.join(OUT_DIR, 'final_failures.csv')
    };
  }

  writeCsv(files.results, allRows);
  writeCsv(files.summary, summaries);
  writeCsv(files.bias, flatBias);
  writeCsv(files.power, powerRows);
  writeCsv(files.runtime, runtimeRows);
  writeCsv(files.failures, failures);

  var cpus = os.cpus() || [];
  var meta = {
    mode: args.mode,
    n_sim: args.nSim,
    seed: args.seed,
    alpha: args.alpha,
    sample_sizes: args.sampleSizes,
    scenarios: args.scenarioIds,
    methods: args.methods,
    covariance_generation: 'compound_symmetry_only',
    node_version: process.version,
    platform: process.platform + ' ' + os.release(),
    cpu_model: cpus[0] ? cpus[0].model : null,
    logical_cpus: cpus.length,
    total_memory_bytes: os.totalmem(),
    r_session: rSession,
    elapsed_seconds: (Date.now() - t0) / 1000,
    output_files: Object.keys(files).reduce(function (acc, k) {
      acc[k] = relPath(files[k]);
      return acc;
    }, {}),
    checksums_sha256: Object.keys(files).reduce(function (acc, k) {
      acc[k] = sha256File(files[k]);
      return acc;
    }, {}),
    distinctions: {
      A: 'bias versus prespecified full-data estimand',
      B: 'bias versus analytic observed-data mixture target',
      C: 'paired JS minus R estimate on identical observed datasets'
    },
    browser_benchmark: 'NOT_RUN',
    notes: [
      'Installed mmrm version is recorded in r_session; original app comments mention 0.3.6/0.3.11.',
      'Normal Wald intervals (±1.96 SE) are used for coverage; JS default p-values use the HTML custom correction.',
      'No absolute local paths are stored in output_files.'
    ]
  };
  var metaPath = path.join(OUT_DIR, 'run_metadata.json');
  fs.writeFileSync(metaPath, JSON.stringify(meta, null, 2));
  console.log('Wrote metadata', relPath(metaPath));
  console.log('Elapsed ' + meta.elapsed_seconds.toFixed(1) + 's');
}

if (require.main === module) {
  main().catch(function (err) {
    console.error(err);
    process.exit(1);
  });
}

module.exports = { parseArgs: parseArgs, runCell: runCell };
