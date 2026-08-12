#!/usr/bin/env node
/**
 * Re-export R mmrm fits at full JSON precision for an existing paired-results CSV.
 * Regenerates identical observed datasets from (scenario, sample_size, seed) and
 * replaces r_estimate / r_se / r_pvalue / differences while keeping JS fields.
 *
 * Example:
 *   node codes/reexport_r_full_precision.js \
 *     --input results/final/.../final_paired_results.csv \
 *     --outdir results/final/..._fullprec \
 *     --workers 6 --batch-size 40
 */
'use strict';

var fs = require('fs');
var path = require('path');
var os = require('os');
var readline = require('readline');
var { spawn } = require('child_process');

var engine = require('./eoss_engine');
var scenarios = require('./scenarios');
var { summarizePaired } = require('./paired_statistics');

var ROOT = path.join(__dirname, '..');
var R_SCRIPT = path.join(__dirname, 'r_mmrm_fit.R');

function parseArgs(argv) {
  var args = {
    input: null,
    outdir: null,
    workers: Math.max(1, Math.min(6, (os.cpus().length || 2) - 1)),
    batchSize: 40,
    scenarios: null,
    sampleSizes: null,
    limit: null
  };
  for (var i = 2; i < argv.length; i++) {
    var a = argv[i];
    if (a === '--input') args.input = argv[++i];
    else if (a === '--outdir') args.outdir = argv[++i];
    else if (a === '--workers') args.workers = parseInt(argv[++i], 10);
    else if (a === '--batch-size') args.batchSize = parseInt(argv[++i], 10);
    else if (a === '--scenarios') args.scenarios = argv[++i].split(',').map(function (x) { return x.trim(); });
    else if (a === '--sample-sizes') {
      args.sampleSizes = argv[++i].split(',').map(function (x) { return parseInt(x.trim(), 10); });
    } else if (a === '--limit') args.limit = parseInt(argv[++i], 10);
  }
  if (!args.input || !args.outdir) {
    throw new Error('Usage: --input <paired_results.csv> --outdir <dir> [--workers N] [--batch-size N]');
  }
  return args;
}

function ensureDir(d) {
  if (!fs.existsSync(d)) fs.mkdirSync(d, { recursive: true });
}

function splitCsvLine(line) {
  var out = [], cur = '', inQ = false;
  for (var i = 0; i < line.length; i++) {
    var ch = line[i];
    if (inQ) {
      if (ch === '"' && line[i + 1] === '"') { cur += '"'; i++; }
      else if (ch === '"') inQ = false;
      else cur += ch;
    } else if (ch === '"') inQ = true;
    else if (ch === ',') { out.push(cur); cur = ''; }
    else cur += ch;
  }
  out.push(cur);
  return out;
}

function csvEscape(v) {
  if (v === null || v === undefined) return 'NA';
  var s = String(v);
  if (/[",\n]/.test(s)) return '"' + s.replace(/"/g, '""') + '"';
  return s;
}

function scenarioById(id) {
  var s = scenarios.SCENARIOS.find(function (x) { return x.id === id; });
  if (!s) throw new Error('Unknown scenario ' + id);
  return s;
}

function callRBatch(trials, alpha, methods) {
  return new Promise(function (resolve, reject) {
    var tmp = path.join(
      os.tmpdir(),
      'eoss_reexport_' + process.pid + '_' + Date.now() + '_' + Math.random().toString(16).slice(2) + '.json'
    );
    fs.writeFileSync(tmp, JSON.stringify({ alpha: alpha, methods: methods, trials: trials }));
    var child = spawn('Rscript', [R_SCRIPT, '--file=' + tmp], {
      encoding: 'utf8',
      maxBuffer: 64 * 1024 * 1024
    });
    var stdout = '', stderr = '';
    child.stdout.on('data', function (d) { stdout += d; });
    child.stderr.on('data', function (d) { stderr += d; });
    child.on('error', function (err) {
      try { fs.unlinkSync(tmp); } catch (e) { /* ignore */ }
      reject(err);
    });
    child.on('close', function (code) {
      try { fs.unlinkSync(tmp); } catch (e) { /* ignore */ }
      if (code !== 0) {
        reject(new Error('Rscript failed: ' + (stderr || stdout || ('status ' + code))));
        return;
      }
      try {
        resolve(JSON.parse(stdout));
      } catch (e) {
        reject(new Error('Failed to parse R JSON: ' + e.message + '\n' + stdout.slice(0, 500)));
      }
    });
  });
}

function readFilteredRows(input, args) {
  return new Promise(function (resolve, reject) {
    var headers = null;
    var rows = [];
    var rl = readline.createInterface({
      input: fs.createReadStream(input, { encoding: 'utf8' }),
      crlfDelay: Infinity
    });
    rl.on('line', function (line) {
      if (!headers) {
        headers = splitCsvLine(line);
        return;
      }
      if (!line) return;
      var cols = splitCsvLine(line);
      var row = {};
      headers.forEach(function (h, i) { row[h] = cols[i]; });
      if (args.scenarios && args.scenarios.indexOf(row.scenario) < 0) return;
      if (args.sampleSizes && args.sampleSizes.indexOf(Number(row.sample_size)) < 0) return;
      rows.push(row);
      if (args.limit != null && rows.length >= args.limit) {
        rl.pause();
        rl.close();
      }
    });
    rl.on('close', function () { resolve({ headers: headers, rows: rows }); });
    rl.on('error', reject);
  });
}

function boolish(x) {
  var s = String(x).toLowerCase();
  return s === 'true' || s === '1';
}

async function processBatch(batch, alpha, method) {
  var trials = [];
  var rebuilt = [];
  for (var i = 0; i < batch.length; i++) {
    var row = batch[i];
    var sc = scenarioById(row.scenario);
    var n = Number(row.sample_size);
    var seed = Number(row.seed);
    var params = scenarios.scenarioParams(sc, n);
    var nf = params.fus.length;
    var observed = engine.toObservedWide(engine.simulateTrialStructured(params, seed), nf);
    trials.push({ id: String(i), nf: nf, rows: observed });
    rebuilt.push(row);
  }
  var out = await callRBatch(trials, alpha, [method]);
  var byId = {};
  (out.results || []).forEach(function (r) { byId[String(r.id)] = r; });
  var updated = rebuilt.map(function (row, idx) {
    var r = byId[String(idx)] || {};
    var jsOk = boolish(row.js_ok);
    var rOk = r.ok === true || r.ok === 'true';
    var jsEst = jsOk ? Number(row.js_estimate) : null;
    var rEst = rOk && r.estimate != null ? Number(r.estimate) : null;
    var jsSe = jsOk ? Number(row.js_se) : null;
    var rSe = rOk && r.standardError != null ? Number(r.standardError) : null;
    var jsP = jsOk ? Number(row.js_pvalue) : null;
    var rP = rOk && r.pValue != null ? Number(r.pValue) : null;
    var jsReject = jsOk && jsP < alpha;
    var rReject = rOk && rP < alpha;
    var next = Object.assign({}, row);
    next.r_estimate = rEst;
    next.r_se = rSe;
    next.r_pvalue = rP;
    next.r_ok = rOk;
    next.r_convergence = rOk ? !!(r.converged === true || r.converged === 'true') : false;
    next.r_failure_reason = rOk ? 'NA' : (r.failureReason || 'r_fail');
    next.runtime_r = r.runtime_seconds != null ? Number(r.runtime_seconds) : null;
    next.estimate_difference = jsOk && rOk ? jsEst - rEst : null;
    next.se_difference = jsOk && rOk ? jsSe - rSe : null;
    next.js_reject = jsReject;
    next.r_reject = rReject;
    next.decision_disagreement = jsOk && rOk ? (jsReject !== rReject) : null;
    next.r_method = method;
    return next;
  });
  return { rows: updated, session: out.session };
}

async function mapPool(items, concurrency, worker) {
  var out = new Array(items.length);
  var next = 0;
  async function run() {
    while (true) {
      var i = next++;
      if (i >= items.length) return;
      out[i] = await worker(items[i], i);
    }
  }
  var runners = [];
  for (var k = 0; k < concurrency; k++) runners.push(run());
  await Promise.all(runners);
  return out;
}

function writeCsv(file, headers, rows) {
  var lines = [headers.join(',')];
  rows.forEach(function (r) {
    lines.push(headers.map(function (h) { return csvEscape(r[h]); }).join(','));
  });
  fs.writeFileSync(file, lines.join('\n') + '\n');
}

function groupByCell(rows) {
  var map = {};
  rows.forEach(function (r) {
    var key = r.scenario + '|' + r.sample_size;
    if (!map[key]) map[key] = [];
    map[key].push(r);
  });
  return map;
}

function toBoolFields(row) {
  // summarizePaired expects booleans
  return Object.assign({}, row, {
    js_ok: boolish(row.js_ok),
    r_ok: boolish(row.r_ok),
    js_reject: boolish(row.js_reject),
    r_reject: boolish(row.r_reject),
    decision_disagreement: (function () {
      var d = row.decision_disagreement;
      if (d === true || d === false) return d;
      var s = String(d).toLowerCase();
      if (s === 'true' || s === '1') return true;
      if (s === 'false' || s === '0') return false;
      return null;
    })(),
    js_estimate: row.js_estimate == null || row.js_estimate === 'NA' ? null : Number(row.js_estimate),
    r_estimate: row.r_estimate == null || row.r_estimate === 'NA' ? null : Number(row.r_estimate),
    estimate_difference: row.estimate_difference == null || row.estimate_difference === 'NA' ? null : Number(row.estimate_difference),
    js_se: row.js_se == null || row.js_se === 'NA' ? null : Number(row.js_se),
    r_se: row.r_se == null || row.r_se === 'NA' ? null : Number(row.r_se),
    se_difference: row.se_difference == null || row.se_difference === 'NA' ? null : Number(row.se_difference),
    js_pvalue: row.js_pvalue == null || row.js_pvalue === 'NA' ? null : Number(row.js_pvalue),
    r_pvalue: row.r_pvalue == null || row.r_pvalue === 'NA' ? null : Number(row.r_pvalue),
    runtime_js: row.runtime_js == null || row.runtime_js === 'NA' ? null : Number(row.runtime_js),
    runtime_r: row.runtime_r == null || row.runtime_r === 'NA' ? null : Number(row.runtime_r),
    true_full_data_target: Number(row.true_full_data_target),
    observed_data_target: Number(row.observed_data_target),
    sample_size: Number(row.sample_size),
    replicate: Number(row.replicate),
    seed: Number(row.seed),
    js_iterations: row.js_iterations == null || row.js_iterations === 'NA' ? null : Number(row.js_iterations),
    js_convergence: boolish(row.js_convergence),
    r_convergence: boolish(row.r_convergence)
  });
}

async function main() {
  var args = parseArgs(process.argv);
  ensureDir(args.outdir);
  console.log('Reading', args.input);
  var loaded = await readFilteredRows(path.resolve(ROOT, args.input), args);
  console.log('Rows to re-export:', loaded.rows.length, 'workers=', args.workers, 'batch=', args.batchSize);

  // smoke fingerprint on first row
  if (loaded.rows.length) {
    var probe = loaded.rows[0];
    var sc = scenarioById(probe.scenario);
    var params = scenarios.scenarioParams(sc, Number(probe.sample_size));
    var observed = engine.toObservedWide(
      engine.simulateTrialStructured(params, Number(probe.seed)),
      params.fus.length
    );
    var jsProbe = engine.fitMMRM(observed, params.fus.length);
    console.log('Probe JS estimate', jsProbe.estimate, 'stored', probe.js_estimate);
    if (!jsProbe.ok || Math.abs(jsProbe.estimate - Number(probe.js_estimate)) > 1e-10) {
      throw new Error('Seed/generator mismatch on first row; aborting.');
    }
  }

  var batches = [];
  for (var i = 0; i < loaded.rows.length; i += args.batchSize) {
    batches.push(loaded.rows.slice(i, i + args.batchSize));
  }

  var t0 = Date.now();
  var done = 0;
  var rSession = null;
  var updatedChunks = await mapPool(batches, args.workers, async function (batch) {
    var res = await processBatch(batch, 0.05, 'kenward-roger');
    if (!rSession && res.session) rSession = res.session;
    done += batch.length;
    var elapsed = (Date.now() - t0) / 1000;
    var rate = done / Math.max(elapsed, 1e-6);
    var eta = (loaded.rows.length - done) / Math.max(rate, 1e-6);
    if (done % (args.batchSize * args.workers) < args.batchSize || done === loaded.rows.length) {
      console.log(
        'progress', done + '/' + loaded.rows.length,
        'rate=' + rate.toFixed(1) + '/s',
        'eta_min=' + (eta / 60).toFixed(1)
      );
    }
    return res.rows;
  });

  var allRows = [];
  updatedChunks.forEach(function (chunk) {
    chunk.forEach(function (r) { allRows.push(r); });
  });

  // stable order
  allRows.sort(function (a, b) {
    var oa = ['Ex1', 'Ex2', 'Ex3', 'Null'].indexOf(a.scenario);
    var ob = ['Ex1', 'Ex2', 'Ex3', 'Null'].indexOf(b.scenario);
    if (oa !== ob) return oa - ob;
    if (Number(a.sample_size) !== Number(b.sample_size)) return Number(a.sample_size) - Number(b.sample_size);
    return Number(a.replicate) - Number(b.replicate);
  });

  var resultsPath = path.join(args.outdir, 'final_paired_results.csv');
  writeCsv(resultsPath, loaded.headers, allRows);
  console.log('Wrote', resultsPath);

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
      A_js_mean_bias_vs_full: s.A_js_bias_vs_full.mean_bias,
      A_js_mcse: s.A_js_bias_vs_full.mcse,
      A_r_mean_bias_vs_full: s.A_r_bias_vs_full.mean_bias,
      A_r_mcse: s.A_r_bias_vs_full.mcse,
      B_js_mean_bias_vs_observed: s.B_js_bias_vs_observed.mean_bias,
      B_js_mcse: s.B_js_bias_vs_observed.mcse,
      B_r_mean_bias_vs_observed: s.B_r_bias_vs_observed.mean_bias,
      B_r_mcse: s.B_r_bias_vs_observed.mcse,
      js_valid_fits: s.js_valid_fits,
      r_valid_fits: s.r_valid_fits,
      js_failed_fits: s.js_failed_fits,
      r_failed_fits: s.r_failed_fits,
      median_runtime_js: s.median_runtime_js,
      median_runtime_r: s.median_runtime_r,
      median_runtime_ratio_r_over_js: s.median_runtime_ratio_r_over_js
    };
  }

  var byCell = groupByCell(allRows);
  var summaries = [];
  Object.keys(byCell).forEach(function (key) {
    var cellRows = byCell[key].map(toBoolFields);
    var sc = scenarioById(cellRows[0].scenario);
    var summary = summarizePaired(cellRows, {
      scenario: sc.id,
      sample_size: cellRows[0].sample_size,
      n_sim: cellRows.length,
      seed: 1234,
      alpha: 0.05,
      true_full_data_target: sc.true_full_data_target,
      analytic_observed_mixture_target: sc.analytic_observed_mixture_target,
      stored_reference_power_pct: scenarios.storedReferencePowerPct(sc.id, cellRows[0].sample_size),
      r_method: 'kenward-roger'
    });
    summaries.push(flattenSummary(summary));
  });
  summaries.sort(function (a, b) {
    var oa = ['Ex1', 'Ex2', 'Ex3', 'Null'].indexOf(a.scenario);
    var ob = ['Ex1', 'Ex2', 'Ex3', 'Null'].indexOf(b.scenario);
    if (oa !== ob) return oa - ob;
    return a.sample_size - b.sample_size;
  });

  var summaryPath = path.join(args.outdir, 'final_paired_summary.csv');
  writeCsv(summaryPath, Object.keys(summaries[0]), summaries);

  // bias and runtime extracts
  var biasRows = summaries.map(function (s) {
    return {
      run_id: 'fullprec-reexport',
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
      note: 'full-precision R re-export'
    };
  });
  writeCsv(path.join(args.outdir, 'final_bias_summary.csv'), Object.keys(biasRows[0]), biasRows);

  var runtimeRows = summaries.map(function (s) {
    return {
      run_id: 'fullprec-reexport',
      scenario: s.scenario,
      sample_size: s.sample_size,
      median_runtime_js: s.median_runtime_js,
      median_runtime_r: s.median_runtime_r,
      median_runtime_ratio_r_over_js: s.median_runtime_ratio_r_over_js,
      cell_runtime_seconds: null,
      note: 'R times from re-export; JS times retained from original run where present'
    };
  });
  writeCsv(path.join(args.outdir, 'final_runtime_summary.csv'), Object.keys(runtimeRows[0]), runtimeRows);

  var meta = {
    run_id: 'fullprec-reexport',
    source_input: path.relative(ROOT, path.resolve(ROOT, args.input)),
    generated_at: new Date().toISOString(),
    n_rows: allRows.length,
    workers: args.workers,
    batch_size: args.batchSize,
    r_session: rSession,
    elapsed_seconds: (Date.now() - t0) / 1000,
    note: 'R estimates re-exported with jsonlite digits=NA; JS fields retained from source run'
  };
  fs.writeFileSync(path.join(args.outdir, 'run_metadata.json'), JSON.stringify(meta, null, 2) + '\n');
  console.log('Done in', meta.elapsed_seconds.toFixed(1), 's');
}

main().catch(function (err) {
  console.error(err);
  process.exit(1);
});
