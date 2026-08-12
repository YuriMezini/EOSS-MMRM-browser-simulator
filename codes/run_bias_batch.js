#!/usr/bin/env node
/**
 * EOSS_MMRM batch model-bias validation runner.
 *
 * Usage:
 *   node run_bias_batch.js --n-sim 100 --sample-sizes 40,100 --workers 4 --seed 1234
 *
 * npm scripts:
 *   npm run smoke | preliminary | final
 */

'use strict';

var fs = require('fs');
var path = require('path');
var os = require('os');
var { Worker } = require('worker_threads');
var XLSX = require('xlsx');

var {
  COMMON,
  DEFAULT_SAMPLE_SIZES,
  PAPER_SAMPLE_SIZES,
  STORED_REFERENCE_POWER,
  SCENARIOS,
  scenarioParams,
  storedReferencePowerPct
} = require('./src/scenarios');
var { summarizeBias, rowToCsv, csvHeader, RESULT_COLUMNS } = require('./src/statistics');

var ROOT = __dirname;
var RESULTS_DIR = path.join(ROOT, 'results');
var CSV_PATH = path.join(RESULTS_DIR, 'model_bias_results.csv');
var JSON_PATH = path.join(RESULTS_DIR, 'model_bias_results.json');
var XLSX_PATH = path.join(RESULTS_DIR, 'model_bias_results.xlsx');
var META_PATH = path.join(RESULTS_DIR, 'run_metadata.json');
var WORKER_PATH = path.join(ROOT, 'src', 'bias_worker.js');

function parseArgs(argv) {
  var args = {
    nSim: 100,
    workers: Math.max(1, (os.cpus() || []).length - 1),
    seed: COMMON.baseSeed,
    alpha: COMMON.alpha,
    sampleSizes: DEFAULT_SAMPLE_SIZES.slice(),
    scenarios: SCENARIOS.map(function (s) { return s.id; }),
    mode: 'custom'
  };

  for (var i = 2; i < argv.length; i++) {
    var a = argv[i];
    if (a === '--n-sim' && argv[i + 1]) args.nSim = parseInt(argv[++i], 10);
    else if (a === '--workers' && argv[i + 1]) args.workers = Math.max(1, parseInt(argv[++i], 10));
    else if (a === '--seed' && argv[i + 1]) args.seed = parseInt(argv[++i], 10);
    else if (a === '--alpha' && argv[i + 1]) args.alpha = parseFloat(argv[++i]);
    else if (a === '--sample-sizes' && argv[i + 1]) {
      args.sampleSizes = argv[++i].split(',').map(function (x) { return parseInt(x.trim(), 10); }).filter(function (n) { return n > 0; });
    } else if (a === '--scenarios' && argv[i + 1]) {
      args.scenarios = argv[++i].split(',').map(function (x) { return x.trim(); });
    } else if (a === '--mode' && argv[i + 1]) args.mode = argv[++i];
    else if (a === '--help' || a === '-h') {
      printHelp();
      process.exit(0);
    }
  }
  return args;
}

function printHelp() {
  console.log([
    'EOSS_MMRM batch model-bias validation',
    '',
    'Options:',
    '  --n-sim N              Simulations per scenario/sample-size (default 100)',
    '  --workers N            Worker threads (default CPU cores - 1)',
    '  --seed N               Base seed (default 1234)',
    '  --alpha A              Significance level (default 0.05)',
    '  --sample-sizes LIST    Comma-separated sample sizes',
    '  --scenarios LIST       Comma-separated scenario ids (Ex1,Ex2,Ex3,Null)',
    '  --mode NAME            Label written to metadata (smoke|preliminary|final|custom)',
    '',
    'Examples:',
    '  node run_bias_batch.js --n-sim 5000 --workers 6 --seed 1234',
    '  node run_bias_batch.js --n-sim 100 --sample-sizes 40,100 --workers 2'
  ].join('\n'));
}

function ensureResultsDir() {
  if (!fs.existsSync(RESULTS_DIR)) fs.mkdirSync(RESULTS_DIR, { recursive: true });
}

function formatPct(x, digits) {
  if (x === null || x === undefined || !Number.isFinite(x)) return 'NA';
  return x.toFixed(digits == null ? 1 : digits) + '%';
}

function formatNum(x, digits) {
  if (x === null || x === undefined || !Number.isFinite(x)) return 'NA';
  return x.toFixed(digits == null ? 4 : digits);
}

function printSmokeTable(rows) {
  console.log('');
  console.log('Smoke / summary table');
  console.log([
    'Scenario'.padEnd(8),
    'SS'.padStart(5),
    'True'.padStart(8),
    'MeanEst'.padStart(9),
    'Bias'.padStart(9),
    'EmpSD'.padStart(8),
    'MeanSE'.padStart(8),
    'SE/SD'.padStart(7),
    'RMSE'.padStart(8),
    'Cov%'.padStart(7),
    'Rej%'.padStart(7),
    'Fail%'.padStart(7)
  ].join(' '));
  rows.forEach(function (r) {
    console.log([
      String(r.scenario).padEnd(8),
      String(r.sample_size).padStart(5),
      formatNum(r.true_effect, 3).padStart(8),
      formatNum(r.mean_estimate, 3).padStart(9),
      formatNum(r.bias, 3).padStart(9),
      formatNum(r.empirical_sd, 3).padStart(8),
      formatNum(r.mean_model_se, 3).padStart(8),
      formatNum(r.se_sd_ratio, 3).padStart(7),
      formatNum(r.rmse, 3).padStart(8),
      formatPct(r.coverage_95_pct, 1).padStart(7),
      formatPct(r.rejection_rate_pct, 1).padStart(7),
      formatPct(r.failure_rate_pct, 1).padStart(7)
    ].join(' '));
  });
  console.log('');
}

/**
 * Partition [0, nSim) into up to nWorkers contiguous blocks.
 * Contiguous blocks + predetermined seeds => worker count cannot change results.
 */
function partitionIndices(nSim, nWorkers) {
  var workers = Math.min(Math.max(1, nWorkers), nSim);
  var base = Math.floor(nSim / workers);
  var rem = nSim % workers;
  var parts = [];
  var start = 0;
  for (var w = 0; w < workers; w++) {
    var len = base + (w < rem ? 1 : 0);
    if (len <= 0) continue;
    parts.push({ startIndex: start, endIndex: start + len });
    start += len;
  }
  return parts;
}

function runReplicatesParallel(params, baseSeed, alpha, nSim, nWorkers) {
  var parts = partitionIndices(nSim, nWorkers);
  if (parts.length === 1) {
    // Inline path avoids worker overhead for tiny jobs / tests
    var { runBlock } = require('./src/bias_worker');
    return Promise.resolve(runBlock({
      params: params,
      alpha: alpha,
      baseSeed: baseSeed,
      startIndex: parts[0].startIndex,
      endIndex: parts[0].endIndex
    }));
  }

  return Promise.all(parts.map(function (part) {
    return new Promise(function (resolve, reject) {
      var worker = new Worker(WORKER_PATH);
      var settled = false;
      worker.on('message', function (msg) {
        settled = true;
        worker.terminate();
        if (!msg.ok) reject(new Error(msg.error || 'worker failed'));
        else resolve(msg.results);
      });
      worker.on('error', function (err) {
        if (!settled) {
          settled = true;
          reject(err);
        }
      });
      worker.on('exit', function (code) {
        if (!settled && code !== 0) reject(new Error('worker exited with code ' + code));
      });
      worker.postMessage({
        params: params,
        alpha: alpha,
        baseSeed: baseSeed,
        startIndex: part.startIndex,
        endIndex: part.endIndex
      });
    });
  })).then(function (chunks) {
    var all = [];
    chunks.forEach(function (chunk) {
      for (var i = 0; i < chunk.length; i++) all.push(chunk[i]);
    });
    all.sort(function (a, b) { return a.index - b.index; });
    return all;
  });
}

function writeIncrementalOutputs(rows) {
  ensureResultsDir();
  var csv = csvHeader() + '\n' + rows.map(rowToCsv).join('\n') + '\n';
  fs.writeFileSync(CSV_PATH, csv, 'utf8');
  fs.writeFileSync(JSON_PATH, JSON.stringify(rows, null, 2), 'utf8');
}

function highlightRow(ws, rowIndex, rowObj) {
  // xlsx (SheetJS community) has limited styling; store flags as notes in a helper sheet instead.
  // We still freeze panes / autofilter below. Conditional visual cues are documented in README.
  void ws; void rowIndex; void rowObj;
}

function buildExcel(rows, meta) {
  var wb = XLSX.utils.book_new();

  // Sheet 1: Model Bias Results
  var resultAoA = [RESULT_COLUMNS];
  rows.forEach(function (r) {
    resultAoA.push(RESULT_COLUMNS.map(function (c) {
      var v = r[c];
      return v === null || v === undefined ? 'NA' : v;
    }));
  });
  var wsResults = XLSX.utils.aoa_to_sheet(resultAoA);
  wsResults['!autofilter'] = { ref: XLSX.utils.encode_range({ s: { r: 0, c: 0 }, e: { r: rows.length, c: RESULT_COLUMNS.length - 1 } }) };
  wsResults['!freeze'] = { xSplit: 0, ySplit: 1 };
  if (!wsResults['!views']) wsResults['!views'] = [{ state: 'frozen', ySplit: 1, topLeftCell: 'A2', activeCell: 'A2' }];
  XLSX.utils.book_append_sheet(wb, wsResults, 'Model Bias Results');

  // Sheet 2: Scenario Definitions
  var scenAoA = [[
    'scenario', 'name', 'target_definition', 'true_effect', 'rejection_label',
    'visit', 'eff', 'effEfu', 'etT', 'etP', 'efuRate',
    'trtRate', 'baselineMean', 'commonSD', 'correlation', 'numberOfVisits', 'baseSeed', 'alpha'
  ]];
  SCENARIOS.forEach(function (s) {
    s.visits.forEach(function (v, idx) {
      scenAoA.push([
        s.id, s.name, s.target_definition, s.true_effect, s.rejection_label,
        idx + 1, v.eff, v.effEfu, v.etT, v.etP, v.efuRate,
        COMMON.trtRate, COMMON.baselineMean, COMMON.commonSD, COMMON.correlation,
        COMMON.numberOfVisits, COMMON.baseSeed, COMMON.alpha
      ]);
    });
  });
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(scenAoA), 'Scenario Definitions');

  // Sheet 3: Paper References
  var paperAoA = [['scenario', 'sample_size', 'stored_reference_power_pct']];
  ['Ex1', 'Ex2', 'Ex3'].forEach(function (id) {
    PAPER_SAMPLE_SIZES.forEach(function (ss, i) {
      paperAoA.push([id, ss, STORED_REFERENCE_POWER[id][i]]);
    });
  });
  paperAoA.push([]);
  paperAoA.push(['Note', 'Power difference versus R/mmrm reference = HTML rejection_rate_pct - stored_reference_power_pct']);
  paperAoA.push(['Note', 'Do not label power difference as estimator bias.']);
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(paperAoA), 'Paper References');

  // Sheet 4: Run Metadata + highlight guidance
  var metaAoA = [['key', 'value']];
  Object.keys(meta).forEach(function (k) {
    var v = meta[k];
    metaAoA.push([k, typeof v === 'object' ? JSON.stringify(v) : v]);
  });
  metaAoA.push([]);
  metaAoA.push(['highlight_rule', 'abs(relative_bias_pct) > 5']);
  metaAoA.push(['highlight_rule', 'coverage_95_pct < 93 or > 97']);
  metaAoA.push(['highlight_rule', 'se_sd_ratio < 0.90 or > 1.10']);
  metaAoA.push(['highlight_rule', 'Null Type I error < 4 or > 6']);
  metaAoA.push(['highlight_rule', 'failure_rate_pct > 1']);
  metaAoA.push(['highlight_rule', 'abs(power_difference_pp) > 2']);
  metaAoA.push(['formatting_note', 'Community SheetJS has limited cell styling; apply conditional formatting in Excel using the rules above.']);
  metaAoA.push(['formatting_note', 'Percentages are stored as numeric percent points (e.g. 73.5 means 73.5%). Effects use 3–4 decimal places in CSV/JSON.']);
  XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(metaAoA), 'Run Metadata');

  XLSX.writeFile(wb, XLSX_PATH);
}

async function runBatch(args) {
  ensureResultsDir();
  var selected = SCENARIOS.filter(function (s) { return args.scenarios.indexOf(s.id) >= 0; });
  if (!selected.length) throw new Error('No valid scenarios selected');
  if (!args.sampleSizes.length) throw new Error('No sample sizes provided');

  var combinations = [];
  selected.forEach(function (sc) {
    args.sampleSizes.forEach(function (ss) {
      combinations.push({ scenario: sc, sampleSize: ss });
    });
  });

  var totalCombos = combinations.length;
  var rows = [];
  var t0 = Date.now();

  console.log('EOSS_MMRM model-bias batch validation');
  console.log('Mode:          ' + args.mode);
  console.log('Simulations:   ' + args.nSim + ' per combination');
  console.log('Scenarios:     ' + selected.map(function (s) { return s.id; }).join(', '));
  console.log('Sample sizes:  ' + args.sampleSizes.join(', '));
  console.log('Workers:       ' + args.workers);
  console.log('Base seed:     ' + args.seed);
  console.log('Combinations:  ' + totalCombos);
  console.log('Total fits:    ' + (totalCombos * args.nSim));
  console.log('');

  for (var ci = 0; ci < combinations.length; ci++) {
    var combo = combinations[ci];
    var params = scenarioParams(combo.scenario, combo.sampleSize);
    var started = Date.now();
    var replicates = await runReplicatesParallel(params, args.seed, args.alpha, args.nSim, args.workers);
    var runtimeSeconds = (Date.now() - started) / 1000;
    var row = summarizeBias({
      scenario: combo.scenario,
      sampleSize: combo.sampleSize,
      nSim: args.nSim,
      seed: args.seed,
      alpha: args.alpha,
      replicates: replicates,
      storedReferencePowerPct: storedReferencePowerPct(combo.scenario.id, combo.sampleSize),
      runtimeSeconds: runtimeSeconds
    });
    rows.push(row);
    writeIncrementalOutputs(rows);

    var done = ci + 1;
    var remaining = totalCombos - done;
    var elapsed = ((Date.now() - t0) / 1000).toFixed(1);
    console.log(
      '[' + done + '/' + totalCombos + '] ' +
      combo.scenario.id + ' n=' + combo.sampleSize +
      '  valid=' + row.valid_fits + '/' + row.n_sim_requested +
      '  bias=' + formatNum(row.bias, 4) +
      '  rej=' + formatPct(row.rejection_rate_pct, 1) +
      '  (' + runtimeSeconds.toFixed(1) + 's)  elapsed=' + elapsed + 's  remaining_combos=' + remaining
    );
  }

  var meta = {
    mode: args.mode,
    n_sim: args.nSim,
    workers: args.workers,
    seed: args.seed,
    alpha: args.alpha,
    sample_sizes: args.sampleSizes,
    scenarios: selected.map(function (s) { return s.id; }),
    combinations: totalCombos,
    total_fits_requested: totalCombos * args.nSim,
    node_version: process.version,
    platform: process.platform + ' ' + os.release(),
    cpus: (os.cpus() || []).length,
    started_at: new Date(t0).toISOString(),
    finished_at: new Date().toISOString(),
    elapsed_seconds: (Date.now() - t0) / 1000,
    engine: 'HTML iterative GLS approximation with EM-type covariance updates',
    reference_method: 'R mmrm with Kenward-Roger (paper Figure 4)',
    seed_rule: 'trialSeed = baseSeed + trialIndex * 31337',
    output_files: {
      csv: CSV_PATH,
      json: JSON_PATH,
      xlsx: XLSX_PATH,
      metadata: META_PATH
    }
  };

  fs.writeFileSync(META_PATH, JSON.stringify(meta, null, 2), 'utf8');
  buildExcel(rows, meta);
  writeIncrementalOutputs(rows);

  printSmokeTable(rows);
  console.log('Wrote:');
  console.log('  ' + CSV_PATH);
  console.log('  ' + JSON_PATH);
  console.log('  ' + XLSX_PATH);
  console.log('  ' + META_PATH);
  console.log('Total elapsed: ' + meta.elapsed_seconds.toFixed(1) + 's');

  return { rows: rows, meta: meta };
}

if (require.main === module) {
  var args = parseArgs(process.argv);
  runBatch(args).catch(function (err) {
    console.error('Batch validation failed:', err);
    process.exit(1);
  });
}

module.exports = {
  parseArgs: parseArgs,
  runBatch: runBatch,
  runReplicatesParallel: runReplicatesParallel,
  partitionIndices: partitionIndices,
  CSV_PATH: CSV_PATH,
  RESULT_COLUMNS: RESULT_COLUMNS
};
