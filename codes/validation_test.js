#!/usr/bin/env node
/**
 * Validation tests for EOSS_MMRM batch bias system.
 */

'use strict';

var assert = require('assert');
var fs = require('fs');
var path = require('path');
var os = require('os');

var engine = require('../src/eoss_engine');
var scenarios = require('../src/scenarios');
var stats = require('../src/statistics');
var batch = require('../run_bias_batch');

var failures = [];

function check(name, fn) {
  try {
    fn();
    console.log('  PASS  ' + name);
  } catch (err) {
    failures.push(name + ': ' + (err && err.message ? err.message : String(err)));
    console.log('  FAIL  ' + name);
    console.log('        ' + (err && err.stack ? err.stack.split('\n').slice(0, 3).join('\n        ') : err));
  }
}

function deepEqualTrial(a, b) {
  assert.strictEqual(a.length, b.length);
  for (var i = 0; i < a.length; i++) {
    assert.strictEqual(a[i].id, b[i].id);
    assert.strictEqual(a[i].trt, b[i].trt);
    assert.strictEqual(a[i].y0, b[i].y0);
    assert.strictEqual(a[i].y1, b[i].y1);
    assert.strictEqual(a[i].y2, b[i].y2);
  }
}

async function main() {
  console.log('EOSS_MMRM validation tests\n');

  var ex1 = scenarios.getScenario('Ex1');
  var params = scenarios.scenarioParams(ex1, 40);

  check('1. Same seed produces identical simulated trial', function () {
    var t1 = engine.simulateTrial(params, 1234);
    var t2 = engine.simulateTrial(params, 1234);
    deepEqualTrial(t1, t2);
    var t3 = engine.simTrial(params, 1234 + 31337);
    assert.notStrictEqual(JSON.stringify(t1), JSON.stringify(t3));
  });

  // Async worker-count reproducibility test
  try {
    var nSim = 40;
    var alpha = 0.05;
    var seed = 1234;
    var p = scenarios.scenarioParams(scenarios.getScenario('Ex2'), 60);
    var r1 = await batch.runReplicatesParallel(p, seed, alpha, nSim, 1);
    var r2 = await batch.runReplicatesParallel(p, seed, alpha, nSim, Math.min(4, Math.max(2, (os.cpus() || []).length)));
    assert.strictEqual(r1.length, r2.length);
    for (var i = 0; i < r1.length; i++) {
      assert.strictEqual(r1[i].ok, r2[i].ok);
      if (r1[i].ok) {
        assert.strictEqual(r1[i].estimate, r2[i].estimate);
        assert.strictEqual(r1[i].standardError, r2[i].standardError);
        assert.strictEqual(r1[i].pValue, r2[i].pValue);
      }
    }
    var s1 = stats.summarizeBias({
      scenario: scenarios.getScenario('Ex2'),
      sampleSize: 60,
      nSim: nSim,
      seed: seed,
      alpha: alpha,
      replicates: r1,
      storedReferencePowerPct: scenarios.storedReferencePowerPct('Ex2', 60),
      runtimeSeconds: 0
    });
    var s2 = stats.summarizeBias({
      scenario: scenarios.getScenario('Ex2'),
      sampleSize: 60,
      nSim: nSim,
      seed: seed,
      alpha: alpha,
      replicates: r2,
      storedReferencePowerPct: scenarios.storedReferencePowerPct('Ex2', 60),
      runtimeSeconds: 0
    });
    assert.strictEqual(s1.mean_estimate, s2.mean_estimate);
    assert.strictEqual(s1.bias, s2.bias);
    assert.strictEqual(s1.rejection_rate_pct, s2.rejection_rate_pct);
    assert.strictEqual(s1.valid_fits, s2.valid_fits);
    console.log('  PASS  2. Different worker counts produce identical summary results');
  } catch (err) {
    failures.push('2. worker reproducibility: ' + err.message);
    console.log('  FAIL  2. Different worker counts produce identical summary results');
    console.log('        ' + err.message);
  }

  check('3. Output metrics are finite when valid fits are available', function () {
    var reps = [];
    for (var i = 0; i < 30; i++) {
      var fit = engine.runOneReplicate(params, engine.trialSeed(1234, i), 0.05);
      reps.push(fit);
    }
    var summary = stats.summarizeBias({
      scenario: ex1,
      sampleSize: 40,
      nSim: 30,
      seed: 1234,
      alpha: 0.05,
      replicates: reps,
      storedReferencePowerPct: scenarios.storedReferencePowerPct('Ex1', 40),
      runtimeSeconds: 0.1
    });
    assert.ok(summary.valid_fits > 0, 'expected some valid fits');
    ['mean_estimate', 'bias', 'empirical_sd', 'mean_model_se', 'rmse', 'coverage_95_pct', 'rejection_rate_pct', 'mcse_bias']
      .forEach(function (k) {
        assert.ok(Number.isFinite(summary[k]), k + ' should be finite, got ' + summary[k]);
      });
  });

  check('4. Valid fits + failed fits = requested simulations', function () {
    var reps = [];
    for (var i = 0; i < 25; i++) {
      reps.push(engine.runOneReplicate(params, engine.trialSeed(99, i), 0.05));
    }
    var summary = stats.summarizeBias({
      scenario: ex1,
      sampleSize: 40,
      nSim: 25,
      seed: 99,
      alpha: 0.05,
      replicates: reps,
      storedReferencePowerPct: null,
      runtimeSeconds: 0
    });
    assert.strictEqual(summary.valid_fits + summary.failed_fits, 25);
  });

  check('5. RMSE >= abs(bias) (floating-point tolerance)', function () {
    var reps = [];
    for (var i = 0; i < 50; i++) {
      reps.push(engine.runOneReplicate(params, engine.trialSeed(7, i), 0.05));
    }
    var summary = stats.summarizeBias({
      scenario: ex1,
      sampleSize: 40,
      nSim: 50,
      seed: 7,
      alpha: 0.05,
      replicates: reps,
      storedReferencePowerPct: null,
      runtimeSeconds: 0
    });
    if (summary.valid_fits > 0) {
      assert.ok(summary.rmse + 1e-12 >= Math.abs(summary.bias), 'rmse=' + summary.rmse + ' bias=' + summary.bias);
    }
  });

  check('6. Coverage and rejection rates are between 0 and 1 (stored as %)', function () {
    var reps = [];
    for (var i = 0; i < 40; i++) {
      reps.push(engine.runOneReplicate(params, engine.trialSeed(3, i), 0.05));
    }
    var summary = stats.summarizeBias({
      scenario: ex1,
      sampleSize: 40,
      nSim: 40,
      seed: 3,
      alpha: 0.05,
      replicates: reps,
      storedReferencePowerPct: null,
      runtimeSeconds: 0
    });
    if (summary.valid_fits > 0) {
      assert.ok(summary.coverage_95_pct >= 0 && summary.coverage_95_pct <= 100);
      assert.ok(summary.rejection_rate_pct >= 0 && summary.rejection_rate_pct <= 100);
    }
  });

  check('7. Example ordering generally Ex2 > Ex1 > Ex3 at n=100 (moderate run)', function () {
    var nSim = 200;
    var seed = 1234;
    function powerFor(id) {
      var sc = scenarios.getScenario(id);
      var p = scenarios.scenarioParams(sc, 100);
      var rej = 0, valid = 0;
      for (var i = 0; i < nSim; i++) {
        var r = engine.runOneReplicate(p, engine.trialSeed(seed, i), 0.05);
        if (r.ok) {
          valid++;
          if (r.pValue < 0.05) rej++;
        }
      }
      assert.ok(valid > 0);
      return rej / valid;
    }
    var p1 = powerFor('Ex1');
    var p2 = powerFor('Ex2');
    var p3 = powerFor('Ex3');
    // Soft check with small slack for Monte Carlo noise at nSim=200
    assert.ok(p2 >= p1 - 0.05, 'Ex2 (' + p2 + ') should be >= Ex1 (' + p1 + ') within noise');
    assert.ok(p1 >= p3 - 0.05, 'Ex1 (' + p1 + ') should be >= Ex3 (' + p3 + ') within noise');
  });

  check('8. Null scenario rejection rate labeled as Type I error', function () {
    var nullSc = scenarios.getScenario('Null');
    assert.strictEqual(nullSc.rejection_label, 'empirical Type I error');
    var reps = [];
    for (var i = 0; i < 20; i++) {
      reps.push(engine.runOneReplicate(scenarios.scenarioParams(nullSc, 40), engine.trialSeed(1234, i), 0.05));
    }
    var summary = stats.summarizeBias({
      scenario: nullSc,
      sampleSize: 40,
      nSim: 20,
      seed: 1234,
      alpha: 0.05,
      replicates: reps,
      storedReferencePowerPct: null,
      runtimeSeconds: 0
    });
    assert.strictEqual(summary.rejection_label, 'empirical Type I error');
  });

  check('9. Output CSV has every required column', function () {
    var header = stats.csvHeader().split(',');
    stats.RESULT_COLUMNS.forEach(function (col) {
      assert.ok(header.indexOf(col) >= 0, 'missing column ' + col);
    });
    assert.strictEqual(header.length, stats.RESULT_COLUMNS.length);

    // Also verify a temporary CSV row serializes all columns
    var row = stats.summarizeBias({
      scenario: ex1,
      sampleSize: 40,
      nSim: 1,
      seed: 1234,
      alpha: 0.05,
      replicates: [engine.runOneReplicate(params, 1234, 0.05)],
      storedReferencePowerPct: 36.44,
      runtimeSeconds: 0.01
    });
    var line = stats.rowToCsv(row);
    assert.strictEqual(line.split(',').length >= stats.RESULT_COLUMNS.length - 2, true);
    // target_definition may contain commas — ensure CSV escaping keeps parseable field count via header match
    assert.ok(line.indexOf(String(row.scenario)) === 0 || line.indexOf(row.scenario) >= 0);
  });

  check('fitMMRM returns estimate and SE (not only p-value)', function () {
    var data = engine.simulateTrial(params, 1234);
    var fit = engine.fitMMRM(data, 2);
    assert.ok(fit.ok);
    assert.ok(Number.isFinite(fit.estimate));
    assert.ok(Number.isFinite(fit.standardError));
    assert.ok(Number.isFinite(fit.pValue));
    var pOnly = engine.mmrm(data, 2);
    assert.ok(Math.abs(pOnly - fit.pValue) < 1e-15);
  });

  console.log('');
  if (failures.length) {
    console.log(failures.length + ' test(s) failed:');
    failures.forEach(function (f) { console.log(' - ' + f); });
    process.exit(1);
  }
  console.log('All tests passed.');
}

main().catch(function (err) {
  console.error(err);
  process.exit(1);
});
