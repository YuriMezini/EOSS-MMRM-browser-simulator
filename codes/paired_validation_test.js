#!/usr/bin/env node
'use strict';

var assert = require('assert');
var engine = require('../src/eoss_engine');
var scenarios = require('../src/scenarios');
var paired = require('../run_paired_batch');

var failures = [];
function check(name, fn) {
  try {
    fn();
    console.log('  PASS', name);
  } catch (e) {
    failures.push(name + ': ' + e.message);
    console.log('  FAIL', name, e.message);
  }
}

async function main() {
  console.log('Paired validation tests\n');
  var sc = scenarios.getScenario('Ex1');
  var params = scenarios.scenarioParams(sc, 40);

  check('LEGACY and STRUCTURED observed data match for same seed', function () {
    for (var seed = 1234; seed < 1234 + 5 * 31337; seed += 31337) {
      var leg = engine.simTrial(params, seed);
      var obs = engine.toObservedWide(engine.simulateTrialStructured(params, seed), 2);
      assert.strictEqual(leg.length, obs.length);
      for (var i = 0; i < leg.length; i++) {
        assert.strictEqual(leg[i].id, obs[i].id);
        assert.strictEqual(leg[i].trt, obs[i].trt);
        assert.strictEqual(leg[i].y0, obs[i].y0);
        assert.strictEqual(leg[i].y1, obs[i].y1);
        assert.strictEqual(leg[i].y2, obs[i].y2);
      }
    }
  });

  check('Observed mixture targets match closed forms', function () {
    assert.ok(Math.abs(scenarios.getScenario('Ex1').analytic_observed_mixture_target + 5.25) < 1e-12);
    assert.ok(Math.abs(scenarios.getScenario('Ex3').analytic_observed_mixture_target + 4.5) < 1e-12);
  });

  check('Paper sample-size grid has 11 values including 140,180,220', function () {
    assert.deepStrictEqual(scenarios.PAPER_SAMPLE_SIZES,
      [40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240]);
    assert.deepStrictEqual(scenarios.DEFAULT_SAMPLE_SIZES, scenarios.PAPER_SAMPLE_SIZES);
  });

  console.log('  .... running mini paired cell (R required)');
  try {
    var cell = paired.runCell(sc, 40, {
      nSim: 4,
      seed: 1234,
      alpha: 0.05,
      methods: ['kenward-roger'],
      rBatchSize: 4
    });
    assert.strictEqual(cell.rows.length, 4);
    cell.rows.forEach(function (r) {
      assert.ok(r.dataset_fingerprint);
      assert.strictEqual(r.scenario, 'Ex1');
      assert.strictEqual(r.sample_size, 40);
    });
    // Same observed dataset => fingerprint unique per replicate, shared by construction
    var fps = cell.rows.map(function (r) { return r.dataset_fingerprint; });
    assert.strictEqual(new Set(fps).size, fps.length);
    assert.ok(cell.summary.paired_n >= 1, 'expected at least one paired successful fit');
    console.log('  PASS mini paired cell (paired_n=' + cell.summary.paired_n + ')');
  } catch (e) {
    failures.push('mini paired: ' + e.message);
    console.log('  FAIL mini paired', e.message);
  }

  console.log('');
  if (failures.length) {
    console.log(failures.length + ' failed');
    process.exit(1);
  }
  console.log('All paired tests passed.');
}

main();
