/**
 * Worker thread: run a contiguous block of replicates for one scenario/sample size.
 * Each replicate uses predetermined seed = baseSeed + trialIndex * 31337.
 */

'use strict';

var { parentPort, workerData } = require('worker_threads');
var { runOneReplicate, trialSeed } = require('./eoss_engine');

function runBlock(job) {
  var params = job.params;
  var alpha = job.alpha;
  var baseSeed = job.baseSeed;
  var startIndex = job.startIndex;
  var endIndex = job.endIndex; // exclusive
  var results = [];

  for (var i = startIndex; i < endIndex; i++) {
    var seed = trialSeed(baseSeed, i);
    var rep = runOneReplicate(params, seed, alpha);
    results.push({
      index: i,
      seed: seed,
      ok: rep.ok,
      estimate: rep.estimate,
      standardError: rep.standardError,
      pValue: rep.pValue,
      failureReason: rep.failureReason
    });
  }
  return results;
}

if (parentPort) {
  parentPort.on('message', function (job) {
    try {
      var results = runBlock(job);
      parentPort.postMessage({ ok: true, results: results });
    } catch (err) {
      parentPort.postMessage({ ok: false, error: err && err.message ? err.message : String(err) });
    }
  });

  // Support one-shot jobs passed via workerData
  if (workerData && workerData.job) {
    try {
      var results = runBlock(workerData.job);
      parentPort.postMessage({ ok: true, results: results });
    } catch (err) {
      parentPort.postMessage({ ok: false, error: err && err.message ? err.message : String(err) });
    }
  }
}

module.exports = { runBlock: runBlock };
