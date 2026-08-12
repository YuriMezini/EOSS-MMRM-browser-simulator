/**
 * Paper scenarios with separated target definitions:
 * - true_full_data_target (A): prespecified analytic full-population estimand
 * - analytic_observed_mixture_target (B): closed-form observed-data mixture
 * - paired JS−R difference (C): computed in paired runner, not stored here
 */

'use strict';

var COMMON = {
  baseSeed: 1234,
  alpha: 0.05,
  trtRate: 0.5,
  baselineMean: 38,
  commonSD: 10.5,
  correlation: 0.5,
  numberOfVisits: 2,
  covarianceStructure: 'compound_symmetry'
};

/** Paper Figure 4 grid (40 to 240 by 20). */
var PAPER_SAMPLE_SIZES = [40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240];
var DEFAULT_SAMPLE_SIZES = PAPER_SAMPLE_SIZES.slice();

var STORED_REFERENCE_POWER = {
  Ex1: [36.44, 52.74, 63.26, 73.64, 81.14, 86.66, 90.52, 94.58, 96.06, 97.26, 98.52],
  Ex2: [36.44, 52.74, 67.26, 75.66, 83.52, 89.28, 92.62, 95.48, 97.14, 98.02, 98.92],
  Ex3: [25.52, 37.42, 48.70, 58.10, 65.78, 73.76, 79.22, 84.14, 87.82, 89.28, 92.88],
  Null: null
};

/**
 * Observed-mixture efficacy (final visit, treatment arm), Ex1/Ex3 style:
 *   p_complete = 1 - (etT1+etT2)
 *   p_obs_et   = (etT1+etT2) * efuRate
 *   mean_eff   = (p_complete*eff2 + p_obs_et*effEfu2) / (p_complete+p_obs_et)
 * Signed target = -mean_eff (same sign convention as analytic full targets).
 */
function observedMixtureTarget(v1, v2) {
  var pEt = v1.etT + v2.etT;
  var pComplete = Math.max(0, 1 - pEt);
  var pObsEt = pEt * v2.efuRate;
  var denom = pComplete + pObsEt;
  if (denom <= 0) return null;
  var meanEff = (pComplete * v2.eff + pObsEt * v2.effEfu) / denom;
  return -meanEff;
}

var SCENARIOS = [
  {
    id: 'Ex1',
    name: 'Example 1 — Treatment Policy',
    target_definition_full: 'Prespecified full-population ET-weighted treatment-policy effect',
    target_definition_observed: 'Analytic observed-data mixture at final visit',
    true_full_data_target: -4.8,
    true_effect: -4.8, // legacy alias for older scripts
    analytic_observed_mixture_target: -5.25,
    rejection_label: 'empirical power',
    visits: [
      { eff: 4, effEfu: 2, etT: 0.20, etP: 0.10, efuRate: 0.50 },
      { eff: 6, effEfu: 3, etT: 0.20, etP: 0.10, efuRate: 0.50 }
    ]
  },
  {
    id: 'Ex2',
    name: 'Example 2 — Hypothetical',
    target_definition_full: 'Final-visit hypothetical full-treatment effect',
    target_definition_observed: 'Observed completers under hypothetical (no EFU)',
    true_full_data_target: -6.0,
    true_effect: -6.0,
    analytic_observed_mixture_target: -6.0,
    rejection_label: 'empirical power',
    visits: [
      { eff: 4, effEfu: 0, etT: 0.20, etP: 0.10, efuRate: 0 },
      { eff: 6, effEfu: 0, etT: 0.20, etP: 0.10, efuRate: 0 }
    ]
  },
  {
    id: 'Ex3',
    name: 'Example 3 — Placebo-like controlled',
    target_definition_full: 'Prespecified full-population ET-weighted placebo-like effect',
    target_definition_observed: 'Analytic observed-data mixture at final visit',
    true_full_data_target: -3.6,
    true_effect: -3.6,
    analytic_observed_mixture_target: -4.5,
    rejection_label: 'empirical power',
    visits: [
      { eff: 4, effEfu: 0, etT: 0.20, etP: 0.10, efuRate: 0.50 },
      { eff: 6, effEfu: 0, etT: 0.20, etP: 0.10, efuRate: 0.50 }
    ]
  },
  {
    id: 'Null',
    name: 'Null — Type I error',
    target_definition_full: 'Null treatment effect',
    target_definition_observed: 'Null observed-data target',
    true_full_data_target: 0,
    true_effect: 0,
    analytic_observed_mixture_target: 0,
    rejection_label: 'empirical Type I error',
    visits: [
      { eff: 0, effEfu: 0, etT: 0.20, etP: 0.10, efuRate: 0.50 },
      { eff: 0, effEfu: 0, etT: 0.20, etP: 0.10, efuRate: 0.50 }
    ]
  }
];

// Verify closed forms
SCENARIOS.forEach(function (s) {
  if (s.id === 'Ex1' || s.id === 'Ex3') {
    var calc = observedMixtureTarget(s.visits[0], s.visits[1]);
    if (Math.abs(calc - s.analytic_observed_mixture_target) > 1e-9) {
      s.analytic_observed_mixture_target = calc;
    }
  }
});

var SENSITIVITY_SCENARIOS = [
  {
    id: 'SensCorr02',
    name: 'Sensitivity correlation 0.2 (Ex1 visits)',
    base: 'Ex1',
    corr: 0.2,
    trtRate: 0.5
  },
  {
    id: 'SensCorr08',
    name: 'Sensitivity correlation 0.8 (Ex1 visits)',
    base: 'Ex1',
    corr: 0.8,
    trtRate: 0.5
  },
  {
    id: 'SensAlloc67',
    name: 'Sensitivity allocation 2:1 (Ex1 visits)',
    base: 'Ex1',
    corr: 0.5,
    trtRate: 2 / 3
  }
];

function scenarioParams(scenario, sampleSize, overrides) {
  overrides = overrides || {};
  return {
    ss: sampleSize,
    trtRate: overrides.trtRate != null ? overrides.trtRate : COMMON.trtRate,
    mu: COMMON.baselineMean,
    sd: COMMON.commonSD,
    corr: overrides.corr != null ? overrides.corr : COMMON.correlation,
    fus: scenario.visits.map(function (v) {
      return {
        eff: v.eff,
        effEfu: v.effEfu,
        etT: overrides.etScale != null ? v.etT * overrides.etScale : v.etT,
        etP: overrides.etScale != null ? v.etP * overrides.etScale : v.etP,
        efuRate: overrides.efuRate != null ? overrides.efuRate : v.efuRate
      };
    })
  };
}

function storedReferencePowerPct(scenarioId, sampleSize) {
  var powers = STORED_REFERENCE_POWER[scenarioId];
  if (!powers) return null;
  var idx = PAPER_SAMPLE_SIZES.indexOf(sampleSize);
  if (idx < 0) return null;
  return powers[idx];
}

function getScenario(id) {
  for (var i = 0; i < SCENARIOS.length; i++) {
    if (SCENARIOS[i].id === id) return SCENARIOS[i];
  }
  return null;
}

module.exports = {
  COMMON: COMMON,
  DEFAULT_SAMPLE_SIZES: DEFAULT_SAMPLE_SIZES,
  PAPER_SAMPLE_SIZES: PAPER_SAMPLE_SIZES,
  STORED_REFERENCE_POWER: STORED_REFERENCE_POWER,
  SCENARIOS: SCENARIOS,
  SENSITIVITY_SCENARIOS: SENSITIVITY_SCENARIOS,
  scenarioParams: scenarioParams,
  storedReferencePowerPct: storedReferencePowerPct,
  getScenario: getScenario,
  observedMixtureTarget: observedMixtureTarget
};
