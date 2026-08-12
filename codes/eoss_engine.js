/**
 * EOSS_MMRM statistical engine
 * Extracted from eoss_mmrm_simulator_bias_validation.html (Web Worker WK_SRC).
 *
 * Modes:
 * - LEGACY: simTrial() — observed rows only; RNG path unchanged from HTML.
 * - STRUCTURED: simulateTrialStructured() — latent outcomes, observation
 *   indicators, ET/EFU flags; observed outcomes match LEGACY bit-for-bit.
 *
 * Fitting is an iterative GLS approximation with EM-type covariance updates
 * (not established as REML). Terminology: "EM-type GLS".
 */

'use strict';

function mkR(s) {
  return function () {
    s = (s + 0x6d2b79f5) | 0;
    var t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function rn(r) {
  var u = 0, v = 0;
  while (!u) u = r();
  while (!v) v = r();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(6.28318530718 * v);
}

function chol(d, sd, co) {
  var L = [], i, j, k, s;
  for (i = 0; i < d; i++) {
    L.push(new Float64Array(d));
    for (j = 0; j <= i; j++) {
      s = i === j ? sd * sd : co * sd * sd;
      for (k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      L[i][j] = i === j ? Math.sqrt(Math.max(s, 0)) : s / L[j][j];
    }
  }
  return L;
}

function mvr(ss, d, mu, sd, co, rng) {
  var L = chol(d, sd, co), r = [], i, j, row, c, z, x;
  for (i = 0; i < ss; i++) {
    z = [];
    for (j = 0; j < d; j++) z.push(rn(rng));
    x = new Float64Array(d);
    for (row = 0; row < d; row++) {
      x[row] = mu;
      for (c = 0; c <= row; c++) x[row] += L[row][c] * z[c];
    }
    r.push(x);
  }
  return r;
}

function scat(p, rng) {
  var u = rng(), c = 0, i;
  for (i = 0; i < p.length; i++) {
    c += p[i];
    if (u < c) return i;
  }
  return p.length - 1;
}

/**
 * Shared generator. mode: 'legacy' | 'structured'
 * RNG consumption order is identical to the original HTML simTrial.
 */
function generateTrial(p, seed, mode) {
  var structured = mode === 'structured';
  var ss = p.ss, nf = p.fus.length, rng = mkR(seed), i, fi, fu, isET, efuDraw, efuOn;
  var dat = mvr(ss, nf + 1, p.mu, p.sd, p.corr, rng);
  var trt = new Uint8Array(ss);
  for (i = 0; i < ss; i++) trt[i] = rng() < p.trtRate ? 1 : 0;
  var tp0 = Math.max(0, 1 - p.fus.reduce(function (a, f) { return a + f.etT; }, 0));
  var pp0 = Math.max(0, 1 - p.fus.reduce(function (a, f) { return a + f.etP; }, 0));
  var tp = [tp0].concat(p.fus.map(function (f) { return f.etT; }));
  var pp = [pp0].concat(p.fus.map(function (f) { return f.etP; }));
  var tS = tp.reduce(function (a, b) { return a + b; }, 0);
  var pS = pp.reduce(function (a, b) { return a + b; }, 0);
  for (i = 0; i < tp.length; i++) {
    tp[i] /= tS;
    pp[i] /= pS;
  }
  var et = new Int32Array(ss);
  for (i = 0; i < ss; i++) et[i] = scat(trt[i] ? tp : pp, rng);
  var rows = [];
  for (i = 0; i < ss; i++) {
    var y0 = Math.round(dat[i][0] * 100) / 100;
    var row = structured
      ? { id: i + 1, trt: trt[i], et: et[i], y0: y0 }
      : { id: i + 1, trt: trt[i], y0: y0 };
    for (fi = 0; fi < nf; fi++) {
      fu = p.fus[fi];
      isET = et[i] >= 1 && et[i] <= fi + 1;
      var v = fi + 1;
      var appliedEff = null;
      var R = 0;
      efuOn = 0;
      if (!isET) {
        appliedEff = fu.eff;
        R = 1;
      } else {
        // Same Bernoulli draw as legacy HTML (efuRate decides observed EFU)
        efuDraw = rng() < fu.efuRate;
        efuOn = efuDraw ? 1 : 0;
        if (efuDraw) {
          appliedEff = fu.effEfu;
          R = 1;
        } else {
          appliedEff = null; // legacy: missing
          R = 0;
        }
      }
      var yObs = null;
      if (appliedEff !== null) {
        yObs = Math.round((dat[i][fi + 1] - appliedEff * trt[i]) * 100) / 100;
      }
      if (!structured) {
        row['y' + v] = yObs;
      } else {
        // Full-data latent under ET uses effEfu even when unobserved (for A targets)
        var latentEff = isET ? fu.effEfu : fu.eff;
        var yLat = Math.round((dat[i][fi + 1] - latentEff * trt[i]) * 100) / 100;
        row['y_lat' + v] = yLat;
        row['R' + v] = R;
        row['efu' + v] = efuOn;
        // Observed column always follows LEGACY rule (yObs)
        row['y' + v] = yObs;
      }
    }
    rows.push(row);
  }
  return rows;
}

/** LEGACY: observed wide trial (HTML-compatible). */
function simTrial(p, seed) {
  return generateTrial(p, seed, 'legacy');
}

/**
 * STRUCTURED trial: latent, R, observed, ET, EFU.
 * Observed columns y1..yf match LEGACY for the same seed.
 */
function simulateTrialStructured(p, seed) {
  return generateTrial(p, seed, 'structured');
}

/** Extract HTML/R-compatible observed data frame from structured rows. */
function toObservedWide(structuredRows, nf) {
  return structuredRows.map(function (r) {
    var o = { id: r.id, trt: r.trt, y0: r.y0 };
    for (var v = 1; v <= nf; v++) o['y' + v] = r['y' + v];
    return o;
  });
}

/**
 * Empirical full-data contrast at final visit from latent outcomes:
 * mean(y_lat_final | trt=1) - mean(y_lat_final | trt=0)
 * (same sign convention as the fitted contrast: treated minus placebo on change-like scale)
 * Note: the analysis contrast is a model estimand; this is a descriptive full-data mean difference.
 */
function empiricalFullDataTarget(structuredRows, nf) {
  var t = [], c = [];
  for (var i = 0; i < structuredRows.length; i++) {
    var y = structuredRows[i]['y_lat' + nf];
    if (y == null || !Number.isFinite(y)) continue;
    if (structuredRows[i].trt === 1) t.push(y);
    else c.push(y);
  }
  if (!t.length || !c.length) return null;
  var mt = t.reduce(function (a, b) { return a + b; }, 0) / t.length;
  var mc = c.reduce(function (a, b) { return a + b; }, 0) / c.length;
  return mt - mc;
}

/** Observed-data mean difference at final visit among subjects with R_final=1. */
function empiricalObservedDataTarget(structuredRows, nf) {
  var t = [], c = [];
  for (var i = 0; i < structuredRows.length; i++) {
    if (structuredRows[i]['R' + nf] !== 1) continue;
    var y = structuredRows[i]['y' + nf];
    if (y == null || !Number.isFinite(y)) continue;
    if (structuredRows[i].trt === 1) t.push(y);
    else c.push(y);
  }
  if (!t.length || !c.length) return null;
  var mt = t.reduce(function (a, b) { return a + b; }, 0) / t.length;
  var mc = c.reduce(function (a, b) { return a + b; }, 0) / c.length;
  return mt - mc;
}

function matInv(A, n) {
  var M = new Float64Array(n * 2 * n), i, j, row, piv, d, f, t, R;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) M[i * (2 * n) + j] = A[i * n + j];
    M[i * (2 * n) + n + i] = 1;
  }
  for (var col = 0; col < n; col++) {
    piv = col;
    for (row = col + 1; row < n; row++) {
      if (Math.abs(M[row * (2 * n) + col]) > Math.abs(M[piv * (2 * n) + col])) piv = row;
    }
    if (piv !== col) {
      for (j = 0; j < 2 * n; j++) {
        t = M[col * (2 * n) + j];
        M[col * (2 * n) + j] = M[piv * (2 * n) + j];
        M[piv * (2 * n) + j] = t;
      }
    }
    d = M[col * (2 * n) + col];
    if (Math.abs(d) < 1e-12) return null;
    for (j = 0; j < 2 * n; j++) M[col * (2 * n) + j] /= d;
    for (row = 0; row < n; row++) {
      if (row !== col) {
        f = M[row * (2 * n) + col];
        for (j = 0; j < 2 * n; j++) M[row * (2 * n) + j] -= f * M[col * (2 * n) + j];
      }
    }
  }
  R = new Float64Array(n * n);
  for (i = 0; i < n; i++) for (j = 0; j < n; j++) R[i * n + j] = M[i * (2 * n) + n + j];
  return R;
}

function solveLS(A, b, n) {
  var M = new Float64Array(n * (n + 1)), i, j, row, piv, d, f, t, x;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) M[i * (n + 1) + j] = A[i * n + j];
    M[i * (n + 1) + n] = b[i];
  }
  for (var col = 0; col < n; col++) {
    piv = col;
    for (row = col + 1; row < n; row++) {
      if (Math.abs(M[row * (n + 1) + col]) > Math.abs(M[piv * (n + 1) + col])) piv = row;
    }
    if (piv !== col) {
      for (j = 0; j <= n; j++) {
        t = M[col * (n + 1) + j];
        M[col * (n + 1) + j] = M[piv * (n + 1) + j];
        M[piv * (n + 1) + j] = t;
      }
    }
    d = M[col * (n + 1) + col];
    if (Math.abs(d) < 1e-12) return null;
    for (j = 0; j <= n; j++) M[col * (n + 1) + j] /= d;
    for (row = 0; row < n; row++) {
      if (row !== col) {
        f = M[row * (n + 1) + col];
        for (j = 0; j <= n; j++) M[row * (n + 1) + j] -= f * M[col * (n + 1) + j];
      }
    }
  }
  x = new Float64Array(n);
  for (i = 0; i < n; i++) x[i] = M[i * (n + 1) + n];
  return x;
}

function normP(z) {
  var u = 1 / (1 + 0.2316419 * z);
  var p = u * (0.31938153 + u * (-0.356563782 + u * (1.781477937 + u * (-1.821255978 + u * 1.330274429))));
  return (Math.exp(-0.5 * z * z) / 2.5066282746) * p;
}

function buildSubjs(data, nf) {
  var k = nf === 2 ? 6 : 9, sm = {}, id, fi, x, yv;
  data.forEach(function (row) {
    id = row.id;
    if (!sm[id]) sm[id] = { obs: [] };
    for (fi = 1; fi <= nf; fi++) {
      yv = row['y' + fi];
      if (yv == null) continue;
      x = new Float64Array(k);
      x[0] = 1;
      x[1] = row.y0;
      x[2] = row.trt;
      if (nf === 2) {
        if (fi === 2) {
          x[3] = 1;
          x[4] = row.trt;
          x[5] = row.y0;
        }
      } else {
        if (fi === 2) {
          x[3] = 1;
          x[5] = row.trt;
          x[7] = row.y0;
        } else if (fi === 3) {
          x[4] = 1;
          x[6] = row.trt;
          x[8] = row.y0;
        }
      }
      sm[id].obs.push({ j: fi - 1, y: yv, x: x });
    }
  });
  return Object.keys(sm)
    .filter(function (id) { return sm[id].obs.length > 0; })
    .map(function (id) { return sm[id]; });
}

function glsStep(subjs, k, sig, nf) {
  var XWX = new Float64Array(k * k), XWy = new Float64Array(k), r, c, a, b, ni, obs, w, wab, Si, SiInv;
  subjs.forEach(function (s) {
    obs = s.obs;
    ni = obs.length;
    if (!sig || ni === 1) {
      obs.forEach(function (o) {
        w = sig ? 1 / sig[o.j * nf + o.j] : 1;
        for (r = 0; r < k; r++) {
          XWy[r] += o.x[r] * o.y * w;
          for (c = 0; c < k; c++) XWX[r * k + c] += o.x[r] * o.x[c] * w;
        }
      });
    } else {
      Si = new Float64Array(ni * ni);
      for (a = 0; a < ni; a++) for (b = 0; b < ni; b++) Si[a * ni + b] = sig[obs[a].j * nf + obs[b].j];
      SiInv = matInv(Si, ni);
      if (!SiInv) {
        obs.forEach(function (o) {
          w = 1 / sig[o.j * nf + o.j];
          for (r = 0; r < k; r++) {
            XWy[r] += o.x[r] * o.y * w;
            for (c = 0; c < k; c++) XWX[r * k + c] += o.x[r] * o.x[c] * w;
          }
        });
        return;
      }
      for (a = 0; a < ni; a++) {
        for (b = 0; b < ni; b++) {
          wab = SiInv[a * ni + b];
          for (r = 0; r < k; r++) {
            XWy[r] += obs[a].x[r] * wab * obs[b].y;
            for (c = 0; c < k; c++) XWX[r * k + c] += obs[a].x[r] * wab * obs[b].x[c];
          }
        }
      }
    }
  });
  return { XWX: XWX, XWy: XWy };
}

function emUpd(subjs, beta, sig, nf) {
  var k = beta.length, SS = new Float64Array(nf * nf), nS = 0, j1, j2, a, b, i1, i2, Ej1, Ej2, cov12, c, wJ1, wJ2;
  subjs.forEach(function (s) {
    nS++;
    var obs = s.obs, ni = obs.length;
    if (!ni) return;
    var rO = new Float64Array(ni), oJ = [];
    obs.forEach(function (o, oi) {
      var yh = 0;
      for (var j = 0; j < k; j++) yh += o.x[j] * beta[j];
      rO[oi] = o.y - yh;
      oJ.push(o.j);
    });
    var Soo = new Float64Array(ni * ni);
    for (a = 0; a < ni; a++) for (b = 0; b < ni; b++) Soo[a * ni + b] = sig[oJ[a] * nf + oJ[b]];
    var SooInv = ni > 1 ? matInv(Soo, ni) : null;
    for (j1 = 0; j1 < nf; j1++) {
      i1 = oJ.indexOf(j1);
      wJ1 = null;
      if (i1 >= 0) {
        Ej1 = rO[i1];
      } else {
        wJ1 = new Float64Array(ni);
        if (ni === 1) wJ1[0] = sig[j1 * nf + oJ[0]] / (sig[oJ[0] * nf + oJ[0]] || 1e-9);
        else {
          for (a = 0; a < ni; a++) for (b = 0; b < ni; b++) wJ1[a] += sig[j1 * nf + oJ[b]] * SooInv[b * ni + a];
        }
        Ej1 = 0;
        for (a = 0; a < ni; a++) Ej1 += wJ1[a] * rO[a];
      }
      for (j2 = 0; j2 < nf; j2++) {
        i2 = oJ.indexOf(j2);
        wJ2 = null;
        if (i2 >= 0) {
          Ej2 = rO[i2];
        } else {
          wJ2 = new Float64Array(ni);
          if (ni === 1) wJ2[0] = sig[j2 * nf + oJ[0]] / (sig[oJ[0] * nf + oJ[0]] || 1e-9);
          else {
            for (a = 0; a < ni; a++) for (b = 0; b < ni; b++) wJ2[a] += sig[j2 * nf + oJ[b]] * SooInv[b * ni + a];
          }
          Ej2 = 0;
          for (a = 0; a < ni; a++) Ej2 += wJ2[a] * rO[a];
        }
        cov12 = 0;
        if (i1 < 0 && i2 < 0) {
          if (!wJ1) {
            wJ1 = new Float64Array(ni);
            if (ni === 1) wJ1[0] = sig[j1 * nf + oJ[0]] / (sig[oJ[0] * nf + oJ[0]] || 1e-9);
            else for (a = 0; a < ni; a++) for (b = 0; b < ni; b++) wJ1[a] += sig[j1 * nf + oJ[b]] * SooInv[b * ni + a];
          }
          c = sig[j1 * nf + j2];
          for (a = 0; a < ni; a++) c -= wJ1[a] * sig[oJ[a] * nf + j2];
          cov12 = c;
        }
        SS[j1 * nf + j2] += Ej1 * Ej2 + cov12;
      }
    }
  });
  var Neff = Math.max(nS - 3, 1), sg = new Float64Array(nf * nf), j, m;
  for (j = 0; j < nf * nf; j++) sg[j] = SS[j] / Neff;
  for (j1 = 0; j1 < nf; j1++) {
    for (j2 = j1 + 1; j2 < nf; j2++) {
      m = (sg[j1 * nf + j2] + sg[j2 * nf + j1]) / 2;
      sg[j1 * nf + j2] = m;
      sg[j2 * nf + j1] = m;
    }
    if (sg[j1 * nf + j1] < 1e-6) sg[j1 * nf + j1] = 1e-6;
  }
  return sg;
}

function sigInit(subjs, beta, nf) {
  var cc = subjs.filter(function (s) { return s.obs.length === nf; }), n = cc.length, j1, j2, res, yh, j;
  if (n < 4) return null;
  var sig = new Float64Array(nf * nf);
  cc.forEach(function (s) {
    res = new Float64Array(nf);
    s.obs.forEach(function (o) {
      yh = 0;
      for (j = 0; j < beta.length; j++) yh += o.x[j] * beta[j];
      res[o.j] = o.y - yh;
    });
    for (j1 = 0; j1 < nf; j1++) for (j2 = 0; j2 < nf; j2++) sig[j1 * nf + j2] += res[j1] * res[j2];
  });
  var Neff = Math.max(n - 3, 1);
  for (j = 0; j < nf * nf; j++) sig[j] /= Neff;
  for (j1 = 0; j1 < nf; j1++) if (sig[j1 * nf + j1] < 1e-6) sig[j1 * nf + j1] = 1e-6;
  return sig;
}

/**
 * Original HTML mmrm: returns p-value only (1 = failure / non-significant sentinel).
 * Preserved for bit-for-bit parity with the interactive app.
 */
function mmrm(data, nf) {
  var fit = fitMMRM(data, nf);
  return fit.ok ? fit.pValue : 1;
}

/**
 * EM-type GLS fit (same numerics as HTML mmrm), returning estimate, SE, p-value.
 * Contrast (unchanged from HTML):
 *   2 visits: treatment main + treatment×Visit2
 *   3 visits: treatment main + treatment×Visit3
 *
 * Also reports iterations and whether the relative-change tolerance was met.
 */
function fitMMRM(data, nf, options) {
  options = options || {};
  var maxIter = options.maxIter != null ? options.maxIter : 20;
  var tol = options.tol != null ? options.tol : 5e-4;
  var pValueMode = options.pValueMode || 'html'; // 'html' | 'normal' | 't'

  function fail(reason) {
    return {
      ok: false,
      estimate: null,
      standardError: null,
      pValue: null,
      failureReason: reason,
      iterations: 0,
      converged: false
    };
  }

  var k = nf === 2 ? 6 : 9, i, r, bt, Vb, cv, eff, ts, df, z, rr, cc, tmp, j;
  var contrast = new Float64Array(k);
  contrast[2] = 1;
  if (nf === 2) contrast[4] = 1;
  else contrast[6] = 1;

  var subjs = buildSubjs(data, nf), ns = subjs.length;
  if (ns < k + 2) return fail('insufficient_subjects');

  r = glsStep(subjs, k, null, nf);
  bt = solveLS(r.XWX, r.XWy, k);
  if (!bt) return fail('ols_solve_failed');

  var sig = sigInit(subjs, bt, nf);

  function studentTTail(zAbs, dfVal) {
    // Regularized incomplete beta via continued fraction (Abramowitz-style)
    // For ablation only; default HTML path does not use this.
    if (!(dfVal > 0) || !Number.isFinite(zAbs)) return 1;
    var x = dfVal / (dfVal + zAbs * zAbs);
    // Approximate 2*(1-F_|T|) via normal when df large
    if (dfVal > 120) return Math.min(1, 2 * normP(zAbs));
    // Simple Cornish-Fisher-ish fallback matching common stats libs poorly —
    // use normal * (1 + z^2/(2 df)) as upper bound alternative; for exact t
    // we use the incomplete-beta transformation.
    var a = dfVal / 2, b = 0.5, btVal = Math.exp(
      lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(Math.max(x, 1e-16)) + b * Math.log(Math.max(1 - x, 1e-16))
    );
    var cont = betacf(a, b, x);
    var ib = btVal * cont / a;
    var cdfCentral = 1 - 0.5 * ib; // P(|T|<=z) roughly for this transform
    // Standard: P(|T|>z) = incomplete beta
    var p = ib; // regularized incomplete beta I_x(df/2, 1/2) = P(|T| > z)
    return Math.min(1, Math.max(0, p));
  }

  function lgamma(xx) {
    // Lanczos approximation
    var cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
    var x = xx, y = xx, tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    var ser = 1.000000000190015;
    for (var j = 0; j < 6; j++) ser += cof[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }

  function betacf(a, b, x) {
    var MAXIT = 100, EPS = 3e-7, FPMIN = 1e-30;
    var qab = a + b, qap = a + 1, qam = a - 1;
    var c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    var h = d;
    for (var m = 1; m <= MAXIT; m++) {
      var m2 = 2 * m;
      var aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d;
      if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c;
      if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      var del = d * c;
      h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }

  function finalize(beta, XWX, nSubj, approxDf, iters, converged) {
    Vb = matInv(XWX, k);
    if (!Vb) return fail('variance_invert_failed');
    cv = 0;
    for (rr = 0; rr < k; rr++) {
      tmp = 0;
      for (cc = 0; cc < k; cc++) tmp += Vb[rr * k + cc] * contrast[cc];
      cv += contrast[rr] * tmp;
    }
    if (cv <= 0) return fail('nonpositive_contrast_variance');
    eff = 0;
    for (j = 0; j < k; j++) eff += contrast[j] * beta[j];
    var se = Math.sqrt(cv);
    if (!(se > 0) || !Number.isFinite(se) || !Number.isFinite(eff)) {
      return fail('nonfinite_estimate_or_se');
    }
    z = Math.abs(eff / se);
    df = Math.max(nSubj - 3, 1);
    var p;
    if (pValueMode === 'normal' || (!approxDf && pValueMode === 'html')) {
      p = Math.min(1, 2 * normP(z));
    } else if (pValueMode === 't') {
      p = studentTTail(z, df);
    } else if (approxDf) {
      // HTML default custom correction
      p = Math.min(1, 2 * normP(z) * (df > 120 ? 1 : 1 + z * z * (0.5 / df)));
    } else {
      p = Math.min(1, 2 * normP(z));
    }
    return {
      ok: true,
      estimate: eff,
      standardError: se,
      pValue: p,
      failureReason: null,
      iterations: iters,
      converged: converged
    };
  }

  // Fallback path when unstructured covariance cannot be initialized (same as HTML)
  if (!sig) {
    return finalize(bt, r.XWX, ns, false, 0, false);
  }

  var converged = false;
  var iters = 0;
  for (i = 0; i < maxIter; i++) {
    iters = i + 1;
    r = glsStep(subjs, k, sig, nf);
    bt = solveLS(r.XWX, r.XWy, k);
    if (!bt) break;
    var sn = emUpd(subjs, bt, sig, nf);
    var maxD = 0, dd;
    for (var q = 0; q < nf * nf; q++) {
      dd = Math.abs(sn[q] - sig[q]) / (Math.abs(sig[q]) + 0.01);
      if (dd > maxD) maxD = dd;
    }
    sig = sn;
    if (maxD < tol) {
      converged = true;
      break;
    }
  }

  r = glsStep(subjs, k, sig, nf);
  bt = solveLS(r.XWX, r.XWy, k);
  if (!bt) return fail('em_solve_failed');
  return finalize(bt, r.XWX, ns, true, iters, converged);
}

function simulateTrial(params, seed) {
  return simTrial(params, seed);
}

function trialSeed(baseSeed, trialIndex) {
  return baseSeed + trialIndex * 31337;
}

function runOneReplicate(params, seed, alpha) {
  var nf = params.fus.length;
  var data = simTrial(params, seed);
  var fit = fitMMRM(data, nf);
  return {
    seed: seed,
    ok: fit.ok,
    estimate: fit.estimate,
    standardError: fit.standardError,
    pValue: fit.pValue,
    rejected: fit.ok ? fit.pValue < alpha : false,
    failureReason: fit.failureReason,
    data: data
  };
}

module.exports = {
  mkR: mkR,
  rn: rn,
  chol: chol,
  mvr: mvr,
  scat: scat,
  simTrial: simTrial,
  simulateTrial: simulateTrial,
  simulateTrialStructured: simulateTrialStructured,
  toObservedWide: toObservedWide,
  empiricalFullDataTarget: empiricalFullDataTarget,
  empiricalObservedDataTarget: empiricalObservedDataTarget,
  generateTrial: generateTrial,
  matInv: matInv,
  solveLS: solveLS,
  normP: normP,
  buildSubjs: buildSubjs,
  glsStep: glsStep,
  emUpd: emUpd,
  sigInit: sigInit,
  mmrm: mmrm,
  fitMMRM: fitMMRM,
  trialSeed: trialSeed,
  runOneReplicate: runOneReplicate
};
