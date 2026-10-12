/*
 * EO / ethylene glycol -> PEG ethoxylation: loss-of-mixing runaway relief model (browser + Node).
 * Physics mirrors eo_runaway_model.py (bulk pool) and adds a two-zone EO pocket, tempered venting at
 * the relieving pressure, and API 520 / Leung omega orifice sizing.  SI internally (K, Pa, mol, kg, s).
 *
 * Sources: Di Serio 1996 (k11, EO Antoine), Zahnd 2026 (dHrxn, EO cp), Loetgering-Lin/BASF 2023
 * (dHrxn, decomposition), Crocco 1959 (decomposition kinetics), Chiu thesis (vessel/plant basis),
 * Gakhar & Carr 2003 (accumulation / interlock philosophy), API 521 s4.2, CCPS 3.2 (homogeneous venting).
 * Perry's Chemical Engineers' Handbook, 9th Ed. (Sec. 2, "Prediction and Correlation of Physical
 * Properties"): EO liquid density & heat of vaporization (DIPPR-105/106, compound #139, CAS 75-21-8),
 * EG liquid density (DIPPR-105, compound #137, CAS 107-21-1) and polyether liquid heat capacity
 * (Ruzicka-Domalski group contribution, Table 2-166). See CLAUDE.md for the PEG density N-trend caveat.
 */
(function (root) {
  'use strict';

  var GAM = (typeof module !== 'undefined' && module.exports)
    ? require('./eo_gamma_data.js') : root.EO_GAMMA;

  var R = 8.314462, RCAL = 1.9872;
  var MW_EO = 0.044053, MW_EG = 0.062068, MW_KOH = 0.056106, MW_N2 = 0.028014, MW_NC = 0.02203;
  var TC = 469.0, TB = 283.85;
  var DH = -92.0e3;              // J/mol EO (Di Serio/Zahnd/Dinh/BASF agree within ~3%)
  var DH_DEC = -132.0e3;         // J/mol EO -> CH4 + CO (BASF, Dinh)
  var K11 = { lnA: 20.52, E: 13000 };       // cm3/mol/s, cal/mol (Di Serio 1996, primary OH)
  var KD = { A: 1e11, E: 42000 };           // 1/s, cal/mol (Crocco 1959; gas phase, extrapolated)
  var ATM = 101325.0;

  // Perry's 9th Ed. DIPPR-105 saturated liquid density: rho [kmol/m3] = C1/C2^(1+(1-T/C3)^C4).
  // DIPPR-106 heat of vaporization: Hvap [J/mol] = 1e4*C1*(1-Tr)^(C2+C3*Tr+C4*Tr^2), Tr = T/Tc.
  // EO (#139): verified vs Di Serio 1996 Table 3 density at 100 C (747 vs 740 kg/m3, 1%) and vs the
  // "from memory" 25.5 kJ/mol latent-heat placeholder it replaces (DIPPR gives 25.57 kJ/mol at Tb, 0.3%).
  var DIPPR105_EO = { C1: 1.836, C2: 0.26024, C3: 469.15, C4: 0.2696, Tmin: 160.65, Tmax: 469.15 };
  var DIPPR106_EO = { Tc: 469.15, C1: 4.4514, C2: 1.1569, C3: -1.2336, C4: 0.50875 };
  // EG (#137): verified vs known EG density at 25 C (1.110 vs 1.113 g/mL, 0.3%).
  var DIPPR105_EG = { C1: 1.315, C2: 0.25125, C3: 720.0, C4: 0.21868, Tmin: 260.15, Tmax: 720.0 };
  function dippr105(T, p) { return p.C1 / Math.pow(p.C2, 1 + Math.pow(1 - T / p.C3, p.C4)); }              // kmol/m3
  function dippr106(T, p) { var Tr = T / p.Tc; return 1e4 * p.C1 * Math.pow(1 - Tr, p.C2 + p.C3 * Tr + p.C4 * Tr * Tr); } // J/mol

  // PEG (EG + N mol EO, N = mol EO per mol EG; N=0 is pure EG, N=1 is diethylene glycol, N=2 is
  // triethylene glycol, ...) liquid density: EG's own DIPPR-105 curve plus an N-dependent rise that
  // saturates toward a high-MW plateau. The plateau shape is calibrated to ONE real data point --
  // triethylene glycol (N=2), specific gravity 1.1252 at 20/20 C, Perry's 9th Ed. solubility/physical-
  // property table -- and a general literature value for the high-MW PEG/PEO plateau density
  // (~1.128 g/mL at 25 C). Simple molar-volume additivity (EG + N x liquid EO volume) was tried first
  // and UNDERSHOT the real TEG density by ~11%, so it was rejected in favor of this saturating fit.
  // Not independently verified beyond N=2 -- flagged in the tool's basis notes.
  var PEG_RHO_INF_OFFSET = 18.0;   // kg/m3, asymptotic density rise above EG (25 C basis)
  var PEG_RHO_NCHAR = 1.82;        // decay constant in N, fit through the TEG point (N=2)

  // ---------------------------------------------------------------- properties
  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }
  function psat(T) { T = clamp(T, 250, TC); return Math.exp(16.74 - 2568.0 / (T - 29.01)) * 133.322; }   // Di Serio 1996 Antoine
  function rhoEO(T) { T = clamp(T, DIPPR105_EO.Tmin, DIPPR105_EO.Tmax - 0.5); return dippr105(T, DIPPR105_EO) * (MW_EO * 1000); }  // kg/m3
  function rhoEG(T) { T = clamp(T, DIPPR105_EG.Tmin, DIPPR105_EG.Tmax - 0.5); return dippr105(T, DIPPR105_EG) * (MW_EG * 1000); }  // kg/m3
  function rhoPoly(T, nAvg) {
    var n = Math.max(nAvg || 0, 0);
    return rhoEG(T) + PEG_RHO_INF_OFFSET * (1 - Math.exp(-n / PEG_RHO_NCHAR));                            // kg/m3
  }
  function cpEO(T) { T = clamp(T, 233, 423); return 0.11017 * T * T - 60.153 * T + 10055.2936; }          // Zahnd 2026 eq 5, J/kg/K
  function Lv(T) { T = clamp(T, 0.35 * DIPPR106_EO.Tc, DIPPR106_EO.Tc - 0.1); return dippr106(T, DIPPR106_EO); }  // J/mol
  // Liquid heat capacity of the EG+N*EO chain (N = mol EO per mol EG; N=0 is pure EG, N=1 is diethylene
  // glycol, ...), Ruzicka-Domalski group contribution, Perry's Table 2-166: 2 diol-OH groups
  // [O}(H,C) diol] + (2N+2) (CH2 next to 1C,1O) [C}(2H,C,O)] + N plain ether-O [O}(2C)] (N=0: pure EG,
  // 2 CH2, 0 ether-O). Verified for N=0 (pure EG): 148.7 vs the real ~149.5 J/(mol K), 0.6%. The ether-O
  // group is tabulated valid only to 350 K; our 373-460 K operating range extrapolates that one term.
  var CPG = { CH2CO: [1.4596, 1.4657, -0.27140], OHdiol: [5.2302, -1.5124, 0.54075], etherO: [5.0312, -1.5718, 0.37860] };
  function cpPolyMolar(T, nAvg) {
    var n = Math.max(nAvg || 0, 0), nCH2 = 2 * n + 2;
    var A = nCH2 * CPG.CH2CO[0] + 2 * CPG.OHdiol[0] + n * CPG.etherO[0];
    var B = nCH2 * CPG.CH2CO[1] + 2 * CPG.OHdiol[1] + n * CPG.etherO[1];
    var D = nCH2 * CPG.CH2CO[2] + 2 * CPG.OHdiol[2] + n * CPG.etherO[2];
    var t = T / 100;
    return R * (A + B * t + D * t * t);                                  // J/mol/K
  }
  function cpPolyMass(T, nAvg) { return cpPolyMolar(T, nAvg) / (MW_EG + Math.max(nAvg || 0, 0) * MW_EO); } // J/kg/K
  function kArr(T) { return Math.exp(K11.lnA) * Math.exp(-K11.E / (RCAL * T)) * 1e-6; }   // m3/mol/s
  function kDec(T) { return KD.A * Math.exp(-KD.E / (RCAL * T)); }

  // Dip-leg nozzle sizing from feed rate. Chiu thesis Table 7.2 gives two real (ID, peak feed rate)
  // pairs -- 40.9 mm / 18 kg/min (QD2, QD3) and 52.5 mm / 28 kg/min (QD1, QD5) -- and both land within
  // 4% of the SAME liquid velocity (0.26 and 0.25 m/s) at their peak design feed rate. Taken as this
  // plant's real design velocity for liquid-EO dip-leg service (not a generic pipeline heuristic), and
  // used to size the pipe continuously with feed rate rather than a 2-point step function.
  var DIP_V_DESIGN = 0.255;   // m/s
  function dipIDForFeed(Fkgps, Tfeed) {
    var A = Fkgps / (rhoEO(Tfeed) * DIP_V_DESIGN);       // m2
    return 1000 * Math.sqrt(4 * A / Math.PI);             // mm
  }

  // gamma_EO(x, T, n) trilinear in (ln(1+n), x, T) from the UNIFAC grid
  var LN = GAM.n.map(function (n) { return Math.log(1 + n); });
  var NX = GAM.x.length, NT = GAM.T.length;
  function loc(arr, v) {
    var m = arr.length;
    if (v <= arr[0]) return [0, 0];
    if (v >= arr[m - 1]) return [m - 2, 1];
    var i = 0; while (arr[i + 1] < v) i++;
    return [i, (v - arr[i]) / (arr[i + 1] - arr[i])];
  }
  function gammaEO(x, T, n) {
    var a = loc(LN, Math.log(1 + Math.max(n, 0))), b = loc(GAM.x, x), c = loc(GAM.T, T);
    var g = GAM.g, s = 0;
    for (var i = 0; i < 2; i++) for (var j = 0; j < 2; j++) for (var k = 0; k < 2; k++) {
      var w = (i ? a[1] : 1 - a[1]) * (j ? b[1] : 1 - b[1]) * (k ? c[1] : 1 - c[1]);
      if (w > 0) s += w * g[(((a[0] + i) * NX) + (b[0] + j)) * NT + (c[0] + k)];
    }
    return s;
  }

  // ---------------------------------------------------------------- VLE flash (bulk liquid + headspace)
  function flash(T, nc, nr, nf, nN2, nNC, V) {
    var nfree = Math.max(nf, 0), nAvg = nc > 0 ? nr / nc : 0;
    var Vpoly = (nc * MW_EG + nr * MW_EO) / rhoPoly(T, nAvg);
    var Ps = psat(T), rE = rhoEO(T);
    function resid(nL) {
      var VL = Vpoly + nL * MW_EO / rE, Vg = Math.max(V - VL, 1e-6), x = nL / (nL + nc);
      return gammaEO(x, T, nAvg) * x * Ps * Vg / (R * T) + nL - nfree;
    }
    var nL;
    if (nfree <= 0) nL = 0;
    else if (resid(nfree) <= 0) nL = nfree;
    else {
      var lo = 0, hi = nfree;
      for (var it = 0; it < 55; it++) { var mid = 0.5 * (lo + hi); if (resid(mid) > 0) hi = mid; else lo = mid; }
      nL = 0.5 * (lo + hi);
    }
    var VL = Vpoly + nL * MW_EO / rE, Vg = Math.max(V - VL, 1e-6), x = nL / (nL + nc);
    var PEO = gammaEO(x, T, nAvg) * x * Ps;
    return { nL: nL, nV: nfree - nL, VL: VL, Vg: Vg, x: x, PEO: PEO, P: PEO + (nN2 + nNC) * R * T / Vg };
  }

  // ---------------------------------------------------------------- batch / normal operation
  function batch(inp) {
    var N = (inp.Mt - 62.068) / 44.053;                                    // mol EO per mol EG
    var mfinal = inp.fillMax * inp.V * rhoPoly(inp.T, N);                  // final (100%-conversion) product density
    var nc = mfinal / (inp.Mt / 1000.0);
    var nr = inp.progress * N * nc;
    var nK = inp.kohWt * mfinal / MW_KOH;
    var Vg0 = inp.V - nc * MW_EG / rhoPoly(298.15, 0);                     // pure-EG starter, before EO dosing
    var nN2 = inp.Ppad * Vg0 / (R * inp.Tpad);
    return { N: N, mfinal: mfinal, nc: nc, nr: nr, nK: nK, nN2: nN2,
             mliq: nc * MW_EG + nr * MW_EO, EOtotal: nc * N * MW_EO, EGmass: nc * MW_EG };
  }

  // normal operating state at (T, P_op): bisection on dissolved EO so that P = P_op
  function normalState(inp, b, T) {
    function P_of(nL) {
      var f = flash(T, b.nc, b.nr, nL, b.nN2, 0, inp.V);       // nfree = nL -> partition fixes nL if subcooled
      return f;
    }
    // search on total free EO nf
    var lo = 0, hi = b.nc * 8.0, f0 = flash(T, b.nc, b.nr, 0, b.nN2, 0, inp.V);
    var nitrogenExceeds = f0.P >= inp.Pop;
    var nf;
    if (nitrogenExceeds) nf = 0;
    else {
      for (var it = 0; it < 70; it++) {
        var mid = 0.5 * (lo + hi);
        if (flash(T, b.nc, b.nr, mid, b.nN2, 0, inp.V).P > inp.Pop) hi = mid; else lo = mid;
      }
      nf = 0.5 * (lo + hi);
    }
    var f = flash(T, b.nc, b.nr, nf, b.nN2, 0, inp.V);
    var rcap = kArr(T) * inp.kScale * f.nL * b.nK / f.VL;      // mol/s
    return { nf: nf, f: f, rcap: rcap, nitrogenExceeds: nitrogenExceeds };
  }

  // temperature at which the kinetic capacity equals the feed rate (bisection on T in [100, 200] C)
  function suggestT(inp) {
    var b = batch(inp), Fm = inp.F / MW_EO, lo = 373.15, hi = 473.15;
    function cap(T) { var bb = batch(Object.assign({}, inp, { T: T })); return normalState(inp, bb, T).rcap; }
    if (cap(hi) < Fm) return { T: hi, feedLimited: true, capHi: cap(hi) * MW_EO };
    if (cap(lo) > Fm) return { T: lo, feedLimited: false, capHi: cap(hi) * MW_EO };
    for (var i = 0; i < 40; i++) { var mid = 0.5 * (lo + hi); if (cap(mid) < Fm) lo = mid; else hi = mid; }
    return { T: 0.5 * (lo + hi), feedLimited: false, capHi: cap(473.15) * MW_EO };
  }

  // ---------------------------------------------------------------- dynamics
  // y = [Tb, nc, nr, nf, nN2, nNC, nK, ndec, Tp, nEp]
  function rhs(t, y, S) {
    var Tb = clamp(y[0], 250, 1500), nc = Math.max(y[1], 1e-9), nr = Math.max(y[2], 0), nf = Math.max(y[3], 0);
    var nN2 = Math.max(y[4], 0), nNC = Math.max(y[5], 0), nK = Math.max(y[6], 0);
    var Tp = clamp(y[8], 250, 1500), nEp = Math.max(y[9], 0);
    var phi = S.phiP, V = S.V;
    var f0 = flash(Tb, nc, nr, nf, nN2, nNC, V);
    var mliq = nc * MW_EG + nr * MW_EO;

    var rb = kArr(Tb) * S.kScale * f0.nL * nK * (1 - phi) / f0.VL;
    var rdec = S.decomp ? kDec(Tb) * f0.nV : 0;

    // pocket
    var nAvg = nc > 0 ? nr / nc : 0;
    var mp = phi * mliq, Vp = mp / rhoPoly(Tp, nAvg) + nEp * MW_EO / rhoEO(Tp);
    var rp = (nEp > 0 && Vp > 0) ? kArr(Tp) * S.kScale * nEp * phi * nK / Vp : 0;
    var cb = f0.nL / f0.VL, kex = 1.0 / S.tauMix;
    var Ex = kex * (nEp - cb * Vp);
    var xp = nEp / (nEp + phi * nc), pp = gammaEO(xp, Tp, nAvg) * xp * psat(Tp);
    var B = S.kBoil * nEp * Math.min(Math.max(pp / f0.P - 1, 0), 3);
    var Cp = mp * cpPolyMass(Tp, nAvg) + nEp * MW_EO * cpEO(Tp);
    var dTp = (-DH * rp - kex * Cp * (Tp - Tb) - Lv(Tp) * B) / Cp;
    var dnEp = -rp - Ex - B;

    // bulk sources (no vent)
    var src = [0, rb + rp, Ex + B - rb - rdec, 0, 2 * rdec, 0];     // nc, nr, nf, nN2, nNC, nK(->idx5)
    var Qcool = S.coolLost ? 0 : S.UA * (Tb - S.Tcw);
    var Q = -DH * rb - DH_DEC * rdec + kex * Cp * (Tp - Tb) - Qcool + Lv(Tb) * B;
    var Cb = (1 - phi) * mliq * cpPolyMass(Tb, nAvg) + f0.nL * MW_EO * cpEO(Tb);
    var Lb = Lv(Tb);

    // sensitivities of nV and P (forward differences) to T and inventory
    var dT = 0.02, fT = flash(Tb + dT, nc, nr, nf, nN2, nNC, V);
    var nV_T = (fT.nV - f0.nV) / dT, P_T = (fT.P - f0.P) / dT;
    var z = [nc, nr, nf, nN2, nNC], nV_j = [0, 0, 0, 0, 0], P_j = [0, 0, 0, 0, 0];
    for (var j = 0; j < 5; j++) {
      var eps = Math.max(1e-6 * Math.abs(z[j]), 1e-3), zz = z.slice(); zz[j] += eps;
      var fj = flash(Tb, zz[0], zz[1], zz[2], zz[3], zz[4], V);
      nV_j[j] = (fj.nV - f0.nV) / eps; P_j[j] = (fj.P - f0.P) / eps;
    }
    var D = Cb + Lb * nV_T, sumS = 0, sumZ = 0, sumPs = 0, sumPz = 0;
    for (j = 0; j < 5; j++) { sumS += nV_j[j] * src[j]; sumZ += nV_j[j] * z[j]; sumPs += P_j[j] * src[j]; sumPz += P_j[j] * z[j]; }
    var A = Q - Lb * sumS;
    var coefT = (Lb * (sumZ - f0.nV) - f0.P * V) / D;
    var g0 = P_T * (A / D) + sumPs, g1 = P_T * coefT - sumPz;

    var omega = 0;
    if (S.vent && g1 < 0) {
      var w = clamp((f0.P - 0.995 * S.Phold) / (0.005 * S.Phold), 0, 1);
      omega = w * Math.max(0, (g0 + Math.max(f0.P - S.Phold, 0) / S.tauP) / (-g1));
    }
    var dTb = A / D + omega * coefT;
    var Mpool = nc * MW_EG + nr * MW_EO + nf * MW_EO + nN2 * MW_N2 + nNC * MW_NC + nK * MW_KOH;
    var dy = [dTb, src[0] - omega * nc, src[1] - omega * nr, src[2] - omega * nf, -omega * nN2,
              src[4] - omega * nNC, -omega * nK, rdec, dTp, dnEp];
    var vapMass = f0.nV * MW_EO + nN2 * MW_N2 + nNC * MW_NC;
    return { dy: dy, P: f0.P, W: omega * Mpool, Mpool: Mpool, quality: vapMass / Mpool, Tb: Tb, Tp: Tp,
             wp: nEp * MW_EO / (nEp * MW_EO + mp + 1e-12), nEp: nEp, rb: rb, rp: rp, rdec: rdec,
             Qrx: -DH * (rb + rp), dTb: dTb, xb: f0.x, nV: f0.nV, nL: f0.nL, Vg: f0.Vg, omega: omega };
  }

  // Dormand-Prince 5(4) adaptive integrator; records diagnostics at every accepted step
  function integrate(S, y0, tEnd, stopFn) {
    var c = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
    var a = [[], [1 / 5], [3 / 40, 9 / 40], [44 / 45, -56 / 15, 32 / 9],
             [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
             [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
             [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84]];
    var b5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];
    var b4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];
    var n = y0.length, t = 0, y = y0.slice(), h = 0.05, out = [], steps = 0;
    var atol = [1e-4, 1e-2, 1e-2, 1e-2, 1e-4, 1e-3, 1e-3, 1e-4, 1e-4, 1e-3], rtol = 2e-5;
    var d0 = rhs(t, y, S); out.push(rec(t, y, d0));
    while (t < tEnd && steps < 60000) {
      h = Math.min(h, tEnd - t);
      var k = [d0.dy], ok = false, ynew, err;
      for (var s = 1; s < 7; s++) {
        var ys = y.slice();
        for (var q = 0; q < s; q++) { var aq = a[s][q]; if (aq) for (var i = 0; i < n; i++) ys[i] += h * aq * k[q][i]; }
        k.push(rhs(t + c[s] * h, ys, S).dy);
      }
      ynew = y.slice(); err = 0;
      for (i = 0; i < n; i++) {
        var s5 = 0, s4 = 0;
        for (q = 0; q < 7; q++) { s5 += b5[q] * k[q][i]; s4 += b4[q] * k[q][i]; }
        ynew[i] = y[i] + h * s5;
        var e = Math.abs(h * (s5 - s4)) / (atol[i] + rtol * Math.max(Math.abs(y[i]), Math.abs(ynew[i])));
        if (e > err) err = e;
      }
      if (!(err <= 1)) { h *= Math.max(0.2, 0.9 * Math.pow(err, -0.2)); if (h < 1e-6) break; steps++; continue; }
      t += h; y = ynew; steps++;
      d0 = rhs(t, y, S); out.push(rec(t, y, d0));
      h *= Math.min(5, Math.max(0.2, 0.9 * Math.pow(Math.max(err, 1e-10), -0.2)));
      if (stopFn && stopFn(t, y, d0)) break;
    }
    return out;
  }
  function rec(t, y, d) {
    return { t: t, W: d.W, P: d.P, Tb: d.Tb, Tp: d.Tp, wp: d.wp, nEp: d.nEp, rb: d.rb, rp: d.rp, Qrx: d.Qrx,
             quality: d.quality, Mpool: d.Mpool, xb: d.xb, dTb: d.dTb, ndec: y[7], nV: d.nV, nL: d.nL,
             Vg: d.Vg, omega: d.omega, y: y.slice() };
  }

  // ---------------------------------------------------------------- Leung omega / API 520 two-phase sizing
  // Open-system (stream) flash of the vessel pool composition: specific volume at the relieving pressure P0 and
  // after an isenthalpic 10% pressure drop (P9 = 0.9 P0); both from the SAME flash so omega is self-consistent.
  function streamVolumes(rec, P0, P9) {
    var y = rec.y, nc = y[1], nr = y[2], nf = y[3], nN2 = y[4], nNC = y[5];
    var T0 = rec.Tb, nAvg = nr / nc, inert = nN2 + nNC;
    function nLat(T, P) {                                  // liquid EO left at (T, P), open-system VLE
      var Ps = psat(T);
      function res(nL) {
        var x = nL / (nL + nc), nV = nf - nL, yEO = nV / (nV + inert);
        return gammaEO(x, T, nAvg) * x * Ps - yEO * P;
      }
      if (res(nf) <= 0) return nf;
      var lo = 0, hi = nf; for (var i = 0; i < 50; i++) { var m = 0.5 * (lo + hi); if (res(m) > 0) hi = m; else lo = m; }
      return 0.5 * (lo + hi);
    }
    function vol(T, nL, P) {
      var Vliq = (nc * MW_EG + nr * MW_EO) / rhoPoly(T, nAvg) + nL * MW_EO / rhoEO(T);
      return Vliq + (nf - nL + inert) * R * T / P;
    }
    var nL0 = nLat(T0, P0), nV0 = nf - nL0, M = nc * MW_EG + nr * MW_EO + nf * MW_EO + nN2 * MW_N2 + nNC * MW_NC + y[6] * MW_KOH;
    var Cl = (nc * MW_EG + nr * MW_EO) * cpPolyMass(T0, nAvg) + nL0 * MW_EO * cpEO(T0);
    var lo = T0 - 80, hi = T0;
    for (var it = 0; it < 45; it++) {
      var T9 = 0.5 * (lo + hi), nV9 = nf - nLat(T9, P9);
      var bal = Cl * (T0 - T9) - Lv(T9) * (nV9 - nV0);      // >0 => over-cooled: root is at higher T9
      if (bal > 0) lo = T9; else hi = T9;
    }
    var T9f = 0.5 * (lo + hi), nL9 = nLat(T9f, P9);
    return { v0: vol(T0, nL0, P0) / M, v9: vol(T9f, nL9, P9) / M, T9: T9f, vesselV: recSafeV(rec) / M };
  }
  function recSafeV(rec) { return rec.Vvessel || 0; }
  function etaC(om) {
    var lo = 1e-6, hi = 1;
    function f(e) { return e * e + (om * om - 2 * om) * Math.pow(1 - e, 2) + 2 * om * om * Math.log(e) + 2 * om * om * (1 - e); }
    for (var i = 0; i < 80; i++) { var m = 0.5 * (lo + hi); if (f(m) > 0) hi = m; else lo = m; }
    return 0.5 * (lo + hi);
  }
  var ORIFICES = [['D', 0.110], ['E', 0.196], ['F', 0.307], ['G', 0.503], ['H', 0.785], ['J', 1.287], ['K', 1.838],
                  ['L', 2.853], ['M', 3.60], ['N', 4.34], ['P', 6.38], ['Q', 11.05], ['R', 16.0], ['T', 26.0]];
  function size(recPk, S, opts) {
    opts = opts || {};
    var Kd = opts.Kd || 0.85, Kc = opts.Kc || 1.0, Pb = opts.Pb || ATM;
    var P0 = recPk.P, sv = streamVolumes(recPk, P0, 0.9 * P0), v0 = sv.v0, v9 = sv.v9;
    var omega = 9 * (v9 / v0 - 1);
    var G, eta = null, regime;
    if (omega < 0.05) {                                                // essentially gas/vapor
      var k = 1.3; G = Math.sqrt(k * (P0 / v0) * Math.pow(2 / (k + 1), (k + 1) / (k - 1))); regime = 'gas (omega<0.05)';
    } else {
      eta = etaC(omega); G = eta * Math.sqrt(P0 / (v0 * omega)); regime = 'omega method, critical flow';
      if (Pb / P0 > eta) regime = 'omega method, SUBCRITICAL (back pressure above critical)';
    }
    var A = recPk.W / (Kd * Kc * G);                                   // m2
    var A_in2 = A * 1550.0031, pick = null;
    for (var i = 0; i < ORIFICES.length; i++) if (ORIFICES[i][1] >= A_in2) { pick = ORIFICES[i]; break; }
    return { omega: omega, eta: eta, G: G, A: A, A_in2: A_in2, regime: regime, v0: v0, v9: v9, T9: sv.T9, vVessel: S.V / recPk.Mpool,
             orifice: pick ? pick[0] : null, orificeArea_in2: pick ? pick[1] : null,
             nT: pick ? 1 : Math.ceil(A_in2 / 26.0), Kd: Kd, Kc: Kc, dEq_mm: 1000 * Math.sqrt(4 * A / Math.PI) };
  }

  // ---------------------------------------------------------------- scenario driver
  // Pocket carrier: the liquid already sitting in the locally stagnant zone at the dip-leg exit when mixing
  // is lost. Sized from real pipe geometry (Chiu thesis Table 7.2: dip-leg ID 40.9 mm for <=8 m3 reactors,
  // 52.5 mm for larger) x a stagnant-zone length (in pipe diameters; no source -- flagged assumption), so the
  // carrier mass scales with the physical dead zone at the nozzle, not with an arbitrary fraction of the whole
  // batch. carrier mass = zone volume x liquid density (i.e. mass/density of the accumulated material itself).
  function pocketGeometry(inp, b) {
    var dID = inp.dipID_mm / 1000, Vzone = (Math.PI / 4) * dID * dID * (inp.zoneLenD * dID);
    var nAvg = b.nc > 0 ? b.nr / b.nc : 0;
    var mc = Vzone * rhoPoly(inp.T, nAvg);
    var phi = clamp(mc / b.mliq, 1e-4, 0.6);          // cap: the zone cannot exceed most of the batch liquid
    return { Vzone: Vzone, mc: mc, phi: phi };
  }

  function run(inp) {
    var b = batch(inp), ns = normalState(inp, b, inp.T);
    var pg = pocketGeometry(inp, b), phi = pg.phi;
    var Fm = inp.F / MW_EO, rb0 = kArr(inp.T) * inp.kScale * ns.f.nL * b.nK * (1 - phi) / ns.f.VL;
    var S = { V: inp.V, phiP: phi, tauMix: inp.tauMix, kScale: inp.kScale, decomp: inp.decomp, vent: true,
              coolLost: inp.coolLost, Phold: inp.MAWP * (1 + inp.acc), tauP: 2.0, kBoil: 0.5,
              Tcw: inp.T - inp.dTcool, UA: (-DH * rb0) / inp.dTcool };
    // Pocket EO = feed rate x response/restart time, PLUS -- when the ESD actually closes the valve
    // (modes A, C) -- the EO already standing in the line between the valve and the nozzle, which
    // keeps draining into the pocket after closure regardless of response time. That segment's length
    // is set by maintenance/bypass valve-clearance requirements, not by the dispersion physics that
    // size the submerged stagnant zone -- a distinct length, no source for a typical value (advanced
    // input, flagged). Mode B (ESD never closes) has no such one-time slug: the line stays topped up
    // and flows normally, so it contributes nothing beyond the steady feed x time already counted.
    var lineMass = inp.esdCloses ? (Math.PI / 4 * Math.pow(inp.dipID_mm / 1000, 2) * (inp.lineLen_m || 0) * rhoEO(inp.Tfeed)) : 0;
    var Mp = inp.F * inp.tFeed + lineMass, nEp0 = Mp / MW_EO, mliq = b.mliq, mc = pg.mc;
    var Ce = Mp * cpEO(inp.Tfeed), Cc = mc * cpPolyMass(inp.T, b.nc > 0 ? b.nr / b.nc : 0);
    var Tp0 = nEp0 > 0 ? (Cc * inp.T + Ce * inp.Tfeed) / (Cc + Ce) : inp.T;
    var y0 = [inp.T, b.nc, b.nr, ns.nf, b.nN2, 0, b.nK, 0, Tp0, nEp0];
    var P0 = ns.f.P;
    var stop = function (t, y, d) {
      return t > 6 * inp.tauMix && y[9] < 1e-4 * Math.max(nEp0, 1) && d.W < 1e-9 && Math.abs(d.dTb) < 0.005 && d.P < S.Phold * 0.9;
    };
    var trace = integrate(S, y0, inp.tEnd || 900, stop);
    var pk = trace[0], Wint = 0, Tmax = 0, Pmax = 0, tVent0 = null, ventMass = 0, decMax = 0, EOvent = 0;
    for (var i = 0; i < trace.length; i++) {
      var r = trace[i];
      if (r.W > pk.W) pk = r;
      if (r.Tb > Tmax) Tmax = r.Tb;
      if (r.P > Pmax) Pmax = r.P;
      if (tVent0 === null && r.W > 1e-6) tVent0 = r.t;
      if (i > 0) { var dt = r.t - trace[i - 1].t; ventMass += 0.5 * (r.W + trace[i - 1].W) * dt; }
      if (r.ndec > decMax) decMax = r.ndec;
    }
    var TpMax = Math.max.apply(null, trace.map(function (r) { return r.Tp; }));
    var out = { b: b, ns: ns, S: S, Mp: Mp, lineMass: lineMass, nEp0: nEp0, Tp0: Tp0, wp0: Mp / (Mp + mc + 1e-12), pocket: pg, trace: trace, peak: pk,
                Tmax: Tmax, TpMax: TpMax, Pmax: Pmax, tVent0: tVent0, ventMass: ventMass, decEO_kg: decMax * MW_EO, P0: P0,
                utilization: Fm / ns.rcap, rb0: rb0, needed: pk.W > 1e-6 };
    out.sizing = out.needed ? size(pk, S, { Kd: inp.Kd, Kc: inp.Kc }) : null;
    return out;
  }

  var api = { run: run, batch: batch, normalState: normalState, suggestT: suggestT, size: size, flash: flash,
              gammaEO: gammaEO, psat: psat, kArr: kArr, C: { MW_EO: MW_EO, MW_EG: MW_EG, DH: DH, ATM: ATM, R: R, TC: TC },
              ORIFICES: ORIFICES, rhoPoly: rhoPoly, rhoEG: rhoEG, rhoEO: rhoEO, Lv: Lv,
              cpPolyMolar: cpPolyMolar, cpPolyMass: cpPolyMass, dipIDForFeed: dipIDForFeed };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.EOModel = api;
})(typeof window !== 'undefined' ? window : this);
