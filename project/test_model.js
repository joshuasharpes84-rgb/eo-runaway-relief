const M = require('./eo_runaway_model.js');
function base(over) {
  return Object.assign({
    V: 6, Mt: 400, progress: 0.5, T: 413.15, Pop: 3.62e5, MAWP: 4.45e5, acc: 0.10, F: 16/60, tFeed: 15, coolLost: false,
    kohWt: 0.003, phiP: 0.05, tauMix: 60, Tfeed: 303.15, fillMax: 0.6, Ppad: 101325, Tpad: 298.15, dTcool: 30,
    decomp: true, kScale: 1.0, Kd: 0.85, Kc: 1.0, tEnd: 900 }, over || {});
}
// 1. gamma check vs Python/thermo (EG-(EO)n, T=410 K x=0.10): n=3 -> ~0.55 ; n=10 -> ~0.23
console.log('gamma n=3 x=.10 T=413:', M.gammaEO(0.10, 413.15, 3).toFixed(3), ' n=10:', M.gammaEO(0.10, 413.15, 10).toFixed(3), ' n=30:', M.gammaEO(0.10, 413.15, 30).toFixed(3));
// 2. batch + normal state
const inp = base();
const b = M.batch(inp);
console.log('N=%s nc=%s mol  mfinal=%s kg  nK=%s mol  nN2=%s mol  EOtotal=%s kg', b.N.toFixed(2), b.nc.toFixed(0), b.mfinal.toFixed(0), b.nK.toFixed(0), b.nN2.toFixed(0), b.EOtotal.toFixed(0));
const ns = M.normalState(inp, b, inp.T);
console.log('normal: nf=%s mol  x=%s  P=%s bar  PEO=%s bar  rcap=%s mol/s (%s kg/min)  N2exceeds=%s', ns.nf.toFixed(1), ns.f.x.toFixed(4), (ns.f.P/1e5).toFixed(3), (ns.f.PEO/1e5).toFixed(3), ns.rcap.toFixed(3), (ns.rcap*44.053e-3*60).toFixed(1), ns.nitrogenExceeds);
const sT = M.suggestT(inp);
console.log('suggested T for F=%s kg/min: %s C  feedLimited=%s  capacity@200C=%s kg/min', (inp.F*60).toFixed(1), (sT.T-273.15).toFixed(1), sT.feedLimited, (sT.capHi*60).toFixed(1));
// 3. progress sweep of suggested T at fixed feed
for (const p of [0.1, 0.3, 0.5, 0.8, 1.0]) {
  const i2 = base({ progress: p }); const s2 = M.suggestT(i2);
  const b2 = M.batch(i2); const n2 = M.normalState(i2, b2, s2.T);
  console.log('progress %s%%: T_sug=%s C feedLimited=%s  P_N2-only? nitrogenExceeds=%s x_ss=%s', (p*100).toFixed(0), (s2.T-273.15).toFixed(1), s2.feedLimited, n2.nitrogenExceeds, n2.f.x.toFixed(3));
}
