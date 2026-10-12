const M = require('./eo_runaway_model.js');
function base(over) {
  return Object.assign({
    V: 6, Mt: 400, progress: 0.5, T: 413.15, Pop: 3.62e5, MAWP: 4.45e5, acc: 0.10, F: 16/60, tFeed: 15, coolLost: false,
    kohWt: 0.003, phiP: 0.05, tauMix: 60, Tfeed: 303.15, fillMax: 0.6, Ppad: 101325, Tpad: 298.15, dTcool: 30,
    decomp: true, kScale: 1.0, Kd: 0.85, Kc: 1.0, tEnd: 900 }, over || {});
}
function show(label, inp) {
  const t0 = Date.now(); const r = M.run(inp); const dt = Date.now() - t0;
  const s = r.sizing;
  console.log(`${label}\n  pocket EO ${r.Mp.toFixed(1)} kg  w0=${r.wp0.toFixed(2)}  Tp0=${(r.Tp0-273.15).toFixed(0)}C | P0=${(r.P0/1e5).toFixed(2)} bar util=${r.utilization.toFixed(2)}` +
    ` | Wpeak=${r.peak.W.toFixed(2)} kg/s @${r.peak.t.toFixed(0)}s  vent starts ${r.tVent0===null?'never':r.tVent0.toFixed(0)+'s'}  ventMass=${r.ventMass.toFixed(0)} kg` +
    ` | Pmax=${(r.Pmax/1e5).toFixed(2)} bar Tb,max=${(r.Tmax-273.15).toFixed(0)}C Tp,max=${(r.TpMax-273.15).toFixed(0)}C dec=${r.decEO_kg.toFixed(3)} kg` +
    (s ? ` | omega=${s.omega.toFixed(1)} G=${s.G.toFixed(0)} A=${(s.A*1e4).toFixed(1)} cm2 (${s.A_in2.toFixed(2)} in2) orifice=${s.orifice||('>T x'+s.nT)}` : ' | no relief needed') +
    ` | steps=${r.trace.length} ${dt} ms`);
}
show('A. ESD works, 15 s, cooling normal', base());
show('A2. ESD 60 s', base({ tFeed: 60 }));
show('B. ESD fails, mixer restart 5 min', base({ tFeed: 300 }));
show('C. Power outage ESD 15 s, cooling lost', base({ coolLost: true }));
show('D. Power outage ESD fails, restart 5 min, cooling lost', base({ tFeed: 300, coolLost: true }));
show('B2. 10 min', base({ tFeed: 600 }));
