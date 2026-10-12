const M = require('./eo_runaway_model.js');
const base = (o) => Object.assign({ V:6, Mt:400, progress:0.5, T:393.15, Pop:3.62e5, MAWP:4.45e5, acc:0.10, F:16/60, tFeed:15, coolLost:false,
  kohWt:0.003, phiP:0.05, tauMix:60, Tfeed:303.15, fillMax:0.6, Ppad:101325, Tpad:298.15, dTcool:30, decomp:true, kScale:1, Kd:0.85, Kc:1, tEnd:900 }, o||{});
function row(l, o){ const r=M.run(base(o)); console.log(l.padEnd(46), 'pocket',r.Mp.toFixed(0).padStart(4),'kg | W',r.peak.W.toFixed(1).padStart(6),'kg/s | Pmax',(r.Pmax/1e5).toFixed(2),'| Tb',(r.Tmax-273.15).toFixed(0),'C | orifice', r.sizing? (r.sizing.orifice||'>T')+' ('+(r.sizing.A*1e4).toFixed(0)+' cm2)':'none'); return r; }
console.log('PEG 400, 6 m3, 50% progress, 120 C, 3.62 bar a, MAWP 4.45, feed 16 kg/min');
row('A  ESD 15 s, cooling on', {tFeed:15}); row('A  ESD 60 s', {tFeed:60});
row('C  outage, ESD 15 s, cooling lost', {tFeed:15, coolLost:true});
for (const m of [2,5,8,10,15]) row('B  restart '+m+' min, cooling on', {tFeed:m*60});
for (const m of [2,5,8,10,15]) row('D  restart '+m+' min, cooling lost', {tFeed:m*60, coolLost:true});
// onset of relief in B: scan delay
let lo=60, hi=1800; for(let i=0;i<14;i++){const mid=.5*(lo+hi); const r=M.run(base({tFeed:mid})); if(r.needed) hi=mid; else lo=mid;} console.log('B: relief begins at restart delay ~', (0.5*(lo+hi)/60).toFixed(1),'min  =', (16*0.5*(lo+hi)/60).toFixed(0),'kg pocket');
