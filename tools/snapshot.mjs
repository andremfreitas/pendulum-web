// ------------------------------------------------------------ |
// Pick a frame of the real simulation for the link preview     |
//                                                              |
// Runs sim.js with the site's default parameters and a fixed   |
// seed, then keeps an extended frame that bends one way (a    |
// hook, not a crumple). Writes the chain state, the tip trail  |
// and the energy spectrum to snapshot.json for og.py.          |
//                                                              |
//   node tools/snapshot.mjs                                    |

import { writeFileSync } from 'node:fs';
import { Chain } from '../sim.js';

const params = { n: 12, nu: 1e-4, typ: 1, amf: 1.0, omf: 1.0, gr: 0.0 };
const SPEED = 2.0;             // site default; sets the trail length
const TRAIL = 1.5 * SPEED;     // simulated time shown in the trail
const T_SKIP = 100.0;          // transient before frames are considered
const T_END = 400.0;
const AVERAGE_TIME = 60.0;     // as app.js

// display taper, as layoutChain in app.js
const n = params.n;
const r = Math.min(0.92, Math.max(0.6, 0.15 ** (1 / (n - 1))));
const disp = Array.from({ length: n }, (_, i) => r ** i);
const tot = disp.reduce((a, b) => a + b, 0);
for (let i = 0; i < n; i++) disp[i] /= tot;

const tip = th => {
    let x = 0, y = 0;
    for (let i = 0; i < n; i++) { x += disp[i] * Math.sin(th[i]); y += disp[i] * Math.cos(th[i]); }
    return [x, y];
};

// net and total turning of the drawn chain, and its sharpest
// joint, in radians
const turning = th => {
    let net = 0, tot = 0, sharp = 0;
    for (let i = 1; i < n; i++) {
        let d = th[i] - th[i - 1];
        d -= 2 * Math.PI * Math.round(d / (2 * Math.PI));
        net += d;
        tot += Math.abs(d);
        sharp = Math.max(sharp, Math.abs(d));
    }
    return [Math.abs(net), tot, sharp];
};

const MAX_JOINT = 50 * Math.PI / 180;   // sharper joints look folded

let seed = 20261103;
const rng = () => (seed = (seed * 16807) % 2147483647) / 2147483647;

const sim = new Chain(params);
sim.randomize(0.3, 0.3, rng);

const spec = new Float64Array(n), avg = new Float64Array(n);
let avgW = 0;
const trail = [];
let best = null;

while (sim.time < T_END) {
    const before = sim.time;
    for (let k = 0; k < 20; k++) {
        sim.step();
        trail.push([sim.time, ...tip(sim.theta)]);
    }
    if (!sim.isFinite()) throw new Error('integration diverged');

    sim.energySpectrum(spec);
    const a = 1 - Math.exp(-(sim.time - before) / AVERAGE_TIME);
    avgW += a * (1 - avgW);
    for (let i = 0; i < n; i++) avg[i] += a * (spec[i] - avg[i]);

    while (sim.time - trail[0][0] > TRAIL) trail.shift();

    if (sim.time < T_SKIP) continue;
    const [x, y] = tip(sim.theta);
    const reach = Math.hypot(x, y);
    if (reach < 0.55 || reach > 0.9) continue;
    const [net, tot, sharp] = turning(sim.theta);
    if (tot > 1.2 * net || sharp > MAX_JOINT) continue;
    const score = net;
    if (!best || score > best.score) {
        best = {
            score, time: sim.time,
            theta: [...sim.theta],
            trail: trail.map(([t, tx, ty]) => [(sim.time - t) / TRAIL, tx, ty]),
            spec: [...spec],
            avg: [...avg].map(v => v / avgW),
        };
    }
}

const out = { params, disp, ...best };
writeFileSync(new URL('./snapshot.json', import.meta.url), JSON.stringify(out));
console.log(`frame at t = ${best.time.toFixed(2)}, score ${best.score.toFixed(2)}`);
