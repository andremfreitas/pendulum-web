import { Chain, FORCING_FEEDBACK } from './sim.js';

// ------------------------------------------------------------ |
// parameters                                                   |

const params = {
    n: 12,
    nu: 1e-4,
    typ: 1,
    amf: 1.0,
    omf: 1.0,
    gr: 0.0,
    speed: 2.0,     // simulated time units per second
};

const INITIAL_SPREAD = 0.3;   // std of the random initial angles and velocities
const WARMUP_TIME    = 40.0;  // simulated before the first frame
const AVERAGE_TIME   = 60.0;  // e-folding time of the spectrum average
const FRAME_BUDGET   = 10.0;  // ms of integration per frame
const MAX_KICK       = 25.0;  // largest angular velocity a flick can set

const $ = id => document.getElementById(id);

const ui = {
    chain: $('chain'),
    spectrum: $('spectrum'),
    play: $('play'),
    sound: $('sound'),
    reset: $('reset'),
    status: $('status'),
    nlabel: $('nlabel'),
    omfLabel: $('omf-label'),
    time: $('r-time'),
    erg: $('r-erg'),
    pow: $('r-pow'),
    dis: $('r-dis'),
};

const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

// ------------------------------------------------------------ |
// number formatting                                            |

const SUP = { '-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹' };
const sup = s => [...String(s)].map(c => SUP[c] ?? c).join('');

function sci(v) {
    if (v === 0) return '0';
    const e = Math.floor(Math.log10(Math.abs(v)));
    const m = v / 10 ** e;
    if (e >= -2 && e <= 2) return v.toPrecision(2).replace(/\.?0+$/, '');
    return `${m.toFixed(1)} × 10${sup(e)}`;
}

function fixed(v, d = 2) {
    if (!Number.isFinite(v)) return '—';
    const a = Math.abs(v);
    if (a !== 0 && (a < 1e-3 || a >= 1e4)) return v.toExponential(1);
    return v.toFixed(d);
}

// ------------------------------------------------------------ |
// simulation state                                             |

let sim, spec, avg, avgW, trail;
let drag = null;   // flick gesture in progress
let pow = 0.0, dis = 0.0;
let running = !reducedMotion;

function build() {
    sim = new Chain({
        n: params.n, nu: params.nu, gr: params.gr,
        typ: params.typ, amf: params.amf, omf: params.omf,
    });
    sim.randomize(INITIAL_SPREAD, INITIAL_SPREAD);

    spec = new Float64Array(params.n);
    avg = new Float64Array(params.n);
    avgW = 0.0;
    trail = [];
    pow = dis = 0.0;
    drag = null;
    resetSound();

    // run a little so the first frame is already in motion
    const t0 = performance.now();
    while (sim.time < WARMUP_TIME && performance.now() - t0 < 150) {
        const before = sim.time;
        for (let k = 0; k < 50; k++) sim.step();
        accumulate(sim.time - before);
    }
    if (!sim.isFinite()) diverged();

    ui.nlabel.textContent = params.n;
    layoutChain();
}

function applyLive() {
    sim.nu = params.nu;
    sim.gr = params.gr;
    sim.typ = params.typ;
    sim.amf = params.amf;
    sim.omf = params.omf;
}

// exponential average of the energy spectrum, weighted by simulated time
function accumulate(dtSim) {
    sim.energySpectrum(spec);
    const a = 1.0 - Math.exp(-dtSim / AVERAGE_TIME);
    avgW += a * (1.0 - avgW);
    for (let i = 0; i < params.n; i++) avg[i] += a * (spec[i] - avg[i]);
    pow += a * (sim.pow - pow);
    dis += a * (sim.dis - dis);
}

function diverged() {
    setStatus('The integration diverged and was restarted. Try a smaller amplitude or a larger viscosity.');
    build();
}

let statusTimer = 0;
function setStatus(msg) {
    ui.status.textContent = msg;
    clearTimeout(statusTimer);
    if (msg) statusTimer = setTimeout(() => { ui.status.textContent = ''; }, 6000);
}

// ------------------------------------------------------------ |
// controls                                                     |

function logRange(input) {
    const [lo, hi] = input.dataset.log.split(' ').map(Number);
    const L = Math.log(lo), H = Math.log(hi);
    return {
        get: () => Math.exp(L + (H - L) * Number(input.value)),
        set: v => { input.value = String((Math.log(v) - L) / (H - L)); },
    };
}

const fmt = {
    n: v => String(v),
    nu: v => sci(v),
    amf: v => v.toFixed(2),
    omf: v => sci(v),
    gr: v => v.toFixed(2),
    speed: v => `${sci(v)} ×`,
};

const sliders = {};

for (const id of ['n', 'nu', 'amf', 'omf', 'gr', 'speed']) {
    const input = $(id);
    const out = $(`${id}-out`);
    const range = input.dataset.log
        ? logRange(input)
        : { get: () => Number(input.value), set: v => { input.value = String(v); } };

    // a log slider stores its position, so seed it from data-value
    if (input.dataset.log) range.set(Number(input.dataset.value));

    const read = () => {
        params[id] = id === 'n' ? Math.round(range.get()) : range.get();
        out.textContent = fmt[id](params[id]);
    };
    read();

    input.addEventListener('input', () => {
        const prev = params[id];
        read();
        if (id === 'n') {
            if (params.n !== prev) build();
        } else {
            applyLive();
        }
    });

    sliders[id] = { input, range, read };
}

function readForcing() {
    const checked = document.querySelector('input[name="typ"]:checked');
    params.typ = Number(checked.value);
    const feedback = params.typ === FORCING_FEEDBACK;
    ui.omfLabel.innerHTML = feedback ? 'Width <i>ε</i>' : 'Frequency <i>ω<sub>f</sub></i>';
    const off = params.typ === 0;
    sliders.amf.input.disabled = off;
    sliders.omf.input.disabled = off;
}

for (const r of document.querySelectorAll('input[name="typ"]')) {
    r.addEventListener('change', () => { readForcing(); applyLive(); });
}
readForcing();

$('controls').addEventListener('submit', e => e.preventDefault());

function setRunning(on) {
    running = on;
    ui.play.textContent = on ? 'Pause' : 'Play';
    ui.play.setAttribute('aria-pressed', String(on));
}

ui.play.addEventListener('click', () => setRunning(!running));
ui.sound.addEventListener('click', () => setSound(!soundOn));
ui.reset.addEventListener('click', () => { setStatus(''); build(); });

document.addEventListener('keydown', e => {
    if (e.target.closest('input, button')) return;
    if (e.code === 'Space') { e.preventDefault(); setRunning(!running); }
});

// ------------------------------------------------------------ |
// theme colours, read from the CSS tokens                      |

let colors = {};

function readColors() {
    const cs = getComputedStyle(document.documentElement);
    for (const k of ['bg', 'surface', 'ink', 'muted', 'line', 'accent']) colors[k] = cs.getPropertyValue(`--${k}`).trim();
    colors.mono = cs.getPropertyValue('--font-mono').trim() || 'monospace';
}

readColors();
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', readColors);
new MutationObserver(readColors).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme', 'class'] });

// ------------------------------------------------------------ |
// canvas sizing                                                |

function fitCanvas(canvas) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, Math.round(canvas.clientWidth * dpr));
    const h = Math.max(1, Math.round(canvas.clientHeight * dpr));
    if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w;
        canvas.height = h;
    }
    return dpr;
}

// ------------------------------------------------------------ |
// chain drawing                                                |
//
// Drawn link lengths taper geometrically from 1 to about 0.15,
// instead of the simulated factor 2 per link.

let disp = null;    // display lengths, sum = 1
let joints = null;  // screen positions of the joints, [x0, y0, x1, y1, ...]

function layoutChain() {
    const n = params.n;
    const r = Math.min(0.92, Math.max(0.6, 0.15 ** (1 / (n - 1))));
    disp = new Float64Array(n);
    let tot = 0.0;
    for (let i = 0; i < n; i++) { disp[i] = r ** i; tot += disp[i]; }
    for (let i = 0; i < n; i++) disp[i] /= tot;
    joints = new Float64Array(2 * (n + 1));
}

// tip position in units of the drawn chain length, pivot at the origin
function tip() {
    let x = 0.0, y = 0.0;
    for (let i = 0; i < params.n; i++) {
        x += disp[i] * Math.sin(sim.theta[i]);
        y += disp[i] * Math.cos(sim.theta[i]);
    }
    return [x, y];
}

function drawChain() {
    const cv = ui.chain;
    const dpr = fitCanvas(cv);
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    ctx.clearRect(0, 0, W, H);

    const n = params.n;
    const R = 0.46 * Math.min(W, H);
    const cx = W / 2, cy = H / 2;

    // joints
    joints[0] = cx; joints[1] = cy;
    for (let i = 0; i < n; i++) {
        joints[2 * i + 2] = joints[2 * i]     + R * disp[i] * Math.sin(sim.theta[i]);
        joints[2 * i + 3] = joints[2 * i + 1] + R * disp[i] * Math.cos(sim.theta[i]);
    }

    // reach circle
    ctx.strokeStyle = colors.line;
    ctx.lineWidth = dpr;
    ctx.setLineDash([2 * dpr, 4 * dpr]);
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);

    // tip trail
    const tnow = sim.time;
    const span = 1.5 * params.speed;
    ctx.strokeStyle = colors.accent;
    ctx.lineCap = 'round';
    ctx.lineWidth = 1.5 * dpr;
    for (let k = 1; k < trail.length; k++) {
        const [t0, x0, y0] = trail[k - 1];
        const [, x1, y1] = trail[k];
        const age = (tnow - t0) / span;
        if (age >= 1) continue;
        ctx.globalAlpha = 0.55 * (1 - age) ** 2;
        ctx.beginPath();
        ctx.moveTo(cx + R * x0, cy + R * y0);
        ctx.lineTo(cx + R * x1, cy + R * y1);
        ctx.stroke();
    }
    ctx.globalAlpha = 1;

    // links
    ctx.strokeStyle = colors.ink;
    for (let i = 0; i < n; i++) {
        const f = n > 1 ? i / (n - 1) : 0;
        ctx.lineWidth = (3.4 - 2.4 * f) * dpr;
        if (drag && drag.k === i) ctx.strokeStyle = colors.accent;
        ctx.beginPath();
        ctx.moveTo(joints[2 * i], joints[2 * i + 1]);
        ctx.lineTo(joints[2 * i + 2], joints[2 * i + 3]);
        ctx.stroke();
        ctx.strokeStyle = colors.ink;
    }

    // masses
    ctx.fillStyle = colors.ink;
    for (let i = 0; i < n; i++) {
        const f = n > 1 ? i / (n - 1) : 0;
        ctx.beginPath();
        ctx.arc(joints[2 * i + 2], joints[2 * i + 3], (3.6 - 2.2 * f) * dpr, 0, 2 * Math.PI);
        ctx.fill();
    }

    // pivot
    ctx.fillStyle = colors.bg;
    ctx.strokeStyle = colors.ink;
    ctx.lineWidth = 1.5 * dpr;
    ctx.beginPath();
    ctx.arc(cx, cy, 4.5 * dpr, 0, 2 * Math.PI);
    ctx.fill();
    ctx.stroke();

    // flick gesture
    if (drag && drag.pts.length) {
        const p = drag.pts[drag.pts.length - 1];
        ctx.strokeStyle = colors.accent;
        ctx.globalAlpha = 0.5;
        ctx.lineWidth = dpr;
        ctx.setLineDash([3 * dpr, 3 * dpr]);
        ctx.beginPath();
        ctx.moveTo(joints[2 * drag.k], joints[2 * drag.k + 1]);
        ctx.lineTo(p.x, p.y);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.globalAlpha = 1;
    }
}

// ------------------------------------------------------------ |
// flicking a link                                              |
//
// Press near a link, swing the pointer around the joint above it
// and release: the link takes the angular velocity of the swing.

function pointerPos(e) {
    const rect = ui.chain.getBoundingClientRect();
    const dpr = ui.chain.width / rect.width;
    return { x: (e.clientX - rect.left) * dpr, y: (e.clientY - rect.top) * dpr };
}

function nearestLink(p) {
    let best = -1, bestD = Infinity;
    for (let i = 0; i < params.n; i++) {
        const ax = joints[2 * i], ay = joints[2 * i + 1];
        const bx = joints[2 * i + 2], by = joints[2 * i + 3];
        const dx = bx - ax, dy = by - ay;
        const L2 = dx * dx + dy * dy || 1;
        const t = Math.max(0, Math.min(1, ((p.x - ax) * dx + (p.y - ay) * dy) / L2));
        const d = Math.hypot(p.x - ax - t * dx, p.y - ay - t * dy);
        if (d < bestD) { bestD = d; best = i; }
    }
    const dpr = ui.chain.width / ui.chain.getBoundingClientRect().width;
    return bestD < 28 * dpr ? best : -1;
}

function swingAngle(p, k) {
    return Math.atan2(p.x - joints[2 * k], p.y - joints[2 * k + 1]);
}

ui.chain.addEventListener('pointerdown', e => {
    const p = pointerPos(e);
    const k = nearestLink(p);
    if (k < 0) return;
    ui.chain.setPointerCapture(e.pointerId);
    ui.chain.classList.add('dragging');
    drag = { k, pts: [{ ...p, t: performance.now(), phi: swingAngle(p, k) }] };
});

ui.chain.addEventListener('pointermove', e => {
    if (!drag) return;
    const p = pointerPos(e);
    drag.pts.push({ ...p, t: performance.now(), phi: swingAngle(p, drag.k) });
    if (drag.pts.length > 64) drag.pts.shift();
});

function endDrag(e) {
    if (!drag) return;
    ui.chain.classList.remove('dragging');

    // angular velocity of the last ~120 ms of the gesture
    const now = performance.now();
    const pts = drag.pts.filter(q => now - q.t < 120);
    if (pts.length >= 2) {
        let dphi = 0.0;
        for (let j = 1; j < pts.length; j++) {
            let d = pts[j].phi - pts[j - 1].phi;
            if (d > Math.PI) d -= 2 * Math.PI;
            if (d < -Math.PI) d += 2 * Math.PI;
            dphi += d;
        }
        const dtReal = (pts[pts.length - 1].t - pts[0].t) / 1000;
        if (dtReal > 0) {
            // phi uses the same convention as theta; convert to simulated time
            const w = Math.max(-MAX_KICK, Math.min(MAX_KICK, dphi / dtReal / params.speed));
            if (Math.abs(w) > 1e-3) sim.omega[drag.k] = w;
        }
    }
    drag = null;
}

ui.chain.addEventListener('pointerup', endDrag);
ui.chain.addEventListener('pointercancel', () => { drag = null; ui.chain.classList.remove('dragging'); });

// ------------------------------------------------------------ |
// spectrum drawing                                             |

let yLo = -4, yHi = 1;

function drawSpectrum() {
    const cv = ui.spectrum;
    const dpr = fitCanvas(cv);
    const ctx = cv.getContext('2d');
    const W = cv.width, H = cv.height;
    ctx.clearRect(0, 0, W, H);

    const n = params.n;
    const norm = avgW > 0 ? 1 / avgW : 1;

    // y range from the average, eased so the axis does not jitter
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < n; i++) {
        const v = Math.log10(avg[i] * norm);
        if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    }
    if (Number.isFinite(lo)) {
        const tLo = Math.floor(lo - 0.3), tHi = Math.ceil(hi + 0.3);
        yLo += 0.08 * (tLo - yLo);
        yHi += 0.08 * (tHi - yHi);
        if (yHi - yLo < 1) yHi = yLo + 1;
    }

    ctx.font = `${11 * dpr}px ${colors.mono}`;
    const padL = 40 * dpr, padR = 6 * dpr, padT = 6 * dpr, padB = 20 * dpr;
    const X = k => padL + (n > 1 ? (k / (n - 1)) : 0.5) * (W - padL - padR);
    const Y = v => padT + (1 - (v - yLo) / (yHi - yLo)) * (H - padT - padB);

    // decades
    ctx.lineWidth = dpr;
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'right';
    const d0 = Math.ceil(yLo), d1 = Math.floor(yHi);
    const every = Math.max(1, Math.ceil((d1 - d0 + 1) / 5));
    for (let d = d0; d <= d1; d++) {
        if ((d - d0) % every) continue;
        const y = Y(d);
        ctx.strokeStyle = colors.line;
        ctx.beginPath();
        ctx.moveTo(padL, y);
        ctx.lineTo(W - padR, y);
        ctx.stroke();
        ctx.fillStyle = colors.muted;
        ctx.fillText(d === 0 ? '1' : `10${sup(d)}`, padL - 6 * dpr, y);
    }

    // link numbers
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const kEvery = n <= 12 ? 1 : n <= 18 ? 2 : 4;
    for (let i = 0; i < n; i++) {
        if ((i + 1) % kEvery && i !== 0) continue;
        ctx.fillStyle = colors.muted;
        ctx.fillText(String(i + 1), X(i), H - padB + 6 * dpr);
    }

    // instantaneous
    ctx.fillStyle = colors.muted;
    ctx.globalAlpha = 0.45;
    for (let i = 0; i < n; i++) {
        const v = Math.log10(spec[i]);
        if (!Number.isFinite(v) || v < yLo - 0.5 || v > yHi + 0.5) continue;
        ctx.beginPath();
        ctx.arc(X(i), Math.max(padT, Math.min(H - padB, Y(v))), 2 * dpr, 0, 2 * Math.PI);
        ctx.fill();
    }
    ctx.globalAlpha = 1;

    // average
    ctx.strokeStyle = colors.accent;
    ctx.fillStyle = colors.accent;
    ctx.lineWidth = 1.5 * dpr;
    ctx.beginPath();
    let started = false;
    for (let i = 0; i < n; i++) {
        const v = Math.log10(avg[i] * norm);
        if (!Number.isFinite(v)) continue;
        if (started) ctx.lineTo(X(i), Y(v)); else ctx.moveTo(X(i), Y(v));
        started = true;
    }
    ctx.stroke();
    for (let i = 0; i < n; i++) {
        const v = Math.log10(avg[i] * norm);
        if (!Number.isFinite(v)) continue;
        ctx.beginPath();
        ctx.arc(X(i), Y(v), 2.6 * dpr, 0, 2 * Math.PI);
        ctx.fill();
    }
}

// ------------------------------------------------------------ |
// sound: the pitch follows how fast the last link spins        |
//
// The map adapts: the running mean of log|omega_N| sits at the
// centre pitch and its spread sets the range, so the tone stays
// audible whatever the parameters. The tip's horizontal position
// sets the stereo pan.

const PITCH_CENTRE = 220.0;  // Hz
const PITCH_SPAN   = 1.5;    // octaves either side of the centre
const SOUND_MEMORY = 3.0;    // s, memory of the adaptive map

let audio = null;
let soundOn = false;
let logW = null, logVar = 1.0;

function resetSound() {
    logW = null;
    logVar = 1.0;
}

function startAudio() {
    const Ctx = window.AudioContext || window.webkitAudioContext;
    if (!Ctx) return false;

    const ctx = new Ctx();
    const osc = ctx.createOscillator();
    osc.type = 'triangle';
    osc.frequency.value = PITCH_CENTRE;

    const filter = ctx.createBiquadFilter();
    filter.type = 'lowpass';
    filter.frequency.value = 1800;

    const gain = ctx.createGain();
    gain.gain.value = 0.0;

    const pan = ctx.createStereoPanner ? ctx.createStereoPanner() : null;

    osc.connect(filter).connect(gain);
    if (pan) gain.connect(pan).connect(ctx.destination);
    else gain.connect(ctx.destination);
    osc.start();

    audio = { ctx, osc, gain, pan };
    return true;
}

function setSound(on) {
    if (on && !audio && !startAudio()) {
        setStatus('This browser cannot play sound.');
        on = false;
    }
    soundOn = on;
    if (audio) {
        if (on) audio.ctx.resume();
        else audio.gain.gain.setTargetAtTime(0.0, audio.ctx.currentTime, 0.05);
    }
    ui.sound.textContent = on ? 'Mute' : 'Listen';
    ui.sound.setAttribute('aria-pressed', String(on));
}

function updateSound(realDt) {
    if (!audio || !soundOn) return;
    const t = audio.ctx.currentTime;

    if (!running) {
        audio.gain.gain.setTargetAtTime(0.0, t, 0.08);
        return;
    }

    const lw = Math.log2(Math.abs(sim.omega[params.n - 1]) + 1e-9);
    if (logW === null) logW = lw;
    const a = 1.0 - Math.exp(-realDt / SOUND_MEMORY);
    logW += a * (lw - logW);
    logVar += a * ((lw - logW) ** 2 - logVar);

    // standardised deviation from the running mean
    const z = (lw - logW) / Math.sqrt(logVar + 1e-6);
    const oct = Math.max(-PITCH_SPAN, Math.min(PITCH_SPAN, 0.6 * z));
    const loud = Math.max(0.0, Math.min(1.0, (z + 1.0) / 3.0));

    audio.osc.frequency.setTargetAtTime(PITCH_CENTRE * 2 ** oct, t, 0.03);
    audio.gain.gain.setTargetAtTime(0.05 + 0.05 * loud, t, 0.05);
    if (audio.pan) {
        const [x] = tip();
        audio.pan.pan.setTargetAtTime(Math.max(-1.0, Math.min(1.0, x)), t, 0.05);
    }
}

// animation frames stop in a background tab, so stop the tone too
document.addEventListener('visibilitychange', () => {
    if (!audio) return;
    if (document.hidden) audio.ctx.suspend();
    else if (soundOn) audio.ctx.resume();
});

// ------------------------------------------------------------ |
// main loop                                                    |

let last = performance.now();
let rate = params.speed;
let frame = 0;

function advance(realDt) {
    const t0 = performance.now();
    const start = sim.time;
    const target = start + params.speed * realDt;
    const span = 1.5 * params.speed;

    while (sim.time < target && performance.now() - t0 < FRAME_BUDGET) {
        const before = sim.time;
        for (let k = 0; k < 20 && sim.time < target; k++) sim.step();
        if (!sim.isFinite()) { diverged(); return; }
        accumulate(sim.time - before);
        const [x, y] = tip();
        trail.push([sim.time, x, y]);
    }

    while (trail.length && sim.time - trail[0][0] > span) trail.shift();

    if (realDt > 0) rate += 0.1 * ((sim.time - start) / realDt - rate);
}

function readout() {
    ui.time.textContent = fixed(sim.time, 1);
    let e = 0.0;
    for (let i = 0; i < params.n; i++) e += spec[i];
    ui.erg.textContent = fixed(e + sim.potential(), 3);
    ui.pow.textContent = fixed(pow, 3);
    ui.dis.textContent = fixed(dis, 3);

    const lagging = running && rate < 0.8 * params.speed;
    if (lagging && !ui.status.textContent) {
        ui.status.textContent = `Running at ${sci(rate)} × instead of ${sci(params.speed)} ×: the small links need very short time steps.`;
    } else if (!lagging && ui.status.textContent.startsWith('Running at')) {
        ui.status.textContent = '';
    }
}

function loop(now) {
    const realDt = Math.min((now - last) / 1000, 0.1);
    last = now;

    if (running) advance(realDt);
    else rate = params.speed;
    updateSound(realDt);

    drawChain();
    drawSpectrum();
    if (frame++ % 6 === 0) readout();

    requestAnimationFrame(loop);
}

build();
setRunning(running);
readout();
requestAnimationFrame(loop);
