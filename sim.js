// ------------------------------------------------------------ |
// N-link pendulum integrator: a port of sim_code/src/f90       |
//                                                              |
// Link i (0-based) has length and end mass 2^-i. Angles are    |
// measured from the downward vertical. Time stepping is the    |
// IMEX-SSP3(4,3,3) scheme of main.f90: dissipation implicit,   |
// everything else explicit.                                    |

// IMEX-SSP3(4,3,3) coefficients (declarations.f90, MODULE solver)
const alpha = 0.24169426078821;
const betaS = alpha / 4.0;
const etaS  = 0.12915286960590;

const ce = [0.0, 0.0, 1.0, 0.5];

const ae = [
    [],
    [0.0],
    [0.0, 1.0],
    [0.0, 0.25, 0.25],
];

const ai = [
    [alpha],
    [-alpha, alpha],
    [0.0, 1.0 - alpha, alpha],
    [betaS, etaS, 0.5 - betaS - etaS - alpha, alpha],
];

const bw = [0.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];

export const FORCING_NONE     = 0;
export const FORCING_PERIODIC = 1;
export const FORCING_FEEDBACK = 2;

export class Chain {

    // n     number of links (>= 2)
    // nu    viscosity
    // gr    gravitational acceleration
    // typ   forcing type: 0 none, 1 periodic, 2 feedback
    // amf   forcing amplitude
    // omf   forcing frequency (periodic) or regularisation (feedback)
    // beta  CFL parameter
    constructor({ n = 12, nu = 1e-6, gr = 0.0, typ = 1, amf = 0.1, omf = 0.01, beta = 0.1 } = {}) {

        if (n < 2) throw new RangeError('the chain needs at least 2 links');

        this.n = n;
        this.nu = nu;
        this.gr = gr;
        this.typ = typ;
        this.amf = amf;
        this.omf = omf;
        this.beta = beta;

        const F = () => new Float64Array(n);

        // --> metric <------------------------------------------ |

        this.ll = F();
        this.mo = F();
        this.mm = F();
        this.gg = F();
        this.ee = F();

        for (let i = 0; i < n; i++) {
            this.ll[i] = 1.0 / 2.0 ** i;
            this.mo[i] = 1.0 / 2.0 ** i;
        }

        this.gg[0] = 1.0 / this.ll[1];
        for (let i = 1; i < n - 1; i++) this.gg[i] = 1.0 / this.ll[i] + 1.0 / this.ll[i + 1];
        this.gg[n - 1] = 1.0 / this.ll[n - 1];

        for (let i = 0; i < n - 1; i++) this.ee[i] = -1.0 / this.ll[i + 1];
        this.ee[n - 1] = 0.0;

        this.mm[n - 1] = this.mo[n - 1];
        for (let i = n - 2; i >= 0; i--) this.mm[i] = this.mm[i + 1] + this.mo[i];

        // --> state <------------------------------------------- |

        this.theta = F();
        this.omega = F();
        this.time = 0.0;
        this.dt = 0.0;

        // power injected by the forcing and dissipated by viscosity
        // during the last step (same definitions as Checks/energy.txt)
        this.pow = 0.0;
        this.dis = 0.0;

        // --> work arrays <------------------------------------- |

        this.ths = [F(), F(), F(), F()];
        this.oms = [F(), F(), F(), F()];
        this.fom = [F(), F(), F(), F()];
        this.gim = [F(), F(), F(), F()];
        this.tau = [F(), F(), F(), F()];

        this.a1 = F();
        this.a2 = F();
        this.a3 = F();

        this.c0 = F(); this.s0 = F();
        this.c1 = F(); this.s1 = F();
        this.al = F(); this.be = F(); this.ga = F();
        this.zz = F(); this.bb = F();
    }

    // Gaussian initial condition with standard deviations a0 (angles)
    // and w0 (angular velocities), as initial.f90 (different RNG)
    randomize(a0 = 0.01, w0 = 0.01, rng = Math.random) {
        const gauss = () => {
            let u = 0.0;
            while (u === 0.0) u = rng();
            return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * rng());
        };
        for (let i = 0; i < this.n; i++) this.omega[i] = w0 * gauss();
        for (let i = 0; i < this.n; i++) this.theta[i] = a0 * gauss();
        this.time = 0.0;
    }

    // ------------------------------------------------------------ |
    // c0 = l cos(th), s0 = l sin(th), c1 = mm c0, s1 = mm s0       |

    trig(th) {
        const { n, ll, mm, c0, s0, c1, s1 } = this;
        for (let i = 0; i < n; i++) {
            c0[i] = ll[i] * Math.cos(th[i]);
            s0[i] = ll[i] * Math.sin(th[i]);
            c1[i] = c0[i] * mm[i];
            s1[i] = s0[i] * mm[i];
        }
    }

    // ------------------------------------------------------------ |
    // centrifugal + gravitational terms (expects trig(th) done)    |
    //   ff(i) = sum_j mm(max(i,j)) l(i) l(j) sin(th(j)-th(i)) om(j)^2
    //         - gr mm(i) l(i) sin(th(i))
    // evaluated in O(n) with prefix / suffix sums

    centrifugal(om, ff) {
        const { n, mm, c0, s0, gr } = this;

        // suffix sums over j > i of mm(j) l(j) (cos, sin)(th(j)) om(j)^2
        let sc = 0.0, ss = 0.0;
        for (let i = n - 1; i >= 0; i--) {
            ff[i] = 0.0;
            this.al[i] = sc;
            this.be[i] = ss;
            const w2 = om[i] * om[i];
            sc += mm[i] * c0[i] * w2;
            ss += mm[i] * s0[i] * w2;
        }

        // prefix sums over j <= i of l(j) (cos, sin)(th(j)) om(j)^2
        let pc = 0.0, ps = 0.0;
        for (let i = 0; i < n; i++) {
            const w2 = om[i] * om[i];
            pc += c0[i] * w2;
            ps += s0[i] * w2;

            const C = mm[i] * pc + this.al[i];
            const S = mm[i] * ps + this.be[i];

            // l(i) (cos(th_i) S - sin(th_i) C) = l(i) sum ... sin(th_j - th_i)
            ff[i] = c0[i] * S - s0[i] * C - gr * mm[i] * s0[i];
        }
    }

    // ------------------------------------------------------------ |
    // forcing torque (forcing.f90)                                 |

    forcing(om, time, tau) {
        tau.fill(0.0);
        if (this.typ === FORCING_PERIODIC) {
            const f = this.amf * Math.cos(this.omf * time);
            tau[0] =  f;
            tau[1] = -f;
        } else if (this.typ === FORCING_FEEDBACK) {
            let f = om[0] - om[1];
            f = this.amf * f / (f * f + this.omf * this.omf);
            tau[0] =  f;
            tau[1] = -f;
        }
    }

    // ------------------------------------------------------------ |
    // solve M(th) xx = bb with the metric matrix (expects trig)    |

    solveLinear(bb, xx) {
        const { n, c0, s0, c1, s1, al, be, zz } = this;

        let uu = 0.0, vv = 0.0, ww = 0.0;
        let pp = 0.0, qq = 0.0;

        for (let i = 0; i < n; i++) {
            const dd = c1[i] * c0[i] + s1[i] * s0[i]
                     - (c1[i] * c1[i] * uu + 2 * c1[i] * s1[i] * vv + s1[i] * s1[i] * ww);

            al[i] = (c0[i] - c1[i] * uu - s1[i] * vv) / dd;
            be[i] = (s0[i] - s1[i] * ww - c1[i] * vv) / dd;

            const yy = bb[i] - c1[i] * pp - s1[i] * qq;
            zz[i] = yy / dd;

            uu += dd * al[i] * al[i];
            vv += dd * al[i] * be[i];
            ww += dd * be[i] * be[i];

            pp += al[i] * yy;
            qq += be[i] * yy;
        }

        let rr = 0.0, ss = 0.0;
        for (let i = n - 1; i >= 0; i--) {
            xx[i] = zz[i] - al[i] * rr - be[i] * ss;
            rr += c1[i] * xx[i];
            ss += s1[i] * xx[i];
        }
    }

    // ------------------------------------------------------------ |
    // solve (M(th) + rho D) xx = M(th) rh, D the tridiagonal       |
    // dissipation matrix (expects trig)                            |

    solveLinearE(rh, rho, xx) {
        const { n, c0, s0, c1, s1, gg, ee, al, be, ga, zz, bb } = this;

        // bb = M rh, with M(i,j) = c1(max) c0(min) + s1(max) s0(min)
        let A = 0.0, B = 0.0;
        for (let i = 0; i < n; i++) {
            A += c0[i] * rh[i];
            B += s0[i] * rh[i];
            bb[i] = c1[i] * A + s1[i] * B;
        }
        let C = 0.0, D = 0.0;
        for (let i = n - 1; i >= 0; i--) {
            bb[i] += c0[i] * C + s0[i] * D;
            C += c1[i] * rh[i];
            D += s1[i] * rh[i];
        }

        let uu = 0.0, vv = 0.0, ww = 0.0;
        let pp = 0.0, qq = 0.0;
        let dp = 0.0, yo = 0.0;
        let alp = 0.0, bep = 0.0, gap = 0.0;

        for (let i = 0; i < n; i++) {
            const dd = c1[i] * c0[i] + s1[i] * s0[i]
                     - (c1[i] * c1[i] * uu + 2 * c1[i] * s1[i] * vv + s1[i] * s1[i] * ww)
                     + rho * gg[i]
                     - dp * gap * (gap + 2 * (c1[i] * alp + s1[i] * bep));

            al[i] = (c0[i] - c1[i] * uu - s1[i] * vv - dp * gap * alp) / dd;
            be[i] = (s0[i] - s1[i] * ww - c1[i] * vv - dp * gap * bep) / dd;

            const yy = bb[i] - c1[i] * pp - s1[i] * qq - gap * yo;
            zz[i] = yy / dd;

            uu += dd * al[i] * al[i];
            vv += dd * al[i] * be[i];
            ww += dd * be[i] * be[i];

            pp += al[i] * yy;
            qq += be[i] * yy;

            ga[i] = rho * ee[i] / dd;

            dp = dd;
            yo = yy;
            alp = al[i]; bep = be[i]; gap = ga[i];
        }

        xx[n - 1] = zz[n - 1];
        let rr = c1[n - 1] * xx[n - 1];
        let ss = s1[n - 1] * xx[n - 1];
        for (let i = n - 2; i >= 0; i--) {
            xx[i] = zz[i] - al[i] * rr - be[i] * ss - ga[i] * xx[i + 1];
            rr += c1[i] * xx[i];
            ss += s1[i] * xx[i];
        }
    }

    // ------------------------------------------------------------ |
    // time step from the CFL conditions (main.f90 and cfl)         |

    timeStep() {
        const { n, omega, theta, a1, a2, a3, beta } = this;

        let wmax = 0.0;
        for (let i = 0; i < n; i++) wmax = Math.max(wmax, Math.abs(omega[i]));
        let dt = beta / (wmax + 1.0);

        this.trig(theta);
        this.centrifugal(omega, a1);
        this.forcing(omega, this.time, a3);
        for (let i = 0; i < n; i++) a1[i] += a3[i];
        this.solveLinear(a1, a2);

        let amax = 0.0;
        for (let i = 0; i < n; i++) amax = Math.max(amax, Math.abs(a2[i]));

        return Math.min(dt, beta / (Math.sqrt(amax) + 1.0));
    }

    // ------------------------------------------------------------ |
    // one IMEX-SSP3(4,3,3) step                                    |

    step() {
        const { n, theta, omega, ths, oms, fom, gim, tau, a1, a2, nu } = this;

        const dt = this.timeStep();

        for (let s = 0; s < 4; s++) {

            const th = ths[s], om = oms[s];

            // explicit predictor for angles and velocities
            for (let i = 0; i < n; i++) {
                let t = theta[i];
                let w = omega[i];
                for (let r = 0; r < s; r++) {
                    t += dt * ae[s][r] * oms[r][i];
                    w += dt * (ae[s][r] * fom[r][i] + ai[s][r] * gim[r][i]);
                }
                th[i] = t;
                a2[i] = w;
            }

            // implicit dissipation
            const rho = dt * ai[s][s] * nu;

            this.trig(th);
            this.solveLinearE(a2, rho, om);

            // explicit terms at the stage values
            this.centrifugal(om, a1);
            this.forcing(om, this.time + dt * ce[s], tau[s]);
            for (let i = 0; i < n; i++) a1[i] += tau[s][i];
            this.solveLinear(a1, fom[s]);

            const inv = 1.0 / (dt * ai[s][s]);
            for (let i = 0; i < n; i++) gim[s][i] = (om[i] - a2[i]) * inv;
        }

        // final update
        for (let i = 0; i < n; i++) {
            let t = 0.0, w = 0.0;
            for (let s = 0; s < 4; s++) {
                t += bw[s] * oms[s][i];
                w += bw[s] * (fom[s][i] + gim[s][i]);
            }
            theta[i] += dt * t;
            omega[i] += dt * w;
        }

        // injection and dissipation rates (Checks/energy.txt)
        let pow = 0.0, dis = 0.0;
        for (let s = 0; s < 4; s++) {
            const om = oms[s];
            let p = 0.0, d = 0.0;
            for (let i = 0; i < n; i++) p += om[i] * tau[s][i];
            for (let i = 1; i < n; i++) d += (om[i] - om[i - 1]) ** 2 / this.ll[i];
            pow += bw[s] * p;
            dis += bw[s] * d;
        }
        this.pow = pow;
        this.dis = nu * dis;

        this.time += dt;
        this.dt = dt;
    }

    // ------------------------------------------------------------ |
    // kinetic energy of each mass, 1/2 mo(i) |v(i)|^2              |
    // (Spectra/kinspc2.txt); returns the total                     |

    energySpectrum(out) {
        const { n, ll, mo, theta, omega } = this;
        let vx = 0.0, vy = 0.0, tot = 0.0;
        for (let i = 0; i < n; i++) {
            vx += ll[i] * omega[i] * Math.cos(theta[i]);
            vy += ll[i] * omega[i] * Math.sin(theta[i]);
            out[i] = 0.5 * mo[i] * (vx * vx + vy * vy);
            tot += out[i];
        }
        return tot;
    }

    // gravitational energy (energy.f90, sw = 2)
    potential() {
        const { n, ll, mm, theta, gr } = this;
        let e = gr * mm[0] * 2.0;
        for (let i = 0; i < n; i++) e -= gr * mm[i] * ll[i] * Math.cos(theta[i]);
        return e;
    }

    // apply an angular-velocity kick dw to link k
    kick(k, dw) {
        this.omega[k] += dw;
    }

    isFinite() {
        for (let i = 0; i < this.n; i++) {
            if (!Number.isFinite(this.theta[i]) || !Number.isFinite(this.omega[i])) return false;
        }
        return true;
    }
}
