/* The physics of the Dust Explorer: a port of dustlib.py, which documents the sources and
   the choices.  validate.py checks this file against the Python on a set of cases
   (testjs.js runs it under node).

   Units as in dustlib: wavelength in um, u_lambda in erg cm^-3 um^-1, powers in erg s^-1,
   cross sections in cm^2, energies in erg. */
(function (root) {
"use strict";

const H = 6.62607015e-27, C = 2.99792458e10, KB = 1.380649e-16, EV = 1.602176634e-12;
const UM = 1e-4, LYMAN_UM = 0.0912;
const XMAX = 80;          // h nu / k T beyond which a Planck term is dropped (< e^-70 of the peak)
let ETOP = H * C / (LYMAN_UM * UM);            // the most energetic photon on the grid

const Phys = { H, C, KB, EV, UM, ready: false };
let NL = 0, LAM, WL, EPH, XT, B0, USTAR, UCMB, DLNLAM, LNE0;
let AG, NA, AGLO, AGDL;                      // size grid
let TGE, LTGE, DLTGE, ESIL, EGRA;            // enthalpy tables
let MSIL, EPS = {}, QT = {}, MATS;

// ------------------------------------------------------------------ data
Phys.init = function (meta, buf) {
  const get = name => {
    const a = meta.arrays.find(x => x.name === name);
    return a.dtype === "f8" ? new Float64Array(buf, a.offset, a.length)
                            : new Float32Array(buf, a.offset, a.length);
  };
  LAM = get("lam"); NL = LAM.length;
  DLNLAM = Math.log(LAM[1] / LAM[0]);
  WL = new Float64Array(NL); EPH = new Float64Array(NL); XT = new Float64Array(NL);
  B0 = new Float64Array(NL); UCMB = new Float64Array(NL);
  for (let i = 0; i < NL; i++) {
    const w = (i === 0 || i === NL - 1) ? 0.5 * DLNLAM : DLNLAM;
    WL[i] = w * LAM[i];                                   // trapezoid in ln lambda
    const l = LAM[i] * UM;
    EPH[i] = H * C / l;
    XT[i] = H * C / (l * KB);                             // h nu / k
    B0[i] = 2 * H * C * C / Math.pow(l, 5) * UM;          // B_lambda = B0 / expm1(XT / T)
  }
  for (let i = 0; i < NL; i++) UCMB[i] = 4 * Math.PI * planck(i, meta.T_CMB) / C;
  LNE0 = Math.log(EPH[NL - 1]);                           // smallest photon energy
  ETOP = EPH[0];
  USTAR = get("ustar");
  MSIL = get("m_sil");
  for (const k of ["pa0.01", "pa0.10", "pe0.01", "pe0.10"]) EPS[k] = get("eps_" + k);
  QT.sil = get("q_sil"); QT.gra = get("q_gra");
  NA = meta.a_grid.n; AGLO = meta.a_grid.lo; AGDL = (meta.a_grid.hi - meta.a_grid.lo) / (NA - 1);
  AG = new Float64Array(NA);
  for (let i = 0; i < NA; i++) AG[i] = Math.pow(10, AGLO + i * AGDL);
  TGE = get("tgrid_e"); ESIL = get("e_sil"); EGRA = get("e_gra");
  LTGE = new Float64Array(TGE.length);
  for (let i = 0; i < TGE.length; i++) LTGE[i] = Math.log(TGE[i]);
  DLTGE = (LTGE[LTGE.length - 1] - LTGE[0]) / (LTGE.length - 1);
  MATS = meta.mat;
  Phys.T_CMB = meta.T_CMB;
  Phys.MRN = meta.mrn;
  Phys.LAM = LAM; Phys.NL = NL; Phys.A_GRID = AG; Phys.WL = WL; Phys.USTAR = USTAR; Phys.UCMB = UCMB;
  Phys.ready = true;
};

function planck(i, T) {
  const x = XT[i] / T;
  return x < XMAX ? B0[i] / Math.expm1(x) : 0;
}
Phys.planck = planck;
function integ(F) { let s = 0; for (let i = 0; i < NL; i++) s += WL[i] * F[i]; return s; }
Phys.integ = integ;

// ------------------------------------------------------------------ materials
function nAtoms(mat, a) { return MATS[mat].atoms_per_cm3 * 4 / 3 * Math.PI * Math.pow(a * UM, 3); }
function grainMass(mat, a) { return MATS[mat].rho * 4 / 3 * Math.PI * Math.pow(a * UM, 3); }
Phys.nAtoms = nAtoms; Phys.grainMass = grainMass;

// principal square root of a complex number
function csqrt(re, im) {
  const r = Math.hypot(re, im);
  if (re >= 0) { const t = Math.sqrt(0.5 * (r + re)); return [t, t > 0 ? im / (2 * t) : 0]; }
  const t = Math.sqrt(0.5 * (r - re));
  return [Math.abs(im) / (2 * t), im >= 0 ? t : -t];
}

// refractive index on the grid: [re, im] arrays; graphite gives two (E par, E perp to c)
function refrIndex(mat, a) {
  if (mat === "sil") {
    const re = new Float64Array(NL), im = new Float64Array(NL);
    for (let i = 0; i < NL; i++) { re[i] = MSIL[2 * i]; im[i] = MSIL[2 * i + 1]; }
    return [[re, im]];
  }
  const w = Math.min(1, Math.max(0, Math.log10(a / 0.01)));
  return ["pa", "pe"].map(o => {
    const e1 = EPS[o + "0.01"], e2 = EPS[o + "0.10"];
    const re = new Float64Array(NL), im = new Float64Array(NL);
    for (let i = 0; i < NL; i++) {
      const m = csqrt((1 - w) * e1[2 * i] + w * e2[2 * i], (1 - w) * e1[2 * i + 1] + w * e2[2 * i + 1]);
      re[i] = m[0]; im[i] = m[1];
    }
    return [re, im];
  });
}

// ------------------------------------------------------------------ Mie theory (BHMIE)
let DBr = new Float64Array(64), DBi = new Float64Array(64);
function mie(x, mr, mi) {
  const nstop = Math.floor(x + 4 * Math.cbrt(x) + 2);
  const yr = mr * x, yi = mi * x, y2 = yr * yr + yi * yi;
  const nmx = Math.floor(Math.max(nstop, Math.sqrt(y2))) + 15;
  if (DBr.length < nstop + 2) { DBr = new Float64Array(2 * nstop + 2); DBi = new Float64Array(2 * nstop + 2); }
  let dr = 0, di = 0;
  for (let n = nmx; n >= 1; n--) {               // D_{n-1} = n/y - 1/(D_n + n/y)
    const nr = n * yr / y2, ni = -n * yi / y2;
    const tr = dr + nr, ti = di + ni, den = tr * tr + ti * ti;
    dr = nr - tr / den; di = ni + ti / den;
    if (n - 1 <= nstop) { DBr[n - 1] = dr; DBi[n - 1] = di; }
  }
  const m2 = mr * mr + mi * mi;
  let psi0 = Math.cos(x), psi1 = Math.sin(x), chi0 = -Math.sin(x), chi1 = Math.cos(x);
  let xi1r = psi1, xi1i = -chi1;
  let qs = 0, qe = 0, gs = 0, a1r = 0, a1i = 0, b1r = 0, b1i = 0;
  for (let n = 1; n <= nstop; n++) {
    const psi = (2 * n - 1) * psi1 / x - psi0;
    const chi = (2 * n - 1) * chi1 / x - chi0;
    const xir = psi, xii = -chi;
    const Dr = DBr[n], Di = DBi[n];
    // ta = D/m + n/x, tb = m D + n/x
    const tar = (Dr * mr + Di * mi) / m2 + n / x, tai = (Di * mr - Dr * mi) / m2;
    const tbr = mr * Dr - mi * Di + n / x, tbi = mr * Di + mi * Dr;
    // a_n = (ta psi - psi1) / (ta xi - xi1)
    let nr = tar * psi - psi1, ni = tai * psi;
    let er = tar * xir - tai * xii - xi1r, ei = tar * xii + tai * xir - xi1i, dd = er * er + ei * ei;
    const anr = (nr * er + ni * ei) / dd, ani = (ni * er - nr * ei) / dd;
    nr = tbr * psi - psi1; ni = tbi * psi;
    er = tbr * xir - tbi * xii - xi1r; ei = tbr * xii + tbi * xir - xi1i; dd = er * er + ei * ei;
    const bnr = (nr * er + ni * ei) / dd, bni = (ni * er - nr * ei) / dd;
    qs += (2 * n + 1) * (anr * anr + ani * ani + bnr * bnr + bni * bni);
    qe += (2 * n + 1) * (anr + bnr);
    gs += (2 * n + 1) / (n * (n + 1)) * (anr * bnr + ani * bni);
    if (n > 1) gs += (n - 1) * (n + 1) / n * (a1r * anr + a1i * ani + b1r * bnr + b1i * bni);
    psi0 = psi1; psi1 = psi; chi0 = chi1; chi1 = chi; xi1r = psi1; xi1i = -chi1;
    a1r = anr; a1i = ani; b1r = bnr; b1i = bni;
  }
  const f = 2 / (x * x);
  return [f * qe, f * qs, x < 1e-3 ? 0 : 2 * gs / qs];
}
Phys.mie = mie;

// Q_abs, Q_sca, g on the grid for one grain, by Mie theory
function optics(mat, a) {
  const ms = refrIndex(mat, a);
  const Qa = new Float64Array(NL), Qs = new Float64Array(NL), g = new Float64Array(NL);
  for (let i = 0; i < NL; i++) {
    const x = 2 * Math.PI * a / LAM[i];
    if (mat === "sil") {
      const r = mie(x, ms[0][0][i], ms[0][1][i]);
      Qa[i] = r[0] - r[1]; Qs[i] = r[1]; g[i] = r[2];
    } else {
      const p = mie(x, ms[0][0][i], ms[0][1][i]), q = mie(x, ms[1][0][i], ms[1][1][i]);
      const qe = p[0] / 3 + 2 * q[0] / 3, qs = p[1] / 3 + 2 * q[1] / 3;
      Qa[i] = qe - qs; Qs[i] = qs; g[i] = (p[1] * p[2] / 3 + 2 * q[1] * q[2] / 3) / qs;
    }
  }
  return { Qabs: Qa, Qsca: Qs, g };
}
Phys.optics = optics;

// Q_abs and Q_sca from the table at grid node i (float32, as shipped)
function tableOptics(mat, i) {
  const T = QT[mat], Qa = new Float64Array(NL), Qs = new Float64Array(NL);
  for (let k = 0; k < NL; k++) { Qa[k] = T[i * NL + k]; Qs[k] = T[(NA + i) * NL + k]; }
  return { Qabs: Qa, Qsca: Qs, g: null };
}
Phys.tableOptics = tableOptics;

// ------------------------------------------------------------------ enthalpy
function Enthalpy(mat, a) {
  this.N = nAtoms(mat, a);
  const e = mat === "sil" ? ESIL : EGRA, n = e.length;
  this.lE = new Float64Array(n);
  for (let i = 0; i < n; i++) this.lE[i] = Math.log((this.N - 2) * e[i]);
}
Enthalpy.prototype.E = function (T) {
  const n = LTGE.length, u = (Math.log(T) - LTGE[0]) / DLTGE;
  if (u <= 0) return Math.exp(this.lE[0]);
  if (u >= n - 1) return Math.exp(this.lE[n - 1]);
  // the grid is uniform in ln T; locate by index and confirm against the stored nodes
  let k = Math.floor(u);
  const lt = Math.log(T);
  while (k > 0 && LTGE[k] > lt) k--;
  while (k < n - 2 && LTGE[k + 1] < lt) k++;
  const t = (lt - LTGE[k]) / (LTGE[k + 1] - LTGE[k]);
  return Math.exp(this.lE[k] + t * (this.lE[k + 1] - this.lE[k]));
};
Enthalpy.prototype.T = function (E) {
  const lE = this.lE, n = lE.length, le = Math.log(E);
  if (le <= lE[0]) return TGE[0];
  if (le >= lE[n - 1]) return TGE[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (lE[m] > le) hi = m; else lo = m; }
  const t = (le - lE[lo]) / (lE[hi] - lE[lo]);
  return Math.exp(LTGE[lo] + t * (LTGE[hi] - LTGE[lo]));
};
Enthalpy.prototype.C = function (T) {
  const d = 1e-4;
  return (this.E(T * (1 + d)) - this.E(T * (1 - d))) / (2 * d * T);
};
Phys.Enthalpy = Enthalpy;

// ------------------------------------------------------------------ one grain in a field
const SIGMA_EQ = 0.01, M0 = 400, K_COOL = 30, DELTA = 0.03, M_MAX = 1500;
const TRIM = 1e-14, TRIM_POW = 1e-10;
Phys.SIGMA_EQ = SIGMA_EQ;

function Grain(mat, a, U, opt) {
  this.mat = mat; this.a = a; this.U = U;
  const o = opt || optics(mat, a);
  this.Qabs = o.Qabs; this.Qsca = o.Qsca; this.g = o.g;
  const geo = Math.PI * Math.pow(a * UM, 2);
  this.geo = geo;
  this.Cabs = new Float64Array(NL);
  this.u = new Float64Array(NL); this.us = new Float64Array(NL);
  this.pabsLam = new Float64Array(NL);
  for (let i = 0; i < NL; i++) {
    this.Cabs[i] = geo * this.Qabs[i];
    this.us[i] = U * USTAR[i];
    this.u[i] = this.us[i] + UCMB[i];
    this.pabsLam[i] = C * this.Cabs[i] * this.u[i];
  }
  this.mass = grainMass(mat, a);
  this.H = new Enthalpy(mat, a);
  this.Pabs = integ(this.pabsLam);
  this.Teq = this.equilibriumT();
}

// emitted power at T; with Emax, only photons less energetic than Emax
Grain.prototype.L = function (T, Emax) {
  let s = 0;
  const lim = Emax === undefined ? Infinity : Emax;
  for (let i = 0; i < NL; i++) {
    if (EPH[i] >= lim) continue;
    const x = XT[i] / T;
    if (x >= XMAX) continue;
    s += WL[i] * this.Cabs[i] * B0[i] / Math.expm1(x);
  }
  return 4 * Math.PI * s;
};

Grain.prototype.equilibriumT = function () {
  let lo = 0, hi = Math.log(5000);
  for (let k = 0; k < 64; k++) {
    const mid = 0.5 * (lo + hi);
    if (this.L(Math.exp(mid)) > this.Pabs) hi = mid; else lo = mid;
  }
  return Math.exp(0.5 * (lo + hi));
};

// starlight photon statistics: absorption rate, power, and int (h nu) dP
Grain.prototype.photonStats = function () {
  if (this._ps) return this._ps;
  let R = 0, S1 = 0, S2 = 0;
  for (let i = 0; i < NL; i++) {
    const p = C * this.Cabs[i] * this.us[i];
    R += WL[i] * p / EPH[i]; S1 += WL[i] * p; S2 += WL[i] * p * EPH[i];
  }
  return (this._ps = [R, S1, S2]);
};

Grain.prototype.sigmaT = function () {
  const [, , S2] = this.photonStats(), T = this.Teq;
  const dL = (this.L(T * 1.001) - this.L(T * 0.999)) / (0.002 * T);
  const Cv = this.H.C(T), tau = Cv / dL, sE = Math.sqrt(0.5 * tau * S2);
  return [sE / Cv / T, sE];
};

// cumulative int f dE and int E f dE, f = C_abs u_E, on the ascending photon-energy grid
Grain.prototype.cumulants = function () {
  if (this._cum) return this._cum;
  const I0 = new Float64Array(NL), I1 = new Float64Array(NL);
  let pg0 = 0, pg1 = 0;
  for (let k = 0; k < NL; k++) {
    const i = NL - 1 - k, E = EPH[i];
    const f = this.Cabs[i] * this.u[i] * LAM[i] / E;
    const g0 = f * E, g1 = f * E * E;
    if (k > 0) { I0[k] = I0[k - 1] + 0.5 * (g0 + pg0) * DLNLAM; I1[k] = I1[k - 1] + 0.5 * (g1 + pg1) * DLNLAM; }
    pg0 = g0; pg1 = g1;
  }
  return (this._cum = [I0, I1]);
};

Grain.prototype.stochastic = function (Mfirst) {
  const M = Mfirst || M0, Hh = this.H;
  const Teq = this.Teq, Eeq = Hh.E(Teq);
  const sE = this.sigmaT()[1];
  const R = this.photonStats()[0];
  const Lcmb = this.L(Phys.T_CMB);
  // bottom of the first range: cool from T_eq for K_COOL photon intervals, or 10 sigma
  const nT = 200, lt0 = Math.log(Teq), lt1 = Math.log(Phys.T_CMB * 1.0005);
  let tc = 0, Tcool = Math.exp(lt1), prev = null;
  for (let k = 0; k < nT; k++) {
    const T = Math.exp(lt0 + (lt1 - lt0) * k / (nT - 1));
    const v = Hh.C(T) * T / Math.max(this.L(T) - Lcmb, 1e-300);
    if (prev !== null) {
      tc += (lt0 - lt1) / (nT - 1) * 0.5 * (v + prev[1]);
      if (tc >= K_COOL / R) { Tcool = T; break; }
    }
    prev = [T, v];
  }
  const Elo = Math.max(Hh.E(Phys.T_CMB), Math.min(Eeq - 10 * sE, Hh.E(Tcool)));
  let Ehi = Eeq + 13.6 * EV + 10 * sE;
  let res, it;
  for (it = 0; it < 8; it++) {
    const Tlo = Hh.T(Elo), Thi = Math.min(Hh.T(Ehi), 4900);
    res = this.solve(logGrid(Tlo, Thi, M));
    if (res.P[M] < 1e-10 || Thi >= 4900) break;
    Ehi = Elo + 1.5 * (Ehi - Elo);
  }
  const [jl, ju] = trim(res.P, res.down);
  const T2 = this.secondGrid(res.T, res.P, res.E, jl, ju, M);
  const out = this.solve(T2);
  out.iterations = it + 2;
  return out;
};

function logGrid(lo, hi, M) {
  const T = new Float64Array(M + 1), a = Math.log(lo), b = Math.log(hi);
  for (let j = 0; j <= M; j++) T[j] = Math.exp(a + (b - a) * j / M);
  return T;
}

function trim(P, down) {
  const M = P.length - 1;
  let pw = 0;
  for (let j = 0; j <= M; j++) pw += P[j] * down[j];
  let jl = 0, c = 0;
  for (let j = 0; j <= M; j++) { c += P[j]; if (c > TRIM) { jl = Math.max(0, j - 1); break; } }
  let hi = 0, hp = 0, ju = 0;
  for (let j = M; j >= 0; j--) {
    hi += P[j]; hp += P[j] * down[j];
    if (hi > TRIM || hp / pw > TRIM_POW) { ju = Math.min(M, j + 1); break; }
  }
  return [jl, ju];
}

Grain.prototype.secondGrid = function (T, P, E, jl, ju, M) {
  const [, S1, S2] = this.photonStats();
  const eps = S2 / S1, n = E.length;
  let k = 0, best = -1;
  for (let j = 0; j < n; j++) {
    // np.gradient: one-sided at the ends, centered inside
    const g = j === 0 ? Math.log(E[1] / E[0]) : j === n - 1 ? Math.log(E[n - 1] / E[n - 2])
                      : 0.5 * Math.log(E[j + 1] / E[j - 1]);
    const v = P[j] / g;
    if (v > best) { best = v; k = j; }
  }
  const Eref = Math.max(E[k], eps);
  const Tref = this.H.T(Eref);
  const slope = Math.log(this.H.E(Tref * 1.001) / this.H.E(Tref / 1.001)) / Math.log(1.001 * 1.001);
  const dlnT = DELTA * eps / (Eref * slope);
  const lo = Math.log(T[jl]), hi = Math.log(T[ju]);
  const nb = Math.min(M_MAX, Math.max(M, Math.ceil((hi - lo) / dlnT)));
  return logGrid(T[jl], T[ju], nb);
};

// The thermal-continuous steady state on the temperature levels T (DL01 Sects. 4, 7.2, 9.2).
// Upward rates are banded (no photon carries more than ETOP) apart from the row into the
// top bin, which is handled as a running sum.
Grain.prototype.solve = function (T) {
  const Hh = this.H, M = T.length - 1;
  const E = new Float64Array(M + 1);
  for (let j = 0; j <= M; j++) E[j] = Hh.E(T[j]);
  const Emin = new Float64Array(M + 1), Emax = new Float64Array(M + 1), dE = new Float64Array(M + 1);
  for (let j = 0; j <= M; j++) {
    Emin[j] = j === 0 ? Math.max(0, E[0] - 0.5 * (E[1] - E[0])) : 0.5 * (E[j] + E[j - 1]);
    Emax[j] = j === M ? E[M] + 0.5 * (E[M] - E[M - 1]) : 0.5 * (E[j] + E[j + 1]);
    dE[j] = Emax[j] - Emin[j];
  }
  const [I0, I1] = this.cumulants();
  const I0top = I0[NL - 1], I1top = I1[NL - 1];
  // ascending photon energies; the cumulants are interpolated linearly in E
  const EgA = this._egA || (this._egA = Float64Array.from({ length: NL }, (_, k) => EPH[NL - 1 - k]));
  const Eg0 = EgA[0], EgTop = EgA[NL - 1];
  const at = (I, k, w) => I[k] + (w - EgA[k]) / (EgA[k + 1] - EgA[k]) * (I[k + 1] - I[k]);
  const locate = w => {                    // k with EgA[k] <= w < EgA[k+1]
    let lo = 0, hi = NL - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (EgA[m] > w) hi = m; else lo = m; }
    return lo;
  };
  const cum = (I, w) => (w <= Eg0 ? 0 : w >= EgTop ? I[NL - 1] : at(I, locate(w), w));
  // band: column l holds A[u, l] for u = l+1 .. ulim[l] (all < M); top[l] = A[M, l]
  const ulim = new Int32Array(M + 1), off = new Int32Array(M + 2);
  let up = 1, wid = 0;
  for (let l = 0; l < M; l++) {
    if (up < l + 1) up = l + 1;
    while (up + 1 < M && Emin[up + 1] - Emax[l] < ETOP) up++;
    ulim[l] = Math.min(up, M - 1);
    off[l + 1] = off[l] + (ulim[l] - l);
    wid = Math.max(wid, ulim[l] - l);
  }
  const band = new Float64Array(off[M] + 1);
  const top = new Float64Array(M + 1);
  // cumulants at the two monotone sequences each column needs:
  // D1[k] = Emin[u] - Emax[l] and D2[k] = Emin[u] - Emin[l], u = l + 1 + k
  const D1 = new Float64Array(wid + 2), D2 = new Float64Array(wid + 2);
  const c01 = new Float64Array(wid + 2), c11 = new Float64Array(wid + 2);
  const c02 = new Float64Array(wid + 2), c12 = new Float64Array(wid + 2);
  const walk = (D, n, o0, o1) => {
    let k = 0;
    for (let i = 0; i < n; i++) {
      const w = D[i];
      if (w <= Eg0) { o0[i] = 0; o1[i] = 0; continue; }
      if (w >= EgTop) { o0[i] = I0top; o1[i] = I1top; continue; }
      while (EgA[k + 1] <= w) k++;
      o0[i] = at(I0, k, w); o1[i] = at(I1, k, w);
    }
  };
  for (let l = 0; l < M; l++) {
    const nu = ulim[l] - l;
    for (let k = 0; k <= nu; k++) { const u = l + 1 + k; D1[k] = Emin[u] - Emax[l]; D2[k] = Emin[u] - Emin[l]; }
    walk(D1, nu + 1, c01, c11);
    walk(D2, nu + 1, c02, c12);
    for (let k = 0; k < nu; k++) {
      const u = l + 1 + k;
      // W1 = D1[k], a2 = D1[k+1], a1 = D2[k], W4 = D2[k+1]; W2/W3 = min/max(a1, a2)
      const W1 = D1[k], W4 = D2[k + 1];
      let W2, c0W2, c1W2, c0W3, c1W3;
      if (D2[k] <= D1[k + 1]) { W2 = D2[k]; c0W2 = c02[k]; c1W2 = c12[k]; c0W3 = c01[k + 1]; c1W3 = c11[k + 1]; }
      else { W2 = D1[k + 1]; c0W2 = c01[k + 1]; c1W2 = c11[k + 1]; c0W3 = c02[k]; c1W3 = c12[k]; }
      const G = ((c1W2 - c11[k]) - W1 * (c0W2 - c01[k]) + Math.min(dE[l], dE[u]) * (c0W3 - c0W2)
                 + W4 * (c02[k + 1] - c0W3) - (c12[k + 1] - c1W3)) / (dE[u] * dE[l]);
      band[off[l] + k] = C * dE[u] / (E[u] - E[l]) * G;
    }
    const W1 = Emin[M] - Emax[l], Wc = Emin[M] - Emin[l], W1c = Math.max(W1, 0);
    const ramp = (cum(I1, Wc) - cum(I1, W1c)) - W1 * (cum(I0, Wc) - cum(I0, W1c));
    top[l] = C / (E[M] - E[l]) * (ramp / (Wc - W1) + (I0top - cum(I0, Wc)));
    // intrabin absorptions in bin l promote to l + 1
    const intra = C / (E[l + 1] - E[l]) * (dE[l] * cum(I0, dE[l]) - cum(I1, dE[l])) / dE[l];
    if (l + 1 === M) top[l] += intra; else band[off[l]] += intra;
  }
  // suffix sums down each column: S[j, l] = sum_{u = j}^{ulim[l]} A[u, l]
  for (let l = 0; l < M; l++) {
    for (let k = off[l + 1] - 2; k >= off[l]; k--) band[k] += band[k + 1];
  }
  // continuous cooling
  const down = new Float64Array(M + 1), Dn = new Float64Array(M + 1);
  for (let u = 0; u <= M; u++) down[u] = this.L(T[u], E[u]);
  for (let u = 1; u <= M; u++) Dn[u] = down[u] / (E[u] - E[u - 1]);
  // Guhathakurta & Draine recursion; lmin[j] = smallest l whose band reaches j
  const X = new Float64Array(M + 1);
  X[0] = 1;
  let topAcc = top[0] * X[0], lmin = 0;
  for (let j = 1; j <= M; j++) {
    while (lmin < j && ulim[lmin] < j) lmin++;
    let s = topAcc;
    if (j < M) for (let l = lmin; l < j; l++) if (ulim[l] >= j) s += band[off[l] + (j - l - 1)] * X[l];
    X[j] = s / Dn[j];
    if (j < M) topAcc += top[j] * X[j];
    if (X[j] > 1e250) {
      for (let i = 0; i <= j; i++) X[i] *= 1e-250;
      topAcc *= 1e-250;
    }
  }
  let sum = 0;
  for (let j = 0; j <= M; j++) sum += X[j];
  const P = new Float64Array(M + 1);
  let pem = 0;
  for (let j = 0; j <= M; j++) { P[j] = X[j] / sum; pem += P[j] * down[j]; }
  return { T, P, E, down, Pemit: pem };
};

// emitted power per um per grain; also the temperatures and probabilities behind it
Grain.prototype.emission = function () {
  const em = new Float64Array(NL);
  const sT = this.sigmaT()[0];
  if (sT < SIGMA_EQ) {
    for (let i = 0; i < NL; i++) em[i] = 4 * Math.PI * this.Cabs[i] * planck(i, this.Teq);
    return { em, T: [this.Teq], P: [1], stochastic: false, sigmaT: sT };
  }
  const r = this.stochastic();
  for (let j = 0; j < r.T.length; j++) {
    const Pj = r.P[j];
    if (Pj < 1e-30) continue;
    const Tj = r.T[j], Ej = r.E[j];
    for (let i = 0; i < NL; i++) {
      if (EPH[i] >= Ej) continue;
      const x = XT[i] / Tj;
      if (x >= XMAX) continue;
      em[i] += Pj * B0[i] / Math.expm1(x);
    }
  }
  for (let i = 0; i < NL; i++) em[i] *= 4 * Math.PI * this.Cabs[i];
  return { em, T: r.T, P: r.P, E: r.E, stochastic: true, sigmaT: sT, Pemit: r.Pemit };
};
Phys.Grain = Grain;

// ------------------------------------------------------------------ size distributions
const NSUB = 8;
function powerlawWeights(mat, q, amin, amax, abund) {
  const step = Math.LN10 * AGDL;
  const n = Math.max(Math.ceil((Math.log(amax) - Math.log(amin)) / step * NSUB) + 1, 2);
  const a = new Float64Array(n), w = new Float64Array(n);
  const l0 = Math.log(amin), dl = (Math.log(amax) - l0) / (n - 1);
  let norm = 0;
  const per = MATS[mat].per_key;
  for (let k = 0; k < n; k++) {
    a[k] = Math.exp(l0 + k * dl);
    const tw = (k === 0 || k === n - 1) ? 0.5 * dl : dl;
    w[k] = Math.pow(a[k], 1 - q) * tw;
    norm += w[k] * nAtoms(mat, a[k]) / per;
  }
  const A = abund * 1e-6 / norm;
  for (let k = 0; k < n; k++) w[k] *= A;
  return { a, w };
}
Phys.powerlawWeights = powerlawWeights;

// fractional index of radius a on the size grid
function aIndex(a) {
  const u = (Math.log10(a) - AGLO) / AGDL;
  let i = Math.floor(u);
  i = Math.max(0, Math.min(NA - 2, i));
  // match numpy's searchsorted on the stored nodes at exact boundaries
  const la = Math.log(a);
  while (i > 0 && Math.log(AG[i]) >= la) i--;
  while (i < NA - 2 && Math.log(AG[i + 1]) < la) i++;
  return [i, (la - Math.log(AG[i])) / (Math.log(AG[i + 1]) - Math.log(AG[i]))];
}

const SIZE_BINS = [[0, 0.01], [0.01, 0.1], [0.1, Infinity]];
Phys.SIZE_BINS = SIZE_BINS;

// Optical depths per H (cm^2) and mass per H for one material
function distOptics(mat, p) {
  const { a, w } = powerlawWeights(mat, p.q, p.amin, p.amax, p.abund);
  const ext = new Float64Array(NL), sca = new Float64Array(NL), T = QT[mat];
  let mass = 0;
  for (let k = 0; k < a.length; k++) {
    const [i, t] = aIndex(a[k]);
    const geo = w[k] * Math.PI * Math.pow(a[k] * UM, 2);
    for (let j = 0; j < NL; j++) {
      const qa = (1 - t) * T[i * NL + j] + t * T[(i + 1) * NL + j];
      const qs = (1 - t) * T[(NA + i) * NL + j] + t * T[(NA + i + 1) * NL + j];
      ext[j] += geo * (qa + qs); sca[j] += geo * qs;
    }
    mass += w[k] * grainMass(mat, a[k]);
  }
  return { ext, sca, mass, a, w };
}
Phys.distOptics = distOptics;

// nodes of the size grid needed to interpolate over [amin, amax]
function nodeRange(amin, amax) {
  let i0 = 0, i1 = NA - 1;
  for (let i = 0; i < NA; i++) if (AG[i] < amin) i0 = i;          // searchsorted(amin) - 1
  for (let i = NA - 1; i >= 0; i--) if (AG[i] >= amax) i1 = i;     // searchsorted(amax)
  if (AG[0] >= amin) i0 = 0;
  return [i0, Math.min(NA - 1, i1)];
}
Phys.nodeRange = nodeRange;

// emission of one grain at node i of the size grid, from the table optics, with a
// summary of its temperatures: T_eq, the 5/50/95th percentiles of P(T), and the
// temperature after absorbing one photon of the mean absorbed energy from T_eq
function nodeEmission(mat, i, U) {
  const g = new Grain(mat, AG[i], U, tableOptics(mat, i));
  const r = g.emission();
  const ps = g.photonStats();
  const out = { em: r.em, pabs: g.Pabs, stochastic: r.stochastic, Teq: g.Teq,
                Tjump: g.H.T(g.H.E(g.Teq) + ps[1] / ps[0]) };
  [out.t05, out.t50, out.t95] = percentiles(r.T, r.P, [0.05, 0.5, 0.95]);
  return out;
}
Phys.nodeEmission = nodeEmission;

// percentiles of a distribution over levels T (interpolated in ln T within a level)
function percentiles(T, P, qs) {
  if (T.length === 1) return qs.map(() => T[0]);
  const out = [];
  let c = 0, j = 0;
  for (const q of qs) {
    while (j < T.length - 1 && c + P[j] < q) { c += P[j]; j++; }
    const f = P[j] > 0 ? Math.min(1, Math.max(0, (q - c) / P[j])) : 0;
    const lo = j > 0 ? Math.sqrt(T[j - 1] * T[j]) : T[0];
    const hi = j < T.length - 1 ? Math.sqrt(T[j] * T[j + 1]) : T[j];
    out.push(lo * Math.pow(hi / lo, f));
  }
  return out;
}
Phys.percentiles = percentiles;

// dn/da = A a^-q (a in um, per um per H) for the normalization of powerlawWeights
function powerlawA(mat, q, amin, amax, abund) {
  const { a, w } = powerlawWeights(mat, q, amin, amax, abund);
  const dl = Math.log(a[1] / a[0]);
  return w[1] / (Math.pow(a[1], 1 - q) * dl);
}
Phys.powerlawA = powerlawA;

// Emission per H for one material, given the per-node emission (array indexed by node)
function distEmission(mat, p, nodes) {
  const { a, w } = powerlawWeights(mat, p.q, p.amin, p.amax, p.abund);
  const em = new Float64Array(NL);
  const bins = SIZE_BINS.map(() => new Float64Array(NL));
  let pabs = 0;
  for (let k = 0; k < a.length; k++) {
    const [i, t] = aIndex(a[k]);
    const n0 = nodes[i], n1 = nodes[i + 1];
    const s0 = (1 - t) / Math.pow(AG[i], 3), s1 = t / Math.pow(AG[i + 1], 3), a3 = Math.pow(a[k], 3);
    const b = SIZE_BINS.findIndex(r => a[k] >= r[0] && a[k] < r[1]);
    for (let j = 0; j < NL; j++) {
      const v = w[k] * a3 * (s0 * n0.em[j] + s1 * n1.em[j]);
      em[j] += v; bins[b][j] += v;
    }
    pabs += w[k] * a3 * (s0 * n0.pabs + s1 * n1.pabs);
  }
  for (let j = 0; j < NL; j++) {
    const f = LAM[j] / (4 * Math.PI);
    em[j] *= f;
    for (const bn of bins) bn[j] *= f;
  }
  return { em, bins, pabs };
}
Phys.distEmission = distEmission;

if (typeof module !== "undefined" && module.exports) module.exports = Phys;
else root.Phys = Phys;
})(typeof window !== "undefined" ? window : this);
