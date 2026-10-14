/* The physics of the Dust Explorer: a port of dustlib.py, which documents the sources and
   the choices.  validate.py checks this file against the Python on a set of cases
   (testjs.js runs it under node).

   Two wavelength grids: the coarse one (80 per decade) for silicate and graphite, and a fine
   one (500 per decade over 3-20 um) for PAHs and for the size-distribution sums.

   Units as in dustlib: wavelength in um, u_lambda in erg cm^-3 um^-1, powers in erg s^-1,
   cross sections in cm^2, energies in erg. */
(function (root) {
"use strict";

const H = 6.62607015e-27, C = 2.99792458e10, KB = 1.380649e-16, EV = 1.602176634e-12;
const AMU = 1.66053906660e-24, UM = 1e-4;
const XMAX = 80;          // h nu / k T beyond which a Planck term is dropped (< e^-70 of the peak)

const Phys = { H, C, KB, EV, UM, ready: false };
let GC = null, GF = null;                    // coarse and fine grids
let AG, NA, AGLO, AGDL;                      // size grid for silicate and graphite
let TGE, LTGE, DLTGE, ESIL, EGRA;            // enthalpy tables
let MSIL, EPS = {}, QT = {}, MATS, PAG;
let TOF_I = null, TOF_T = null;              // coarse -> fine interpolation

// ------------------------------------------------------------------ grids
// A wavelength grid with trapezoid weights in ln(lambda): int F dlambda = sum(wl * F)
function makeGrid(lam, ustar, Tcmb) {
  const n = lam.length, G = { n, lam, ustar };
  G.wl = new Float64Array(n); G.eph = new Float64Array(n); G.xt = new Float64Array(n);
  G.b0 = new Float64Array(n); G.ucmb = new Float64Array(n); G.lnl = new Float64Array(n);
  for (let i = 0; i < n; i++) G.lnl[i] = Math.log(lam[i]);
  for (let i = 0; i < n; i++) {
    let w = 0;
    if (i > 0) w += 0.5 * (G.lnl[i] - G.lnl[i - 1]);
    if (i < n - 1) w += 0.5 * (G.lnl[i + 1] - G.lnl[i]);
    G.wl[i] = w * lam[i];
    const l = lam[i] * UM;
    G.eph[i] = H * C / l;
    G.xt[i] = H * C / (l * KB);                           // h nu / k
    G.b0[i] = 2 * H * C * C / Math.pow(l, 5) * UM;        // B_lambda = b0 / expm1(xt / T)
  }
  for (let i = 0; i < n; i++) G.ucmb[i] = 4 * Math.PI * planck(G, i, Tcmb) / C;
  G.egA = Float64Array.from({ length: n }, (_, k) => G.eph[n - 1 - k]);   // ascending
  G.etop = G.eph[0];
  return G;
}
function planck(G, i, T) {
  const x = G.xt[i] / T;
  return x < XMAX ? G.b0[i] / Math.expm1(x) : 0;
}
function integ(G, F) { let s = 0; for (let i = 0; i < G.n; i++) s += G.wl[i] * F[i]; return s; }
Phys.planck = planck; Phys.integ = integ;

// a coarse-grid array onto the fine grid, linear in ln(lambda) (as numpy's interp)
function toFine(arr) {
  const out = new Float64Array(GF.n);
  for (let j = 0; j < GF.n; j++) out[j] = (1 - TOF_T[j]) * arr[TOF_I[j]] + TOF_T[j] * arr[TOF_I[j] + 1];
  return out;
}
Phys.toFine = toFine;

// ------------------------------------------------------------------ data
Phys.init = function (meta, buf) {
  const get = name => {
    const a = meta.arrays.find(x => x.name === name);
    return a.dtype === "f8" ? new Float64Array(buf, a.offset, a.length)
                            : new Float32Array(buf, a.offset, a.length);
  };
  GC = makeGrid(get("lam"), get("ustar"), meta.T_CMB);
  GF = makeGrid(get("lamf"), get("ustarf"), meta.T_CMB);
  TOF_I = new Int32Array(GF.n); TOF_T = new Float64Array(GF.n);
  for (let j = 0, k = 0; j < GF.n; j++) {
    const x = GF.lnl[j];
    if (x <= GC.lnl[0]) { TOF_I[j] = 0; TOF_T[j] = 0; continue; }
    if (x >= GC.lnl[GC.n - 1]) { TOF_I[j] = GC.n - 2; TOF_T[j] = 1; continue; }
    while (GC.lnl[k + 1] < x) k++;
    TOF_I[j] = k; TOF_T[j] = (x - GC.lnl[k]) / (GC.lnl[k + 1] - GC.lnl[k]);
  }
  MSIL = get("m_sil");
  for (const k of ["pa0.01", "pa0.10", "pe0.01", "pe0.10"]) EPS[k] = get("eps_" + k);
  QT.sil = get("q_sil"); QT.gra = get("q_gra");
  NA = meta.a_grid.n; AGLO = meta.a_grid.lo; AGDL = (meta.a_grid.hi - meta.a_grid.lo) / (NA - 1);
  AG = new Float64Array(NA);
  for (let i = 0; i < NA; i++) AG[i] = Math.pow(10, AGLO + i * AGDL);
  PAG = get("pah_grid");
  TGE = get("tgrid_e"); ESIL = get("e_sil"); EGRA = get("e_gra");
  LTGE = new Float64Array(TGE.length);
  for (let i = 0; i < TGE.length; i++) LTGE[i] = Math.log(TGE[i]);
  DLTGE = (LTGE[LTGE.length - 1] - LTGE[0]) / (LTGE.length - 1);
  MATS = meta.mat;
  Phys.T_CMB = meta.T_CMB;
  Phys.MRN = meta.mrn;
  Phys.GC = GC; Phys.GF = GF;
  Phys.LAM = GC.lam; Phys.LAMF = GF.lam; Phys.A_GRID = AG; Phys.PAH_GRID = PAG;
  Phys.ready = true;
};

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

// refractive index on the coarse grid: [re, im]; graphite gives two (E par, E perp to c)
function refrIndex(mat, a) {
  const n = GC.n;
  if (mat === "sil") {
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let i = 0; i < n; i++) { re[i] = MSIL[2 * i]; im[i] = MSIL[2 * i + 1]; }
    return [[re, im]];
  }
  const w = Math.min(1, Math.max(0, Math.log10(a / 0.01)));
  return ["pa", "pe"].map(o => {
    const e1 = EPS[o + "0.01"], e2 = EPS[o + "0.10"];
    const re = new Float64Array(n), im = new Float64Array(n);
    for (let i = 0; i < n; i++) {
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

// Q_abs, Q_sca, g on the coarse grid for one silicate or graphite grain, by Mie theory
function optics(mat, a) {
  const n = GC.n, ms = refrIndex(mat, a);
  const Qa = new Float64Array(n), Qs = new Float64Array(n), g = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const x = 2 * Math.PI * a / GC.lam[i];
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

// Q_abs and Q_sca from the table at size-grid node i (float32, as shipped)
function tableOptics(mat, i) {
  const n = GC.n, T = QT[mat], Qa = new Float64Array(n), Qs = new Float64Array(n);
  for (let k = 0; k < n; k++) { Qa[k] = T[i * n + k]; Qs[k] = T[(NA + i) * n + k]; }
  return { Qabs: Qa, Qsca: Qs, g: null };
}
Phys.tableOptics = tableOptics;

// ------------------------------------------------------------------ PAHs
// Li & Draine (2001) Eqs. (5)-(11) with the band parameters of Draine & Li (2007) Table 1:
// [lambda_j (um), gamma_j, sigma_int neutral, ionized (1e-20 cm/C), scaled by H/C?].
// See dustlib.PAH_BANDS for the notes.
const PAH_UV = [[0.0722, 0.195, 7.97e7], [0.2175, 0.217, 1.23e7]];
const PAH_BANDS = [
  [1.050, 0.055, 0.0, 2.0e4, 0], [1.260, 0.11, 0.0, 0.078, 0], [1.905, 0.09, 0.0, -146.5, 0],
  [3.300, 0.012, 394.0, 89.4, 1], [5.270, 0.034, 2.5, 20.0, 0], [5.700, 0.035, 4.0, 32.0, 0],
  [6.220, 0.030, 29.4, 235.0, 0], [6.690, 0.070, 7.35, 59.0, 0], [7.417, 0.126, 20.8, 181.0, 0],
  [7.598, 0.044, 18.1, 163.0, 0], [7.850, 0.053, 21.9, 197.0, 0], [8.330, 0.052, 6.94, 48.0, 1],
  [8.610, 0.039, 27.8, 194.0, 1], [10.68, 0.020, 0.3, 0.3, 1], [11.23, 0.012, 18.9, 17.7, 1],
  [11.33, 0.032, 52.0, 49.0, 1], [11.99, 0.045, 24.2, 20.5, 1], [12.62, 0.042, 35.0, 31.0, 1],
  [12.69, 0.013, 1.3, 1.3, 1], [13.48, 0.040, 8.0, 8.0, 1], [14.19, 0.025, 0.45, 0.45, 0],
  [15.90, 0.020, 0.04, 0.04, 0], [16.45, 0.014, 0.5, 0.5, 0], [17.04, 0.065, 2.22, 2.22, 0],
  [17.375, 0.012, 0.11, 0.11, 0], [17.87, 0.016, 0.067, 0.067, 0], [18.92, 0.10, 0.10, 0.17, 0],
  [15.0, 0.8, 50.0, 50.0, 0]];
const PAH_AMIN = 3.55e-4;                      // um; N_C = 20
const HC_OVER_K = 1.4387769, CH_MODES = [886.0, 1161.0, 3030.0];
const M_C = 12.011 * AMU, M_H = 1.008 * AMU;
Phys.PAH_AMIN = PAH_AMIN;

function ncPah(a) { return 460.0 * Math.pow(a / 1e-3, 3); }
function hcRatio(NC) { return NC <= 25 ? 0.5 : NC <= 100 ? 0.5 * Math.sqrt(25.0 / Math.max(NC, 1)) : 0.25; }
function drude(l, lj, gj, sig) {
  const d = l / lj - lj / l;
  return 2 / Math.PI * gj * lj * UM * sig * 1e-20 / (d * d + gj * gj);
}
Phys.ncPah = ncPah; Phys.hcRatio = hcRatio;

// C_abs / N_C (cm^2) on the fine grid for a PAH with NC carbon atoms
function pahCabsPerC(NC, charge) {
  const n = GF.n, out = new Float64Array(n), ion = charge === "ion", hc = hcRatio(NC);
  const M = (NC < 40 ? 0.3 : 0.4) * NC;
  const lc = ion ? 1 / (2.282 / Math.sqrt(M) + 0.889) : 1 / (3.804 / Math.sqrt(M) + 1.052);
  for (let i = 0; i < n; i++) {
    const l = GF.lam[i], x = 1 / l;
    let v = 0;
    if (x > 10 && x < 15) v = drude(l, ...PAH_UV[0]) + (-3.0 + 1.35 * x) * 1e-18;
    else if (x > 7.7 && x <= 10) v = (66.302 - 24.367 * x + 2.950 * x * x - 0.1057 * x * x * x) * 1e-18;
    else if (x > 5.9 && x <= 7.7) {
      const z = x - 5.9;
      v = drude(l, ...PAH_UV[1]) + (1.8687 + 0.1905 * x + 0.4175 * z * z + 0.04370 * z * z * z) * 1e-18;
    } else if (x > 3.3 && x <= 5.9) v = drude(l, ...PAH_UV[1]) + (1.8687 + 0.1905 * x) * 1e-18;
    else if (x <= 3.3) {
      const y = lc / l;
      const cut = Math.atan(1e3 * Math.pow(y - 1, 3) / y) / Math.PI + 0.5;
      v = 34.58e-18 * Math.pow(10, -3.431 / x) * cut;
      for (const b of PAH_BANDS) v += drude(l, b[0], b[1], (ion ? b[3] : b[2]) * (b[4] ? hc : 1));
    }
    if (ion) v += 3.5e-19 * Math.pow(10, -1.45 / x) * Math.exp(-0.1 * x * x);
    out[i] = v;
  }
  return out;
}
Phys.pahCabsPerC = pahCabsPerC;

function xiGra(a) { return a <= 0.005 ? 0.01 : 0.01 + 0.99 * (1 - Math.pow(0.005 / a, 3)); }

// graphite C_abs, C_sca, g on the fine grid: Mie on the coarse grid at max(a, 0.001 um),
// scaled to a in the small-particle limit, then interpolated
function graphiteOnFine(a) {
  const ae = Math.max(a, 1e-3), o = optics("gra", ae), geo = Math.PI * Math.pow(ae * UM, 2), s = a / ae;
  const ca = new Float64Array(GC.n), cs = new Float64Array(GC.n);
  for (let i = 0; i < GC.n; i++) { ca[i] = geo * o.Qabs[i]; cs[i] = geo * o.Qsca[i]; }
  const fa = toFine(ca), fs = toFine(cs), fg = toFine(o.g), s3 = s * s * s, s6 = s3 * s3;
  for (let j = 0; j < GF.n; j++) { fa[j] *= s3; fs[j] *= s6; }
  return [fa, fs, fg];
}

// Q_abs, Q_sca, g on the fine grid for a PAH-graphite grain (floored at zero; see dustlib)
function pahOptics(a, charge) {
  const NC = ncPah(a), [ca, cs, g] = graphiteOnFine(a), cper = pahCabsPerC(NC, charge);
  const xi = xiGra(a), geo = Math.PI * Math.pow(a * UM, 2);
  const Qa = new Float64Array(GF.n), Qs = new Float64Array(GF.n);
  for (let j = 0; j < GF.n; j++) {
    Qa[j] = Math.max((1 - xi) * NC * cper[j] + xi * ca[j], 0) / geo;
    Qs[j] = cs[j] / geo;
  }
  return { Qabs: Qa, Qsca: Qs, g };
}
Phys.pahOptics = pahOptics;

// ------------------------------------------------------------------ enthalpy
// energy of the lowest vibrational mode (erg); see dustlib.e_first_mode
function eFirstMode(mat, N) {
  if (mat === "sil") {
    let best = Infinity;
    for (const [n, th, Nm] of [[2, 500.0, 2 * (N - 2)], [3, 1500.0, N - 2]]) {
      const beta = (Math.pow(Nm, 1 - n / 3) - 1) / (2 * Nm - 1);
      best = Math.min(best, KB * th * Math.pow((1 - beta) * 0.5 / Nm + beta, 1 / n));
    }
    return best;
  }
  const NC = N, Nm = N - 2;
  const beta = NC <= 54 ? 0 : NC <= 102 ? (NC - 54) / 52 / (2 * Nm - 1)
             : ((NC - 2) / 52 * Math.pow(102 / NC, 2 / 3) - 1) / (2 * Nm - 1);
  return KB * 863.0 * Math.sqrt((1 - beta) * 0.5 / Nm + beta);
}
Phys.eFirstMode = eFirstMode;
function tabE(lE, T) {
  const n = LTGE.length, u = (Math.log(T) - LTGE[0]) / DLTGE;
  if (u <= 0) return Math.exp(lE[0]);
  if (u >= n - 1) return Math.exp(lE[n - 1]);
  let k = Math.floor(u);
  const lt = Math.log(T);
  while (k > 0 && LTGE[k] > lt) k--;
  while (k < n - 2 && LTGE[k + 1] < lt) k++;
  const t = (lt - LTGE[k]) / (LTGE[k + 1] - LTGE[k]);
  return Math.exp(lE[k] + t * (lE[k + 1] - lE[k]));
}
function Enthalpy(mat, a) {
  const n = TGE.length;
  this.lE = new Float64Array(n);
  if (mat === "pah") {
    // N_H Einstein oscillators per C-H mode plus the Debye C-C modes (DL01)
    this.NC = ncPah(a); this.NH = hcRatio(this.NC) * this.NC; this.N = this.NC + this.NH;
    this.E1 = eFirstMode("pah", this.NC);
    this.mass = this.NC * M_C + this.NH * M_H;
    for (let i = 0; i < n; i++) {
      const T = TGE[i];
      let ech = 0;
      for (const w of CH_MODES) ech += KB * HC_OVER_K * w / Math.expm1(Math.min(HC_OVER_K * w / T, 700));
      this.lE[i] = Math.log(this.NH * ech + (this.NC - 2) * EGRA[i]);
    }
  } else {
    this.N = nAtoms(mat, a);
    this.E1 = eFirstMode(mat, this.N);
    const e = mat === "sil" ? ESIL : EGRA;
    for (let i = 0; i < n; i++) this.lE[i] = Math.log((this.N - 2) * e[i]);
  }
}
Enthalpy.prototype.E = function (T) { return tabE(this.lE, T); };
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

// mat "sil", "gra" (coarse grid) or "pah" (fine grid, with charge "neu" or "ion")
function Grain(mat, a, U, opt, charge) {
  this.mat = mat; this.a = a; this.U = U; this.charge = charge || "neu";
  const G = this.G = mat === "pah" ? GF : GC, n = G.n;
  const o = opt || (mat === "pah" ? pahOptics(a, this.charge) : optics(mat, a));
  this.Qabs = o.Qabs; this.Qsca = o.Qsca; this.g = o.g;
  const geo = Math.PI * Math.pow(a * UM, 2);
  this.geo = geo;
  this.Cabs = new Float64Array(n);
  this.u = new Float64Array(n); this.us = new Float64Array(n);
  this.pabsLam = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    this.Cabs[i] = geo * this.Qabs[i];
    this.us[i] = U * G.ustar[i];
    this.u[i] = this.us[i] + G.ucmb[i];
    this.pabsLam[i] = C * this.Cabs[i] * this.u[i];
  }
  this.H = new Enthalpy(mat, a);
  this.mass = mat === "pah" ? this.H.mass : grainMass(mat, a);
  this.Pabs = integ(G, this.pabsLam);
  this.Teq = this.equilibriumT();
}

// emitted power at T; with Emax, only photons less energetic than Emax
Grain.prototype.L = function (T, Emax) {
  const G = this.G;
  let s = 0;
  const lim = Emax === undefined ? Infinity : Emax;
  for (let i = 0; i < G.n; i++) {
    if (G.eph[i] >= lim) continue;
    const x = G.xt[i] / T;
    if (x >= XMAX) continue;
    s += G.wl[i] * this.Cabs[i] * G.b0[i] / Math.expm1(x);
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

// the spectrum the grain would emit at its equilibrium temperature
Grain.prototype.eqEmission = function () {
  const G = this.G, em = new Float64Array(G.n);
  for (let i = 0; i < G.n; i++) em[i] = 4 * Math.PI * this.Cabs[i] * planck(G, i, this.Teq);
  return em;
};

// starlight photon statistics: absorption rate, power, and int (h nu) dP
Grain.prototype.photonStats = function () {
  if (this._ps) return this._ps;
  const G = this.G;
  let R = 0, S1 = 0, S2 = 0;
  for (let i = 0; i < G.n; i++) {
    const p = C * this.Cabs[i] * this.us[i];
    R += G.wl[i] * p / G.eph[i]; S1 += G.wl[i] * p; S2 += G.wl[i] * p * G.eph[i];
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
  const G = this.G, n = G.n;
  const I0 = new Float64Array(n), I1 = new Float64Array(n);
  let pg0 = 0, pg1 = 0, pE = 0;
  for (let k = 0; k < n; k++) {
    const i = n - 1 - k, E = G.eph[i];
    const f = this.Cabs[i] * this.u[i] * G.lam[i] / E;
    const g0 = f * E, g1 = f * E * E;
    if (k > 0) {
      const dl = Math.log(E / pE);
      I0[k] = I0[k - 1] + 0.5 * (g0 + pg0) * dl; I1[k] = I1[k - 1] + 0.5 * (g1 + pg1) * dl;
    }
    pg0 = g0; pg1 = g1; pE = E;
  }
  return (this._cum = [I0, I1]);
};

Grain.prototype.stochastic = function (Mfirst) {
  const M = Mfirst || M0, Hh = this.H, G = this.G;
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
  // never below the lowest vibrational mode, nor two quanta of the longest grid wavelength;
  // energies under E1 belong to the ground state, which the bottom bin stands for (see dustlib)
  const Erad = 2 * G.eph[G.n - 1];
  const Elo = Math.max(Hh.E(Phys.T_CMB), Erad, Hh.E1, Math.min(Eeq - 10 * sE, Hh.E(Tcool)));
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
// Upward rates are banded (no photon carries more than the top grid energy) apart from the
// row into the top bin, which is handled as a running sum.
Grain.prototype.solve = function (T) {
  const Hh = this.H, M = T.length - 1, G = this.G, NG = G.n, ETOP = G.etop;
  const E = new Float64Array(M + 1);
  for (let j = 0; j <= M; j++) E[j] = Hh.E(T[j]);
  const Emin = new Float64Array(M + 1), Emax = new Float64Array(M + 1), dE = new Float64Array(M + 1);
  for (let j = 0; j <= M; j++) {
    Emin[j] = j === 0 ? Math.max(0, E[0] - 0.5 * (E[1] - E[0])) : 0.5 * (E[j] + E[j - 1]);
    Emax[j] = j === M ? E[M] + 0.5 * (E[M] - E[M - 1]) : 0.5 * (E[j] + E[j + 1]);
    dE[j] = Emax[j] - Emin[j];
  }
  const [I0, I1] = this.cumulants();
  const I0top = I0[NG - 1], I1top = I1[NG - 1];
  // ascending photon energies; the cumulants are interpolated linearly in E
  const EgA = G.egA;
  const Eg0 = EgA[0], EgTop = EgA[NG - 1];
  const at = (I, k, w) => I[k] + (w - EgA[k]) / (EgA[k + 1] - EgA[k]) * (I[k + 1] - I[k]);
  const locate = w => {                    // k with EgA[k] <= w < EgA[k+1]
    let lo = 0, hi = NG - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (EgA[m] > w) hi = m; else lo = m; }
    return lo;
  };
  const cum = (I, w) => (w <= Eg0 ? 0 : w >= EgTop ? I[NG - 1] : at(I, locate(w), w));
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
      const Gk = ((c1W2 - c11[k]) - W1 * (c0W2 - c01[k]) + Math.min(dE[l], dE[u]) * (c0W3 - c0W2)
                 + W4 * (c02[k + 1] - c0W3) - (c12[k + 1] - c1W3)) / (dE[u] * dE[l]);
      band[off[l] + k] = C * dE[u] / (E[u] - E[l]) * Gk;
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
  for (let u = 1; u <= M; u++) Dn[u] = Math.max(down[u], 1e-300) / (E[u] - E[u - 1]);
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
  for (let j = 0; j <= M; j++) { P[j] = X[j] / sum; if (j > 0) pem += P[j] * down[j]; }
  // the bottom bin has no downward transition: the ground state (or empty); it radiates nothing
  return { T, P, E, down, Pemit: pem, ground: E[0] <= 1.0001 * Hh.E1 ? P[0] : 0 };
};

// emitted power per um per grain; also the temperatures and probabilities behind it
Grain.prototype.emission = function () {
  const G = this.G, NG = G.n;
  const sT = this.sigmaT()[0];
  if (sT < SIGMA_EQ) return { em: this.eqEmission(), T: [this.Teq], P: [1], stochastic: false, sigmaT: sT };
  const em = new Float64Array(NG);
  const r = this.stochastic();
  for (let j = 1; j < r.T.length; j++) {                  // the bottom bin radiates nothing
    const Pj = r.P[j];
    if (Pj < 1e-30) continue;
    const Tj = r.T[j], Ej = r.E[j];
    for (let i = 0; i < NG; i++) {
      if (G.eph[i] >= Ej) continue;
      const x = G.xt[i] / Tj;
      if (x >= XMAX) continue;
      em[i] += Pj * G.b0[i] / Math.expm1(x);
    }
  }
  for (let i = 0; i < NG; i++) em[i] *= 4 * Math.PI * this.Cabs[i];
  return { em, T: r.T, P: r.P, E: r.E, down: r.down, ground: r.ground, stochastic: true, sigmaT: sT, Pemit: r.Pemit };
};
Phys.Grain = Grain;

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

// temperatures below which a fraction q of the power is radiated (levels j >= 1, weights P_j L_j)
function powerPercentiles(r, qs) {
  if (!r.stochastic) return qs.map(() => r.T[0]);
  const w = Array.from(r.P, (p, j) => (j === 0 ? 0 : p * r.down[j]));
  const tot = w.reduce((a, b) => a + b, 0);
  return percentiles(r.T, w.map(v => v / tot), qs);
}
Phys.powerPercentiles = powerPercentiles;

// ------------------------------------------------------------------ size distributions
const SIZE_BINS = [[0, 0.01], [0.01, 0.1], [0.1, Infinity]];
Phys.SIZE_BINS = SIZE_BINS;
const binOf = a => SIZE_BINS.findIndex(r => a >= r[0] && a < r[1]);

// power laws (silicate, graphite): dn/da = A a^-q on [amin, amax], normalized to abund ppm
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

// dn/da = A a^-q (a in um, per um per H) for the normalization of powerlawWeights
function powerlawA(mat, q, amin, amax, abund) {
  const { a, w } = powerlawWeights(mat, q, amin, amax, abund);
  const dl = Math.log(a[1] / a[0]);
  return w[1] / (Math.pow(a[1], 1 - q) * dl);
}
Phys.powerlawA = powerlawA;

// PAHs: log-normal dn/da = (n0/a) exp[-ln^2(a/a0) / 2 sigma^2] on [PAH_AMIN, pahAmax]
function pahAmax(a0, sigma) { return Math.min(0.05, Math.max(2 * a0, a0 * Math.exp(4 * sigma))); }
function lognormalWeights(a0, sigma, abund) {
  const lo = Math.log(PAH_AMIN), hi = Math.log(pahAmax(a0, sigma));
  const n = Math.max(Math.ceil((hi - lo) / Math.LN10 * 120) + 1, 2);
  const a = new Float64Array(n), w = new Float64Array(n), dl = (hi - lo) / (n - 1);
  let norm = 0;
  for (let k = 0; k < n; k++) {
    const lx = lo + k * dl;
    a[k] = Math.exp(lx);
    const tw = (k === 0 || k === n - 1) ? 0.5 * dl : dl;
    const z = (lx - Math.log(a0)) / sigma;
    w[k] = Math.exp(-0.5 * z * z) * tw;
    norm += w[k] * ncPah(a[k]);
  }
  const A = abund * 1e-6 / norm;
  for (let k = 0; k < n; k++) w[k] *= A;
  return { a, w, n0: A };
}
Phys.lognormalWeights = lognormalWeights; Phys.pahAmax = pahAmax;
Phys.pahMass = a => { const NC = ncPah(a); return NC * M_C + hcRatio(NC) * NC * M_H; };

// fractional position of radius a among nodes (numpy searchsorted, side left)
function nodeIndex(nodes, a, clampEnds) {
  const n = nodes.length;
  let la = Math.log(a);
  if (clampEnds) la = Math.min(Math.max(la, Math.log(nodes[0])), Math.log(nodes[n - 1]));
  let lo = 0, hi = n;
  while (lo < hi) { const m = (lo + hi) >> 1; if (Math.log(nodes[m]) < la) lo = m + 1; else hi = m; }
  const i = Math.max(0, Math.min(n - 2, lo - 1));
  const l0 = Math.log(nodes[i]), l1 = Math.log(nodes[i + 1]);
  return [i, (la - l0) / (l1 - l0)];
}
const aIndex = a => nodeIndex(AG, a, false);

// Silicate or graphite: extinction and scattering per H (cm^2) by material and size bin,
// on the fine grid, and the mass per H
function distOptics(mat, p) {
  const n = GC.n, { a, w } = powerlawWeights(mat, p.q, p.amin, p.amax, p.abund);
  const ext = new Float64Array(n), sca = new Float64Array(n), T = QT[mat];
  const extB = SIZE_BINS.map(() => new Float64Array(n)), scaB = SIZE_BINS.map(() => new Float64Array(n));
  let mass = 0;
  for (let k = 0; k < a.length; k++) {
    const [i, t] = aIndex(a[k]);
    const geo = w[k] * Math.PI * Math.pow(a[k] * UM, 2);
    const eb = extB[binOf(a[k])], sb = scaB[binOf(a[k])];
    for (let j = 0; j < n; j++) {
      const qa = (1 - t) * T[i * n + j] + t * T[(i + 1) * n + j];
      const qs = (1 - t) * T[(NA + i) * n + j] + t * T[(NA + i + 1) * n + j];
      ext[j] += geo * (qa + qs); sca[j] += geo * qs;
      eb[j] += geo * (qa + qs); sb[j] += geo * qs;
    }
    mass += w[k] * grainMass(mat, a[k]);
  }
  return { ext: toFine(ext), sca: toFine(sca), mass, extB: extB.map(toFine), scaB: scaB.map(toFine) };
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

// Emission of one grain at a node with a summary of its temperatures: T_eq, the 5/50/95th
// percentiles of P(T), and the temperature after absorbing one photon of the mean energy.
// mat "sil"/"gra" use the size grid and table optics; "pah" uses the PAH nodes and a charge.
function nodeEmission(mat, i, U, charge) {
  const g = mat === "pah" ? new Grain("pah", PAG[i], U, pahNodeOptics(i, charge), charge)
                          : new Grain(mat, AG[i], U, tableOptics(mat, i));
  const r = g.emission();
  const ps = g.photonStats();
  const out = { em: r.em, pabs: g.Pabs, stochastic: r.stochastic, Teq: g.Teq,
                Tjump: g.H.T(g.H.E(g.Teq) + ps[1] / ps[0]), ground: r.ground || 0 };
  [out.w05, out.w95] = powerPercentiles(r, [0.05, 0.95]);
  return out;
}
Phys.nodeEmission = nodeEmission;

// Silicate or graphite: emission per H on the fine grid, from the per-node emission
function distEmission(mat, p, nodes) {
  const n = GC.n, { a, w } = powerlawWeights(mat, p.q, p.amin, p.amax, p.abund);
  const em = new Float64Array(n);
  const bins = SIZE_BINS.map(() => new Float64Array(n));
  let pabs = 0;
  for (let k = 0; k < a.length; k++) {
    const [i, t] = aIndex(a[k]);
    const n0 = nodes[i], n1 = nodes[i + 1];
    const s0 = (1 - t) / Math.pow(AG[i], 3), s1 = t / Math.pow(AG[i + 1], 3), a3 = Math.pow(a[k], 3);
    const bn = bins[binOf(a[k])];
    for (let j = 0; j < n; j++) {
      const v = w[k] * a3 * (s0 * n0.em[j] + s1 * n1.em[j]);
      em[j] += v; bn[j] += v;
    }
    pabs += w[k] * a3 * (s0 * n0.pabs + s1 * n1.pabs);
  }
  for (let j = 0; j < n; j++) {
    const f = GC.lam[j] / (4 * Math.PI);
    em[j] *= f;
    for (const b of bins) b[j] *= f;
  }
  return { em: toFine(em), bins: bins.map(toFine), pabs };
}
Phys.distEmission = distEmission;

// PAH optics at the PAH nodes (independent of U), cached
const PNO = {};
function pahNodeOptics(i, charge) {
  const k = charge + i;
  return PNO[k] || (PNO[k] = pahOptics(PAG[i], charge));
}
Phys.pahNodeOptics = pahNodeOptics;

// PAH extinction and scattering per H on the fine grid, by size bin; p = {a0, sigma, abund, fion}
function distOpticsPAH(p) {
  const n = GF.n, { a, w } = lognormalWeights(p.a0, p.sigma, p.abund), f = p.fion;
  const ext = new Float64Array(n), sca = new Float64Array(n);
  const extB = SIZE_BINS.map(() => new Float64Array(n)), scaB = SIZE_BINS.map(() => new Float64Array(n));
  let mass = 0;
  for (let k = 0; k < a.length; k++) {
    const [i, t] = nodeIndex(PAG, a[k], true);
    const on = [pahNodeOptics(i, "neu"), pahNodeOptics(i, "ion"), pahNodeOptics(i + 1, "neu"), pahNodeOptics(i + 1, "ion")];
    const g0 = Math.PI * Math.pow(PAG[i] * UM, 2), g1 = Math.PI * Math.pow(PAG[i + 1] * UM, 2);
    const s0 = (1 - t) / Math.pow(PAG[i], 3), s1 = t / Math.pow(PAG[i + 1], 3), wa3 = w[k] * Math.pow(a[k], 3);
    const eb = extB[binOf(a[k])], sb = scaB[binOf(a[k])];
    for (let j = 0; j < n; j++) {
      const ca = s0 * g0 * ((1 - f) * on[0].Qabs[j] + f * on[1].Qabs[j]) + s1 * g1 * ((1 - f) * on[2].Qabs[j] + f * on[3].Qabs[j]);
      const cs = s0 * g0 * on[0].Qsca[j] + s1 * g1 * on[2].Qsca[j];
      const e = wa3 * (ca + cs), s = wa3 * cs;
      ext[j] += e; sca[j] += s; eb[j] += e; sb[j] += s;
    }
    mass += w[k] * Phys.pahMass(a[k]);
  }
  return { ext, sca, mass, extB, scaB };
}
Phys.distOpticsPAH = distOpticsPAH;

// PAH emission per H on the fine grid; nodes = {neu: [...], ion: [...]} indexed by PAH node
function distEmissionPAH(p, nodes) {
  const n = GF.n, { a, w } = lognormalWeights(p.a0, p.sigma, p.abund), f = p.fion;
  const em = new Float64Array(n), bins = SIZE_BINS.map(() => new Float64Array(n));
  let pabs = 0;
  for (let k = 0; k < a.length; k++) {
    const [i, t] = nodeIndex(PAG, a[k], true);
    const s0 = (1 - t) / Math.pow(PAG[i], 3), s1 = t / Math.pow(PAG[i + 1], 3), wa3 = w[k] * Math.pow(a[k], 3);
    const n0 = nodes.neu[i], i0 = nodes.ion[i], n1 = nodes.neu[i + 1], i1 = nodes.ion[i + 1];
    const bn = bins[binOf(a[k])];
    for (let j = 0; j < n; j++) {
      const v = wa3 * (s0 * ((1 - f) * n0.em[j] + f * i0.em[j]) + s1 * ((1 - f) * n1.em[j] + f * i1.em[j]));
      em[j] += v; bn[j] += v;
    }
    pabs += wa3 * (s0 * ((1 - f) * n0.pabs + f * i0.pabs) + s1 * ((1 - f) * n1.pabs + f * i1.pabs));
  }
  for (let j = 0; j < n; j++) {
    const fl = GF.lam[j] / (4 * Math.PI);
    em[j] *= fl;
    for (const b of bins) b[j] *= fl;
  }
  return { em, bins, pabs };
}
Phys.distEmissionPAH = distEmissionPAH;

// PAH nodes needed for a log-normal: the bracket of [PAH_AMIN, pahAmax]
function pahNodeRange(a0, sigma) {
  const [i1] = nodeIndex(PAG, pahAmax(a0, sigma), true);
  return [0, Math.min(PAG.length - 1, i1 + 1)];
}
Phys.pahNodeRange = pahNodeRange;

// ------------------------------------------------------------------ reference extinction curve
// Cardelli, Clayton & Mathis (1989), Eqs. (1)-(5): A_lambda / A_V = a(x) + b(x)/R_V for
// 0.3 <= x <= 10 um^-1 (NaN outside); per H scaled with N_H / E(B-V) = 5.8e21 cm^-2 mag^-1
// (Bohlin, Savage & Drake 1978).  See dustlib.ccm89.
Phys.NH_EBV = 5.8e21;
function ccm89(x, Rv) {
  let a, b;
  if (x >= 0.3 && x < 1.1) { a = 0.574 * Math.pow(x, 1.61); b = -0.527 * Math.pow(x, 1.61); }
  else if (x >= 1.1 && x < 3.3) {
    const y = x - 1.82;
    const pa = [1.0, 0.17699, -0.50447, -0.02427, 0.72085, 0.01979, -0.77530, 0.32999];
    const pb = [0.0, 1.41338, 2.28305, 1.07233, -5.38434, -0.62251, 5.30260, -2.09002];
    a = 0; b = 0;
    for (let k = 7; k >= 0; k--) { a = a * y + pa[k]; b = b * y + pb[k]; }
  } else if (x >= 3.3 && x <= 8) {
    const fa = x >= 5.9 ? -0.04473 * Math.pow(x - 5.9, 2) - 0.009779 * Math.pow(x - 5.9, 3) : 0;
    const fb = x >= 5.9 ? 0.2130 * Math.pow(x - 5.9, 2) + 0.1207 * Math.pow(x - 5.9, 3) : 0;
    a = 1.752 - 0.316 * x - 0.104 / (Math.pow(x - 4.67, 2) + 0.341) + fa;
    b = -3.090 + 1.825 * x + 1.206 / (Math.pow(x - 4.62, 2) + 0.263) + fb;
  } else if (x > 8 && x <= 10) {
    const z = x - 8;
    a = -1.073 - 0.628 * z + 0.137 * z * z - 0.070 * z * z * z;
    b = 13.670 + 4.257 * z - 0.420 * z * z + 0.374 * z * z * z;
  } else return NaN;
  return a + b / Rv;
}
Phys.ccm89 = ccm89;

if (typeof module !== "undefined" && module.exports) module.exports = Phys;
else root.Phys = Phys;
})(typeof window !== "undefined" ? window : this);
