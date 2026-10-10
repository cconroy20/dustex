/* The page: state, controls and drawing.  All of the physics is in physics.js. */
(function () {
"use strict";

window.addEventListener("error", ev => {
  const e = document.getElementById("err");
  e.hidden = false;
  e.textContent = String(ev.message) + (ev.filename ? "\n" + ev.filename + ":" + ev.lineno : "");
});

const P = window.Phys;
const FONT = 'Charter, "PT Serif", Palatino, Georgia, "DejaVu Serif", serif';
// colorblind-safe pairs (dataviz validator: protan/deutan dE >= 9.7 for every pair used together)
const C = { ink: "#1a1a1a", lab: "#3b3733", muted: "#6a6560", faint: "#a39c92", frame: "#cfc8bd",
            abs: "#004C8C", sca: "#00897B", emi: "#C2571A",
            sil: "#B5533C", gra: "#2F5D8C", size: ["#9D7BC9", "#6A3D9A", "#3B1F5E"] };
const MATNAME = { sil: "Silicate", gra: "Graphite" };
const $ = id => document.getElementById(id);
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const L10 = v => (v > 0 ? Math.log10(v) : -Infinity);
const LN10 = Math.LN10;
const MH = 1.6735575e-24;          // g, mass of an H atom
const KEY_UM = { V: 0.55, B: 0.44 };

// slider ranges (log10 where noted)
const R = { U: [-2, 4], a: [Math.log10(0.003), Math.log10(3)], q: [2.5, 4.5],
            amin: [Math.log10(0.003), Math.log10(0.1)], amax: [Math.log10(0.03), Math.log10(3)],
            ab: { sil: [0, 60], gra: [0, 400] } };
const T_NODES = [11, 88];          // size-grid nodes spanning the radius slider (0.003-3 um)

let S = null, DEF = null;

function mrnDist() {
  const d = {};
  for (const m of ["sil", "gra"])
    d[m] = { on: true, q: 3.5, amin: Math.log10(0.005), amax: Math.log10(0.25), abund: +P.MRN[m].toFixed(1) };
  return d;
}
function defaults() {
  return { tab: "grain", U: 0, mat: "sil", a: -1, dist: mrnDist(),
           elog: false, plin: false, sd: "mass", xn: "H", xs: "mat", sp: "mat" };
}
const clone = o => JSON.parse(JSON.stringify(o));

// ------------------------------------------------------------------ text
function expStr(e) { return String(e).replace("-", "−"); }
function sciParts(v, d) {
  if (v === 0) return [{ t: "0" }];
  if (!isFinite(v)) return [{ t: "∞" }];
  let e = Math.floor(Math.log10(Math.abs(v)));
  let ms = (v / Math.pow(10, e)).toFixed(d - 1);
  if (Math.abs(parseFloat(ms)) >= 10) { e += 1; ms = (v / Math.pow(10, e)).toFixed(d - 1); }
  if (e === 0) return [{ t: ms }];
  return [{ t: ms + " × 10" }, { t: expStr(e), sup: true }];
}
function html(parts) {
  return parts.map(p => {
    let s = p.t.replace(/</g, "&lt;");
    if (p.it) s = `<i>${s}</i>`;
    return p.sup ? `<sup>${s}</sup>` : p.sub ? `<sub>${s}</sub>` : s;
  }).join("");
}
function numHTML(v, d = 3) {
  const a = Math.abs(v);
  if (a >= 0.01 && a < 1e4) {
    let s = v.toFixed(Math.max(0, d - 1 - Math.floor(Math.log10(a))));
    if (s.includes(".")) s = s.replace(/0+$/, "").replace(/\.$/, "");
    return s;
  }
  return html(sciParts(v, d));
}
function aHTML(a) {                 // a radius in um
  return (a < 0.01 ? a.toFixed(4) : a < 0.1 ? a.toFixed(3) : a < 1 ? a.toFixed(2) : a.toFixed(2)) + " µm";
}
function tHTML(T) { return (T < 10 ? T.toFixed(1) : T < 1000 ? T.toFixed(T < 100 ? 1 : 0) : T.toFixed(0)) + " K"; }
function durHTML(t) {
  if (!isFinite(t)) return "∞";
  if (t < 1) return numHTML(t, 2) + " s";
  if (t < 120) return t.toFixed(t < 10 ? 1 : 0) + " s";
  if (t < 7200) return (t / 60).toFixed(t < 600 ? 1 : 0) + " min";
  if (t < 2 * 86400) return (t / 3600).toFixed(t < 36000 ? 1 : 0) + " hr";
  if (t < 3.156e7) return (t / 86400).toFixed(t < 864000 ? 1 : 0) + " days";
  return numHTML(t / 3.156e7, 2) + " yr";
}
function uHTML(u) { return numHTML(Math.pow(10, u), 2); }

// a tiny markup for axis labels: *...* italic, ^{...} superscript, _{...} subscript
function mk(s) {
  const parts = [];
  let it = false, mode = "", buf = "";
  const push = () => { if (buf) parts.push({ t: buf, it, sup: mode === "sup", sub: mode === "sub" }); buf = ""; };
  for (let i = 0; i < s.length; i++) {
    const ch = s[i];
    if ((ch === "^" || ch === "_") && s[i + 1] === "{") { push(); mode = ch === "^" ? "sup" : "sub"; i++; continue; }
    if (ch === "}" && mode) { push(); mode = ""; continue; }
    if (ch === "*") { push(); it = !it; continue; }
    buf += ch;
  }
  push();
  return parts;
}

// draw rich text; returns its width
function rich(ctx, parts, x, y, align, size, color) {
  const small = Math.round(size * 0.72);
  const fnt = p => `${p.it ? "italic " : ""}${p.sup || p.sub ? small : size}px ${FONT}`;
  const ws = parts.map(p => { ctx.font = fnt(p); return ctx.measureText(p.t).width; });
  const w = ws.reduce((s, v) => s + v, 0);
  let cx = align === "center" ? x - w / 2 : align === "right" ? x - w : x;
  ctx.textAlign = "left"; ctx.textBaseline = "alphabetic"; ctx.fillStyle = color;
  parts.forEach((p, i) => {
    ctx.font = fnt(p);
    ctx.fillText(p.t, cx, y + (p.sup ? -size * 0.42 : p.sub ? size * 0.24 : 0));
    cx += ws[i];
  });
  return w;
}

// ------------------------------------------------------------------ canvas
function fit(cv) {
  const dpr = window.devicePixelRatio || 1;
  const w = cv.clientWidth, h = cv.clientHeight;
  const W = Math.round(w * dpr), H = Math.round(h * dpr);
  if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  return { ctx, w, h };
}
// axes in plot units (log10 on log axes)
function axes(box, x0, x1, y0, y1) {
  const px = u => box.l + (u - x0) / (x1 - x0) * (box.r - box.l);
  const py = u => clamp(box.b - (u - y0) / (y1 - y0) * (box.b - box.t), box.t - 4000, box.b + 4000);
  const ux = p => x0 + (p - box.l) / (box.r - box.l) * (x1 - x0);
  const uy = p => y0 + (box.b - p) / (box.b - box.t) * (y1 - y0);
  return { box, x0, x1, y0, y1, px, py, ux, uy };
}
function frame(ctx, A) {
  const b = A.box;
  ctx.strokeStyle = C.frame; ctx.lineWidth = 1; ctx.setLineDash([]);
  ctx.strokeRect(b.l + 0.5, b.t + 0.5, b.r - b.l, b.b - b.t);
}
function decParts(e, plain) {
  if (plain && e >= -3 && e <= 3) return [{ t: e < 0 ? Math.pow(10, e).toFixed(-e) : String(Math.pow(10, e)) }];
  if (e === 0) return [{ t: "1" }];
  if (e === 1) return [{ t: "10" }];
  return [{ t: "10" }, { t: expStr(e), sup: true }];
}
function subLabel(m, e) {
  const v = m * Math.pow(10, e);
  return [{ t: e < 0 ? v.toFixed(-e) : String(Math.round(v)) }];
}
// log ticks on axis "x" or "y"; labels every decade, with 2 and 5 (or all digits) on short spans
function logTicks(ctx, A, ax, o = {}) {
  const b = A.box, lo = ax === "x" ? A.x0 : A.y0, hi = ax === "x" ? A.x1 : A.y1, span = hi - lo;
  const every = o.every || 1, labels = o.labels !== false, plain = !!o.plain;
  const subs = plain ? (span < 0.6 ? [2, 3, 4, 5, 6, 7, 8, 9] : span < 1.6 ? [2, 5] : []) : [];
  ctx.strokeStyle = C.faint; ctx.lineWidth = 1;
  for (let e = Math.floor(lo); e <= Math.ceil(hi); e++) {
    for (let m = 1; m < 10; m++) {
      const u = e + Math.log10(m);
      if (u < lo - 1e-9 || u > hi + 1e-9) continue;
      const len = m === 1 ? 6 : 3;
      if (ax === "x") {
        const x = Math.round(A.px(u)) + 0.5;
        ctx.beginPath(); ctx.moveTo(x, b.b); ctx.lineTo(x, b.b - len); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(x, b.t); ctx.lineTo(x, b.t + len); ctx.stroke();
        if (!labels) continue;
        if (m === 1 && ((e % every) + every) % every === 0) rich(ctx, decParts(e, plain), x, b.b + 19, "center", 14, C.lab);
        else if (subs.includes(m)) rich(ctx, subLabel(m, e), x, b.b + 19, "center", 13, C.muted);
      } else {
        const y = Math.round(A.py(u)) + 0.5;
        ctx.beginPath(); ctx.moveTo(b.l, y); ctx.lineTo(b.l + len, y); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(b.r, y); ctx.lineTo(b.r - len, y); ctx.stroke();
        if (!labels) continue;
        if (m === 1 && ((e % every) + every) % every === 0) rich(ctx, decParts(e, plain), b.l - 7, y + 5, "right", 14, C.lab);
        else if (subs.includes(m)) rich(ctx, subLabel(m, e), b.l - 7, y + 5, "right", 13, C.muted);
      }
    }
  }
}
function niceStep(span, n) {
  const raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
}
function fmtTick(v, step) {
  const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return v.toFixed(d).replace("-", "−");
}
function linTicks(ctx, A, ax, step, labels = true) {
  const b = A.box, lo = ax === "x" ? A.x0 : A.y0, hi = ax === "x" ? A.x1 : A.y1;
  ctx.strokeStyle = C.faint; ctx.lineWidth = 1;
  for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) {
    if (ax === "x") {
      const x = Math.round(A.px(v)) + 0.5;
      ctx.beginPath(); ctx.moveTo(x, b.b); ctx.lineTo(x, b.b - 6); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x, b.t); ctx.lineTo(x, b.t + 6); ctx.stroke();
      if (labels) rich(ctx, [{ t: fmtTick(v, step) }], x, b.b + 19, "center", 14, C.lab);
    } else {
      const y = Math.round(A.py(v)) + 0.5;
      ctx.beginPath(); ctx.moveTo(b.l, y); ctx.lineTo(b.l + 6, y); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(b.r, y); ctx.lineTo(b.r - 6, y); ctx.stroke();
      if (labels) rich(ctx, [{ t: fmtTick(v, step) }], b.l - 7, y + 5, "right", 14, C.lab);
    }
  }
}
function xlabel(ctx, A, parts) { rich(ctx, parts, (A.box.l + A.box.r) / 2, A.box.b + 41, "center", 15, C.muted); }
function ylabel(ctx, A, parts, dx = 52) {
  ctx.save();
  ctx.translate(A.box.l - dx, (A.box.t + A.box.b) / 2);
  ctx.rotate(-Math.PI / 2);
  rich(ctx, parts, 0, 0, "center", 15, C.muted);
  ctx.restore();
}
function clip(ctx, A) {
  const b = A.box;
  ctx.save(); ctx.beginPath(); ctx.rect(b.l, b.t, b.r - b.l, b.b - b.t); ctx.clip();
}
// polyline through (u, v) in plot units, broken where v is not finite
function curve(ctx, A, us, vs, color, width, dash, alpha) {
  clip(ctx, A);
  ctx.globalAlpha = alpha === undefined ? 1 : alpha;
  ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []);
  ctx.lineJoin = "round"; ctx.lineCap = dash && dash.length ? "butt" : "round";
  ctx.beginPath();
  let on = false;
  for (let i = 0; i < us.length; i++) {
    const v = vs[i];
    if (!isFinite(v) || !isFinite(us[i])) { on = false; continue; }
    const x = A.px(us[i]), y = A.py(v);
    if (on) ctx.lineTo(x, y); else ctx.moveTo(x, y);
    on = true;
  }
  ctx.stroke();
  ctx.restore();
}
// filled region between v (top) and a baseline in plot units
function fillTo(ctx, A, us, vs, base, color, alpha) {
  clip(ctx, A);
  ctx.globalAlpha = alpha; ctx.fillStyle = color;
  ctx.beginPath();
  let started = false, last = null;
  for (let i = 0; i < us.length; i++) {
    const v = vs[i];
    if (!isFinite(v)) continue;
    const x = A.px(us[i]), y = A.py(v);
    if (!started) { ctx.moveTo(x, A.py(base)); started = true; }
    ctx.lineTo(x, y); last = x;
  }
  if (started) { ctx.lineTo(last, A.py(base)); ctx.closePath(); ctx.fill(); }
  ctx.restore();
}
function band(ctx, A, us, lo, hi, color, alpha) {
  clip(ctx, A);
  ctx.globalAlpha = alpha; ctx.fillStyle = color;
  ctx.beginPath();
  for (let i = 0; i < us.length; i++) { const x = A.px(us[i]), y = A.py(hi[i]); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
  for (let i = us.length - 1; i >= 0; i--) ctx.lineTo(A.px(us[i]), A.py(lo[i]));
  ctx.closePath(); ctx.fill();
  ctx.restore();
}
function vline(ctx, A, u, color, dash, width = 1) {
  const x = Math.round(A.px(u)) + 0.5, b = A.box;
  ctx.save(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []);
  ctx.beginPath(); ctx.moveTo(x, b.t); ctx.lineTo(x, b.b); ctx.stroke(); ctx.restore();
}
function cursorHandle(ctx, A, u) {
  const x = Math.round(A.px(u)) + 0.5, b = A.box;
  ctx.fillStyle = C.muted;
  ctx.beginPath(); ctx.moveTo(x, b.b - 8); ctx.lineTo(x - 5, b.b); ctx.lineTo(x + 5, b.b); ctx.closePath(); ctx.fill();
}
function legend(el, items) {
  el.innerHTML = items.map(it => {
    const cls = it.band ? "sw band" : "sw" + (it.dash ? " dash" : "") + (it.dot ? " dot" : "");
    const st = it.band ? `background:${it.color};opacity:${it.alpha || 0.3}` : `border-top-color:${it.color}`;
    return `<span title="${it.tip || ""}"><span class="${cls}" style="${st}"></span><span>${it.label}</span></span>`;
  }).join("");
}
function local(cv, ev) { const r = cv.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; }

// ------------------------------------------------------------------ wavelength helpers
let LAMLOG = null;           // log10 of the wavelength grid
function lamIndex(l10) {
  const n = LAMLOG.length, u = (l10 - LAMLOG[0]) / (LAMLOG[1] - LAMLOG[0]);
  return clamp(Math.round(u), 0, n - 1);
}
function atLam(arr, lam) {     // log-log interpolation at a wavelength in um
  const l = Math.log10(lam), n = LAMLOG.length;
  const u = (l - LAMLOG[0]) / (LAMLOG[1] - LAMLOG[0]);
  const k = clamp(Math.floor(u), 0, n - 2), t = u - k;
  if (arr[k] > 0 && arr[k + 1] > 0) return Math.exp((1 - t) * Math.log(arr[k]) + t * Math.log(arr[k + 1]));
  return (1 - t) * arr[k] + t * arr[k + 1];
}

// ------------------------------------------------------------------ single grain
let G1 = null;
function computeGrain() {
  const a = Math.pow(10, S.a), U = Math.pow(10, S.U);
  const opt = P.optics(S.mat, a);
  const g = new P.Grain(S.mat, a, U, opt);
  const r = g.emission();
  const [Rph, S1] = g.photonStats();
  const n = P.NL, eq = new Float64Array(n);
  for (let i = 0; i < n; i++) eq[i] = 4 * Math.PI * g.Cabs[i] * P.planck(i, g.Teq);
  const Eeq = g.H.E(g.Teq), hnu = S1 / Rph;
  const pct = P.percentiles(r.T, r.P, [0.05, 0.95]);
  G1 = { a, U, mat: S.mat, opt, g, r, eq, Rph, hnu, Eeq, Tjump: g.H.T(Eeq + hnu), pct };
}

// ------------------------------------------------------------------ per-size nodes (shared)
// Emission of one grain at each node of the size grid, for a material and field, cached.
// The temperature panel and the size-distribution emission both draw on these.  They are
// computed in time slices so the page stays responsive while U is dragged.
const NODE = new Map();
let queue = [], running = false;
const ukey = () => S.U.toFixed(2);
const nkey = (m, i, uk) => m + "|" + i + "|" + uk;
function neededNodes() {
  const list = [];
  if (S.tab === "grain") {
    const ia = Math.round((S.a + 3) * 25);
    for (let i = T_NODES[0]; i <= T_NODES[1]; i++) list.push([S.mat, i]);
    list.sort((p, q) => Math.abs(p[1] - ia) - Math.abs(q[1] - ia));
  } else {
    for (const m of ["sil", "gra"]) {
      const d = S.dist[m];
      if (!d.on) continue;
      const [i0, i1] = P.nodeRange(Math.pow(10, d.amin), Math.pow(10, d.amax));
      for (let i = i0; i <= i1; i++) list.push([m, i]);
    }
  }
  return list;
}
function requestNodes() {
  const uk = ukey();
  queue = neededNodes().filter(([m, i]) => !NODE.has(nkey(m, i, uk))).map(([m, i]) => [m, i, uk]);
  if (NODE.size > 4000) for (const k of [...NODE.keys()]) if (!k.endsWith("|" + uk)) NODE.delete(k);
  if (queue.length && !running) { running = true; setTimeout(work, 0); }
  pendingNote();
}
function work() {
  const t0 = performance.now();
  while (queue.length && performance.now() - t0 < 14) {
    const [m, i, uk] = queue.shift();
    const k = nkey(m, i, uk);
    if (!NODE.has(k)) NODE.set(k, P.nodeEmission(m, i, Math.pow(10, parseFloat(uk))));
  }
  if (queue.length) setTimeout(work, 0); else running = false;
  if (S.tab === "grain") drawT(); else if (!queue.length) { computeDistEm(); drawI(); readDist(); }
  pendingNote();
}
function pendingNote() {
  const busy = queue.length > 0;
  $("tSub").textContent = busy && S.tab === "grain" ? "computing…" : "";
  $("iSub").textContent = busy && S.tab === "dist" ? "computing…" : "";
}

// ------------------------------------------------------------------ size distribution
let DO = null, DE = null;
function distParams(m) {
  const d = S.dist[m];
  return { q: d.q, amin: Math.pow(10, d.amin), amax: Math.pow(10, d.amax), abund: d.abund };
}
const onMats = () => ["sil", "gra"].filter(m => S.dist[m].on && S.dist[m].abund > 0);
function computeDistOptics() {
  const out = { mats: onMats(), ext: {}, sca: {}, mass: {} };
  const n = P.NL;
  out.extT = new Float64Array(n); out.scaT = new Float64Array(n);
  out.extB = P.SIZE_BINS.map(() => new Float64Array(n)); out.scaB = P.SIZE_BINS.map(() => new Float64Array(n));
  for (const m of out.mats) {
    const o = P.distOptics(m, distParams(m));
    out.ext[m] = o.ext; out.sca[m] = o.sca; out.mass[m] = o.mass;
    for (let j = 0; j < n; j++) {
      out.extT[j] += o.ext[j]; out.scaT[j] += o.sca[j];
      for (let k = 0; k < out.extB.length; k++) { out.extB[k][j] += o.extB[k][j]; out.scaB[k][j] += o.scaB[k][j]; }
    }
  }
  DO = out;
}
// emission per H, if every node it needs is ready; otherwise the last result stays, marked stale
function computeDistEm() {
  const uk = ukey(), n = P.NL, mats = onMats();
  const res = { key: "", mats, em: {}, bins: [0, 1, 2].map(() => new Float64Array(n)), tot: new Float64Array(n), pabs: 0 };
  for (const m of mats) {
    const p = distParams(m), [i0, i1] = P.nodeRange(p.amin, p.amax), nodes = [];
    for (let i = i0; i <= i1; i++) {
      const v = NODE.get(nkey(m, i, uk));
      if (!v) { if (DE) DE.stale = true; return; }
      nodes[i] = v;
    }
    const e = P.distEmission(m, p, nodes);
    res.em[m] = e.em; res.pabs += e.pabs;
    for (let j = 0; j < n; j++) { res.tot[j] += e.em[j]; for (let k = 0; k < 3; k++) res.bins[k][j] += e.bins[k][j]; }
  }
  res.stale = false;
  DE = res;
}

// ------------------------------------------------------------------ hover state
let HOV = null;        // {grp, u} in plot x units of the hovered panel group
const tip = $("tip");
function showTip(ev, rows) {
  tip.innerHTML = rows.map(r => r.head ? `<div style="margin-bottom:2px">${r.head}</div>`
    : `<div class="r"><span>${r.sw ? `<span class="sw" style="border-top-color:${r.sw}${r.dash ? ";border-top-style:dashed" : ""}"></span>` : ""}${r.k}</span><b>${r.v}</b></div>`).join("");
  tip.hidden = false;
  const w = tip.offsetWidth, h = tip.offsetHeight;
  let x = ev.clientX + 16, y = ev.clientY + 14;
  if (x + w > window.innerWidth - 8) x = ev.clientX - w - 16;
  if (y + h > window.innerHeight - 8) y = ev.clientY - h - 14;
  tip.style.left = x + "px"; tip.style.top = y + "px";
}
function hideTip() { tip.hidden = true; }

// ================================================================== single-grain panels
const LX = [Math.log10(0.09), Math.log10(20)];  // wavelength axis, efficiencies (to 20 um)
const LXE = [Math.log10(0.09), Math.log10(3000)];
let AQ = null, AW = null, AE = null, AP = null, AT = null;

function drawQ() {
  const { ctx, w, h } = fit($("cvQ"));
  const A = axes({ l: 74, r: w - 20, t: 10, b: h - 6 }, LX[0], LX[1], -7, 1);
  AQ = A;
  const o = G1.opt, n = P.NL;
  const qe = new Float64Array(n), qa = new Float64Array(n), qs = new Float64Array(n);
  for (let i = 0; i < n; i++) { qa[i] = L10(o.Qabs[i]); qs[i] = L10(o.Qsca[i]); qe[i] = L10(o.Qabs[i] + o.Qsca[i]); }
  frame(ctx, A);
  logTicks(ctx, A, "x", { labels: false });
  logTicks(ctx, A, "y", { every: 2 });
  ylabel(ctx, A, mk("Efficiency *Q*"));
  if (HOV && HOV.grp === "q") vline(ctx, A, HOV.u, C.faint);
  curve(ctx, A, LAMLOG, qs, C.sca, 2, [7, 4]);
  curve(ctx, A, LAMLOG, qa, C.abs, 2);
  curve(ctx, A, LAMLOG, qe, C.ink, 2.2);
}
function drawW() {
  const { ctx, w, h } = fit($("cvW"));
  const A = axes({ l: 74, r: w - 20, t: 6, b: h - 46 }, LX[0], LX[1], 0, 1);
  AW = A;
  const o = G1.opt, n = P.NL, al = new Float64Array(n);
  for (let i = 0; i < n; i++) al[i] = o.Qsca[i] / (o.Qabs[i] + o.Qsca[i]);
  frame(ctx, A);
  logTicks(ctx, A, "x", { plain: true });
  linTicks(ctx, A, "y", 0.5);
  xlabel(ctx, A, mk("Wavelength λ (µm)"));
  ylabel(ctx, A, mk("Albedo, *g*"));
  if (HOV && HOV.grp === "q") vline(ctx, A, HOV.u, C.faint);
  curve(ctx, A, LAMLOG, o.g, C.ink, 1.6, [2, 3]);
  curve(ctx, A, LAMLOG, al, C.sca, 2);
}
function hoverQ(ev, A, cv) {
  const [x, y] = local(cv, ev);
  if (!A || x < A.box.l || x > A.box.r) { HOV = null; hideTip(); drawQ(); drawW(); return; }
  const u = A.ux(x), i = lamIndex(u), o = G1.opt;
  HOV = { grp: "q", u: LAMLOG[i] };
  drawQ(); drawW();
  const qe = o.Qabs[i] + o.Qsca[i];
  showTip(ev, [{ head: `λ = ${numHTML(P.LAM[i], 3)} µm` },
    { k: "<i>Q</i><sub>ext</sub>", v: numHTML(qe, 3), sw: C.ink }, { k: "<i>Q</i><sub>abs</sub>", v: numHTML(o.Qabs[i], 3), sw: C.abs },
    { k: "<i>Q</i><sub>sca</sub>", v: numHTML(o.Qsca[i], 3), sw: C.sca, dash: true },
    { k: "albedo", v: (o.Qsca[i] / qe).toFixed(3), sw: C.sca }, { k: "<i>g</i>", v: o.g[i].toFixed(3), sw: C.ink, dash: true }]);
  void y;
}

function drawE() {
  const { ctx, w, h } = fit($("cvE"));
  const n = P.NL, g = G1.g, m = g.mass;
  const ab = new Float64Array(n), em = new Float64Array(n), eq = new Float64Array(n);
  let mx = 0;
  for (let i = 0; i < n; i++) {
    const l = P.LAM[i];
    ab[i] = l * g.pabsLam[i] / m; em[i] = l * G1.r.em[i] / m; eq[i] = l * G1.eq[i] / m;
    mx = Math.max(mx, ab[i], em[i]);
  }
  let A, ya, ye, yq;
  if (S.elog) {
    const top = Math.ceil(Math.log10(mx) * 2 + 0.3) / 2;
    A = axes({ l: 80, r: w - 14, t: 10, b: h - 46 }, LXE[0], LXE[1], top - 6, top);
    ya = Array.from(ab, L10); ye = Array.from(em, L10); yq = Array.from(eq, L10);
  } else {
    let e = Math.floor(Math.log10(mx));
    if (e >= -1 && e <= 2) e = 0;
    const sc = Math.pow(10, e);
    A = axes({ l: 80, r: w - 14, t: 10, b: h - 46 }, LXE[0], LXE[1], 0, 1.08 * mx / sc);
    A.scale = e;
    ya = Array.from(ab, v => v / sc); ye = Array.from(em, v => v / sc); yq = Array.from(eq, v => v / sc);
  }
  AE = A;
  frame(ctx, A);
  logTicks(ctx, A, "x", { plain: true });
  if (S.elog) { logTicks(ctx, A, "y"); ylabel(ctx, A, mk("λ*P*_{λ} per gram (erg s^{−1} g^{−1})"), 58); }
  else {
    linTicks(ctx, A, "y", niceStep(A.y1, 5));
    ylabel(ctx, A, A.scale === 0 ? mk("λ*P*_{λ} per gram (erg s^{−1} g^{−1})")
      : [...mk("λ*P*_{λ} per gram (10"), { t: expStr(A.scale), sup: true }, ...mk(" erg s^{−1} g^{−1})")], 58);
  }
  xlabel(ctx, A, mk("Wavelength λ (µm)"));
  if (HOV && HOV.grp === "e") vline(ctx, A, HOV.u, C.faint);
  const base = S.elog ? A.y0 : 0;
  fillTo(ctx, A, LAMLOG, ya, base, C.abs, 0.10);
  fillTo(ctx, A, LAMLOG, ye, base, C.emi, 0.12);
  if (G1.r.stochastic) curve(ctx, A, LAMLOG, yq, C.emi, 1.6, [6, 4]);
  curve(ctx, A, LAMLOG, ya, C.abs, 2);
  curve(ctx, A, LAMLOG, ye, C.emi, 2.2);
  const items = [{ color: C.abs, label: "absorbed", tip: "starlight and CMB absorbed, c u_λ C_abs" },
                 { color: C.emi, label: G1.r.stochastic ? "emitted (stochastic heating)" : "emitted (at <i>T</i><sub>eq</sub>)" }];
  if (G1.r.stochastic) items.push({ color: C.emi, dash: true, label: "emitted if held at <i>T</i><sub>eq</sub>",
    tip: "the spectrum the same grain would emit at the single temperature T_eq" });
  legend($("lgE"), items);
}
function hoverE(ev) {
  const cv = $("cvE"), [x] = local(cv, ev), A = AE;
  if (!A || x < A.box.l || x > A.box.r) { HOV = null; hideTip(); drawE(); return; }
  const i = lamIndex(A.ux(x)), g = G1.g, l = P.LAM[i];
  HOV = { grp: "e", u: LAMLOG[i] };
  drawE();
  const rows = [{ head: `λ = ${numHTML(l, 3)} µm` },
    { k: "absorbed", v: numHTML(l * g.pabsLam[i] / g.mass, 3), sw: C.abs },
    { k: "emitted", v: numHTML(l * G1.r.em[i] / g.mass, 3), sw: C.emi }];
  if (G1.r.stochastic) rows.push({ k: "at <i>T</i><sub>eq</sub>", v: numHTML(l * G1.eq[i] / g.mass, 3), sw: C.emi, dash: true });
  rows.push({ head: '<span style="color:#6a6560">erg s<sup>−1</sup> g<sup>−1</sup>, as λ<i>P</i><sub>λ</sub></span>' });
  showTip(ev, rows);
}

function drawP() {
  const { ctx, w, h } = fit($("cvP"));
  const r = G1.r, Teq = G1.g.Teq;
  let lt0, lt1, A;
  const box = { l: 74, r: w - 14, t: 10, b: h - 46 };
  if (!r.stochastic) {
    lt0 = Math.log10(Teq) - 0.3; lt1 = Math.log10(Teq) + 0.3;
    A = axes(box, lt0, lt1, S.plin ? 0 : -6, S.plin ? 1.08 : 0.3);
  } else {
    const nT = r.T.length, dl = Math.log(r.T[1] / r.T[0]);
    const dens = Array.from(r.P, p => p / dl);
    const pk = Math.max(...dens);
    let j0 = 0, j1 = nT - 1;
    const thr = S.plin ? 1e-3 : 1e-6;
    while (j0 < nT - 1 && dens[j0] < thr * pk) j0++;
    while (j1 > 0 && dens[j1] < thr * pk) j1--;
    lt0 = Math.log10(Math.min(r.T[j0], Teq)) - 0.06; lt1 = Math.log10(Math.max(r.T[j1], Teq)) + 0.06;
    if (lt1 - lt0 < 0.3) { const c = (lt0 + lt1) / 2; lt0 = c - 0.15; lt1 = c + 0.15; }
    A = S.plin ? axes(box, lt0, lt1, 0, 1.08 * pk) : axes(box, lt0, lt1, Math.log10(pk) - 6, Math.log10(pk) + 0.3);
    A.dens = dens;
  }
  AP = A;
  frame(ctx, A);
  logTicks(ctx, A, "x", { plain: true });
  if (r.stochastic) {
    if (S.plin) linTicks(ctx, A, "y", niceStep(A.y1, 4)); else logTicks(ctx, A, "y");
    ylabel(ctx, A, mk("d*P* / d ln *T*"));
  }
  xlabel(ctx, A, mk("Temperature *T* (K)"));
  vline(ctx, A, Math.log10(Teq), C.ink, [4, 3]);
  if (!r.stochastic) {
    vline(ctx, A, Math.log10(Teq), C.emi, [], 3);
  } else {
    const us = Array.from(r.T, t => Math.log10(t));
    const vs = S.plin ? A.dens : A.dens.map(L10);
    fillTo(ctx, A, us, vs, S.plin ? 0 : A.y0, C.emi, 0.12);
    curve(ctx, A, us, vs, C.emi, 2);
  }
  if (HOV && HOV.grp === "p") vline(ctx, A, HOV.u, C.faint);
  legend($("lgP"), r.stochastic
    ? [{ color: C.emi, label: "P(T)" }, { color: C.ink, dash: true, label: "<i>T</i><sub>eq</sub>" }]
    : [{ color: C.emi, label: "one steady temperature, <i>T</i><sub>eq</sub> (fluctuations below 1%)" }]);
}
function hoverP(ev) {
  const cv = $("cvP"), [x] = local(cv, ev), A = AP;
  if (!A || x < A.box.l || x > A.box.r || !G1.r.stochastic) { if (HOV && HOV.grp === "p") { HOV = null; drawP(); } hideTip(); return; }
  const lt = A.ux(x), r = G1.r;
  let j = 0, best = Infinity;
  for (let k = 0; k < r.T.length; k++) { const d = Math.abs(Math.log10(r.T[k]) - lt); if (d < best) { best = d; j = k; } }
  HOV = { grp: "p", u: Math.log10(r.T[j]) };
  drawP();
  let cum = 0;
  for (let k = j; k < r.T.length; k++) cum += r.P[k];
  showTip(ev, [{ head: `<i>T</i> = ${tHTML(r.T[j])}` }, { k: "d<i>P</i>/d ln <i>T</i>", v: numHTML(A.dens[j], 3), sw: C.emi },
    { k: "time spent above", v: numHTML(100 * cum, 2) + "%" }]);
}

function drawT() {
  const { ctx, w, h } = fit($("cvT"));
  const uk = ukey(), us = [], te = [], lo = [], hi = [], tj = [];
  for (let i = T_NODES[0]; i <= T_NODES[1]; i++) {
    const v = NODE.get(nkey(S.mat, i, uk));
    us.push(Math.log10(P.A_GRID[i]));
    te.push(v ? Math.log10(v.Teq) : NaN); lo.push(v ? Math.log10(v.t05) : NaN);
    hi.push(v ? Math.log10(v.t95) : NaN); tj.push(v ? Math.log10(v.Tjump) : NaN);
  }
  let ymin = Infinity, ymax = -Infinity;
  for (let k = 0; k < us.length; k++) {
    if (isFinite(lo[k])) ymin = Math.min(ymin, lo[k], te[k]);
    if (isFinite(tj[k])) ymax = Math.max(ymax, tj[k], hi[k]);
  }
  if (!isFinite(ymin)) { ymin = Math.log10(G1.g.Teq) - 0.5; ymax = Math.log10(G1.g.Teq) + 1; }
  ymin = Math.max(Math.log10(2), ymin - 0.08); ymax = Math.min(Math.log10(5000), ymax + 0.08);
  const A = axes({ l: 74, r: w - 14, t: 10, b: h - 46 }, R.a[0], R.a[1], ymin, ymax);
  AT = A;
  A.data = { us, te, lo, hi, tj };
  frame(ctx, A);
  logTicks(ctx, A, "x", { plain: true });
  logTicks(ctx, A, "y", { plain: true });
  xlabel(ctx, A, mk("Grain radius *a* (µm)"));
  ylabel(ctx, A, mk("Temperature *T* (K)"));
  // the band only where the nodes are ready
  let k0 = 0;
  while (k0 < us.length) {
    if (!isFinite(lo[k0])) { k0++; continue; }
    let k1 = k0;
    while (k1 + 1 < us.length && isFinite(lo[k1 + 1])) k1++;
    band(ctx, A, us.slice(k0, k1 + 1), lo.slice(k0, k1 + 1), hi.slice(k0, k1 + 1), C.emi, 0.2);
    k0 = k1 + 1;
  }
  curve(ctx, A, us, tj, C.ink, 1.8, [2, 3]);
  curve(ctx, A, us, te, C.emi, 2.2);
  vline(ctx, A, S.a, C.muted);
  cursorHandle(ctx, A, S.a);
  legend($("lgT"), [{ color: C.emi, label: "<i>T</i><sub>eq</sub>" },
    { band: true, color: C.emi, alpha: 0.25, label: "middle 90% of P(T)" },
    { color: C.ink, dot: true, label: "after one photon of mean energy" }]);
}
function hoverT(ev) {
  const cv = $("cvT"), [x] = local(cv, ev), A = AT;
  if (!A || x < A.box.l || x > A.box.r) { hideTip(); return; }
  const d = A.data, k = clamp(Math.round((A.ux(x) - d.us[0]) / (d.us[1] - d.us[0])), 0, d.us.length - 1);
  if (!isFinite(d.te[k])) { hideTip(); return; }
  const p = v => tHTML(Math.pow(10, v));
  showTip(ev, [{ head: `<i>a</i> = ${aHTML(Math.pow(10, d.us[k]))}` },
    { k: "<i>T</i><sub>eq</sub>", v: p(d.te[k]), sw: C.emi }, { k: "middle 90%", v: p(d.lo[k]) + " – " + p(d.hi[k]) },
    { k: "after one photon", v: p(d.tj[k]), sw: C.ink, dash: true }]);
}
let tDrag = false;
function setAfromT(ev) {
  const [x] = local($("cvT"), ev);
  if (!AT) return;
  S.a = clamp(AT.ux(x), R.a[0], R.a[1]);
  $("sA").value = S.a;
  changed("grain");
}

function readGrain() {
  const g = G1.g, o = G1.opt;
  const qV = atLam(Float64Array.from(o.Qabs, (v, i) => v + o.Qsca[i]), KEY_UM.V);
  const sV = atLam(o.Qsca, KEY_UM.V);
  const rows = [
    ["Atoms in the grain", numHTML(g.H.N, 2), "N, from the bulk density (MgFeSiO4 at 3.5 g/cm3; graphite at 2.24 g/cm3)"],
    ["<i>Q</i><sub>ext</sub> at 0.55 µm", numHTML(qV, 3), "extinction efficiency in the V band"],
    ["Albedo at 0.55 µm", (sV / qV).toFixed(3), "fraction of the extinction in the V band that is scattering"],
    ["Equilibrium <i>T</i><sub>eq</sub>", tHTML(g.Teq), "the temperature at which emission balances absorption"],
    ["Middle 90% of <i>T</i>", G1.r.stochastic ? tHTML(G1.pct[0]) + " – " + tHTML(G1.pct[1]) : "steady",
      "5th to 95th percentile of the temperature distribution"],
    ["Absorbed power", numHTML(g.Pabs, 3) + " erg s<sup>−1</sup>", "per grain; equal to the emitted power"],
    ["Starlight photon absorbed every", durHTML(1 / G1.Rph), "the mean interval between absorptions of starlight photons"],
    ["Mean photon energy absorbed", numHTML(G1.hnu / P.EV, 3) + " eV", "absorbed starlight power divided by the photon absorption rate"],
    ["Thermal energy at <i>T</i><sub>eq</sub>", numHTML(G1.Eeq / P.EV, 3) + " eV", "vibrational energy content E(T_eq), Debye heat capacities of Draine & Li (2001)"],
    ["Radiative cooling time", durHTML(G1.Eeq / g.Pabs), "time to radiate the thermal energy at T_eq: E(T_eq) / P_abs. Compare with the photon interval."],
    ["<i>T</i> after one such photon", tHTML(G1.Tjump), "the temperature reached by adding the mean photon energy to E(T_eq)"],
  ];
  $("readG").innerHTML = rows.map((r, k) =>
    `<div class="rd${k === 3 || k === 5 ? " gap" : ""}" title="${r[2]}"><span>${r[0]}</span><b>${r[1]}</b></div>`).join("");
}

// ================================================================== distribution panels
const LXA = [Math.log10(0.002), Math.log10(5)];
let AS = null, AX = null, AA = null, AI = null;

function sizeCurve(m) {
  const d = distParams(m), A = P.powerlawA(m, d.q, d.amin, d.amax, d.abund);
  const n = 160, us = [], vs = [];
  const l0 = Math.log10(d.amin), l1 = Math.log10(d.amax);
  for (let k = 0; k < n; k++) {
    const l = l0 + (l1 - l0) * k / (n - 1), a = Math.pow(10, l);
    const dn = A * Math.pow(a, -d.q);                          // per um per H
    const per = S.sd === "mass" ? P.grainMass(m, a) : Math.PI * Math.pow(a * P.UM, 2);
    us.push(l); vs.push(LN10 * a * dn * per);
  }
  return { us, vs };
}
function drawS() {
  const { ctx, w, h } = fit($("cvS"));
  const cs = {};
  let mx = 0;
  for (const m of onMats()) { cs[m] = sizeCurve(m); mx = Math.max(mx, ...cs[m].vs); }
  const top = mx > 0 ? Math.ceil(Math.log10(mx) * 2 + 0.4) / 2 : -26;
  const A = axes({ l: 80, r: w - 14, t: 10, b: h - 46 }, LXA[0], LXA[1], top - 4, top);
  AS = A; A.cs = cs;
  frame(ctx, A);
  logTicks(ctx, A, "x", { plain: true });
  logTicks(ctx, A, "y");
  xlabel(ctx, A, mk("Grain radius *a* (µm)"));
  ylabel(ctx, A, S.sd === "mass" ? mk("Mass per decade of *a*, per H (g)") : mk("Area per decade of *a*, per H (cm^{2})"), 58);
  if (HOV && HOV.grp === "s") vline(ctx, A, HOV.u, C.faint);
  for (const m of onMats()) {
    const c = cs[m], ys = c.vs.map(L10);
    // drop to the axis at both ends: the distribution is cut off there
    const us = [c.us[0], ...c.us, c.us[c.us.length - 1]], vs = [A.y0 - 1, ...ys, A.y0 - 1];
    fillTo(ctx, A, c.us, ys, A.y0, C[m], 0.08);
    curve(ctx, A, us, vs, C[m], 2.2, m === "gra" ? [7, 4] : []);
  }
  legend($("lgS"), onMats().map(m => ({ color: C[m], dash: m === "gra", label: MATNAME[m] })));
}
function hoverS(ev) {
  const cv = $("cvS"), [x] = local(cv, ev), A = AS;
  if (!A || x < A.box.l || x > A.box.r) { if (HOV && HOV.grp === "s") { HOV = null; drawS(); } hideTip(); return; }
  const l = A.ux(x);
  HOV = { grp: "s", u: l };
  drawS();
  const rows = [{ head: `<i>a</i> = ${aHTML(Math.pow(10, l))}` }];
  for (const m of onMats()) {
    const c = A.cs[m];
    let v = 0;
    if (l >= c.us[0] && l <= c.us[c.us.length - 1]) {
      const k = clamp(Math.round((l - c.us[0]) / (c.us[1] - c.us[0])), 0, c.us.length - 1);
      v = c.vs[k];
    }
    rows.push({ k: MATNAME[m], v: v > 0 ? numHTML(v, 3) + (S.sd === "mass" ? " g" : " cm<sup>2</sup>") : "none", sw: C[m], dash: m === "gra" });
  }
  showTip(ev, rows);
}

const XI = [0, 11];                     // inverse wavelength axis (um^-1)
const SIZELAB = ["<i>a</i> < 0.01 µm", "0.01–0.1 µm", "> 0.1 µm"];
// legend entry for a component: silicate dotted, graphite dashed, size bins solid
const compLegend = c => ({ color: c.color, dash: c.key === "gra", dot: c.key === "sil", label: c.label });
function extArrays() {
  const mag = 2.5 / Math.LN10;
  const AV = mag * atLam(DO.extT, KEY_UM.V);
  const norm = S.xn === "H" ? 1e-21 : AV;
  const f = v => mag * v / norm;
  const inv = Array.from(P.LAM, l => 1 / l);
  const tot = Array.from(DO.extT, f);
  const alb = Array.from(DO.scaT, (v, i) => v / DO.extT[i]);
  const comps = S.xs === "mat"
    ? DO.mats.map(m => ({ key: m, label: MATNAME[m], color: C[m], dash: m === "gra" ? [7, 4] : [2, 3], ext: DO.ext[m], sca: DO.sca[m] }))
    : DO.extB.map((e, k) => ({ key: "b" + k, label: SIZELAB[k], color: C.size[k], dash: [], ext: e, sca: DO.scaB[k] }))
        .filter(c => c.ext.some(v => v > 0));
  for (const c of comps) {
    c.y = Array.from(c.ext, f);
    c.alb = Array.from(c.sca, (v, i) => (c.ext[i] > 0 ? v / c.ext[i] : NaN));
  }
  return { inv, tot, alb, comps, AV };
}
function drawX() {
  const { ctx, w, h } = fit($("cvX"));
  const box = { l: 74, r: w - 18, t: 10, b: h - 6 };
  if (!DO.mats.length) { const A = axes(box, XI[0], XI[1], 0, 1); frame(ctx, A); AX = null; legend($("lgX"), []); return; }
  const e = extArrays();
  const mx = Math.max(...e.tot);
  const A = axes(box, XI[0], XI[1], 0, 1.08 * mx);
  AX = A; A.e = e;
  frame(ctx, A);
  linTicks(ctx, A, "x", 1, false);
  linTicks(ctx, A, "y", niceStep(A.y1, 5));
  ylabel(ctx, A, S.xn === "H" ? mk("*A*_{λ} / *N*_{H} (10^{−21} mag cm^{2})") : mk("*A*_{λ} / *A*_{V}"));
  if (HOV && HOV.grp === "x") vline(ctx, A, HOV.u, C.faint);
  for (const c of e.comps) curve(ctx, A, e.inv, c.y, c.color, 1.8, c.dash);
  curve(ctx, A, e.inv, e.tot, C.ink, 2.4);
  legend($("lgX"), [{ color: C.ink, label: "total" }, ...e.comps.map(compLegend)]);
}
function drawA() {
  const { ctx, w, h } = fit($("cvA"));
  const A = axes({ l: 74, r: w - 18, t: 6, b: h - 46 }, XI[0], XI[1], 0, 1);
  AA = A;
  frame(ctx, A);
  linTicks(ctx, A, "x", 1);
  linTicks(ctx, A, "y", 0.5);
  xlabel(ctx, A, mk("Inverse wavelength 1/λ (µm^{−1})"));
  ylabel(ctx, A, mk("Albedo"));
  if (!DO.mats.length) return;
  const e = AX ? AX.e : extArrays();
  if (HOV && HOV.grp === "x") vline(ctx, A, HOV.u, C.faint);
  for (const c of e.comps) curve(ctx, A, e.inv, c.alb, c.color, 1.6, c.dash);
  curve(ctx, A, e.inv, e.alb, C.ink, 2.2);
}
function hoverX(ev, cv, A) {
  const [x] = local(cv, ev);
  if (!A || !AX || x < A.box.l || x > A.box.r) { if (HOV && HOV.grp === "x") { HOV = null; drawX(); drawA(); } hideTip(); return; }
  const inv = clamp(A.ux(x), 1 / P.LAM[P.NL - 1], 1 / P.LAM[0]);
  const i = lamIndex(-Math.log10(inv)), e = AX.e;
  HOV = { grp: "x", u: e.inv[i] };
  drawX(); drawA();
  const unit = S.xn === "H" ? " × 10<sup>−21</sup>" : "";
  const rows = [{ head: `λ = ${numHTML(P.LAM[i], 3)} µm &nbsp; (1/λ = ${numHTML(e.inv[i], 3)} µm<sup>−1</sup>)` },
    { k: S.xn === "H" ? "<i>A</i><sub>λ</sub>/<i>N</i><sub>H</sub> (mag cm<sup>2</sup>)" : "<i>A</i><sub>λ</sub>/<i>A</i><sub>V</sub>", v: numHTML(e.tot[i], 3) + unit, sw: C.ink }];
  for (const c of e.comps) rows.push({ k: "&nbsp;&nbsp;" + c.label, v: numHTML(c.y[i], 3) + unit, sw: c.color, dash: c.key === "gra" });
  rows.push({ k: "albedo", v: e.alb[i].toFixed(3), sw: C.ink });
  for (const c of e.comps) if (isFinite(c.alb[i])) rows.push({ k: "&nbsp;&nbsp;" + c.label, v: c.alb[i].toFixed(3), sw: c.color, dash: c.key === "gra" });
  showTip(ev, rows);
}

const LXI = [0, 3];
function drawI() {
  const { ctx, w, h } = fit($("cvI"));
  const box = { l: 80, r: w - 22, t: 10, b: h - 46 };
  if (!DE || !DE.mats.length || !onMats().length) {
    const A = axes(box, LXI[0], LXI[1], -30, -26); frame(ctx, A); AI = null; legend($("lgI"), []); return;
  }
  const series = S.sp === "mat"
    ? DE.mats.map(m => ({ y: DE.em[m], color: C[m], dash: m === "gra" ? [7, 4] : [2, 3], label: MATNAME[m], key: m }))
    : DE.bins.map((y, k) => ({ y, color: C.size[k], dash: [], label: SIZELAB[k], key: "b" + k })).filter(sr => sr.y.some(v => v > 0));
  let mx = 0;
  for (let j = 0; j < P.NL; j++) if (P.LAM[j] >= 1) mx = Math.max(mx, DE.tot[j]);
  const top = Math.ceil(Math.log10(mx) * 2 + 0.4) / 2;
  const A = axes(box, LXI[0], LXI[1], top - 6, top);
  AI = A; A.series = series;
  const alpha = DE.stale ? 0.35 : 1;
  frame(ctx, A);
  logTicks(ctx, A, "x", { plain: true });
  logTicks(ctx, A, "y");
  xlabel(ctx, A, mk("Wavelength λ (µm)"));
  ylabel(ctx, A, mk("λ*I*_{λ} / *N*_{H} (erg s^{−1} sr^{−1} H^{−1})"), 58);
  if (HOV && HOV.grp === "i") vline(ctx, A, HOV.u, C.faint);
  for (const s of series) curve(ctx, A, LAMLOG, Array.from(s.y, L10), s.color, 1.9, s.dash, alpha);
  curve(ctx, A, LAMLOG, Array.from(DE.tot, L10), C.ink, 2.4, [], alpha);
  legend($("lgI"), [{ color: C.ink, label: "total" }, ...series.map(s => ({ color: s.color, dash: s.dash.length && s.key === "gra", dot: s.key === "sil", label: s.label }))]);
}
function hoverI(ev) {
  const cv = $("cvI"), [x] = local(cv, ev), A = AI;
  if (!A || x < A.box.l || x > A.box.r) { if (HOV && HOV.grp === "i") { HOV = null; drawI(); } hideTip(); return; }
  const i = lamIndex(A.ux(x));
  HOV = { grp: "i", u: LAMLOG[i] };
  drawI();
  const rows = [{ head: `λ = ${numHTML(P.LAM[i], 3)} µm` }, { k: "total", v: numHTML(DE.tot[i], 3), sw: C.ink }];
  for (const s of A.series) rows.push({ k: s.label, v: numHTML(s.y[i], 3), sw: s.color, dash: s.key === "gra" });
  rows.push({ head: '<span style="color:#6a6560">erg s<sup>−1</sup> sr<sup>−1</sup> H<sup>−1</sup></span>' });
  showTip(ev, rows);
}

function readDist() {
  const rows = [];
  if (!DO.mats.length) { $("readD").innerHTML = '<div class="rd"><span>No grains</span><b></b></div>'; return; }
  const mag = 2.5 / Math.LN10;
  const AV = mag * atLam(DO.extT, KEY_UM.V), AB = mag * atLam(DO.extT, KEY_UM.B);
  const alb = atLam(DO.scaT, KEY_UM.V) / atLam(DO.extT, KEY_UM.V);
  const mass = DO.mats.reduce((s, m) => s + DO.mass[m], 0);
  rows.push(["<i>A</i><sub>V</sub> / <i>N</i><sub>H</sub>", numHTML(AV, 3) + " mag cm<sup>2</sup>", "visual extinction per H nucleus, at 0.55 um"]);
  rows.push(["<i>R</i><sub>V</sub> = <i>A</i><sub>V</sub> / <i>E</i>(<i>B</i>−<i>V</i>)", (AV / (AB - AV)).toFixed(2), "evaluated at 0.44 and 0.55 um"]);
  rows.push(["Albedo at 0.55 µm", alb.toFixed(3), "scattering / extinction in the V band"]);
  rows.push(["Dust mass / H mass", numHTML(mass / MH, 3), "grain mass per H nucleus divided by the mass of an H atom"]);
  if (DE && !DE.stale) {
    let pw = 0, pshort = 0, psmall = 0;
    for (let j = 0; j < P.NL; j++) {
      const wj = P.WL[j] / P.LAM[j];                 // d ln lambda weight
      pw += wj * DE.tot[j]; psmall += wj * DE.bins[0][j];
      if (P.LAM[j] < 50) pshort += wj * DE.tot[j];
    }
    rows.push(["Power absorbed", numHTML(DE.pabs, 3) + " erg s<sup>−1</sup>", "per H nucleus; equal to the power emitted, the integral of the infrared spectrum"]);
    rows.push(["IR power from <i>a</i> < 0.01 µm", (100 * psmall / pw).toFixed(1) + "%", "share of the emitted power from the smallest grains"]);
    rows.push(["IR power at λ < 50 µm", (100 * pshort / pw).toFixed(1) + "%", "share of the emitted power shortward of 50 um"]);
  } else {
    rows.push(["Power absorbed", "…", ""]);
  }
  $("readD").innerHTML = rows.map((r, k) =>
    `<div class="rd${k === 4 ? " gap" : ""}" title="${r[2]}"><span>${r[0]}</span><b>${r[1]}</b></div>`).join("");
}

// ================================================================== controls
function seg(ids, val, keys) { ids.forEach((id, k) => $(id).setAttribute("aria-pressed", String(keys[k] === val))); }
function buildDistControls() {
  for (const m of ["sil", "gra"]) {
    const el = $(m === "sil" ? "dSil" : "dGra"), key = m === "sil" ? "Si" : "C";
    el.innerHTML = `
      <label class="mh"><input type="checkbox" id="${m}On"><span class="sw" style="border-top-color:${C[m]}${m === "gra" ? ";border-top-style:dashed" : ""}"></span>${MATNAME[m]}</label>
      <div class="ctl" title="Slope of dn/da proportional to a^-q. MRN has q = 3.5: the number of grains rises steeply toward small sizes, the mass sits in the large grains for q &lt; 4."><label for="${m}Q"><span>Slope, <i>q</i></span><span class="val" id="${m}Qv"></span></label>
        <input type="range" id="${m}Q" min="${R.q[0]}" max="${R.q[1]}" step="0.05"></div>
      <div class="ctl" title="Smallest grain radius."><label for="${m}A0"><span>Smallest, <i>a</i><sub>min</sub></span><span class="val" id="${m}A0v"></span></label>
        <input type="range" id="${m}A0" min="${R.amin[0]}" max="${R.amin[1]}" step="0.005"></div>
      <div class="ctl" title="Largest grain radius."><label for="${m}A1"><span>Largest, <i>a</i><sub>max</sub></span><span class="val" id="${m}A1v"></span></label>
        <input type="range" id="${m}A1" min="${R.amax[0]}" max="${R.amax[1]}" step="0.005"></div>
      <div class="ctl" title="${key} atoms locked in these grains per million H nuclei.${m === "sil" ? " One Si per MgFeSiO4." : ""} For scale, the solar abundances are roughly 32 ppm of Si and 290 ppm of C."><label for="${m}X"><span>${key} in grains (ppm of H)</span><span class="val" id="${m}Xv"></span></label>
        <input type="range" id="${m}X" min="${R.ab[m][0]}" max="${R.ab[m][1]}" step="${m === "sil" ? 0.1 : 1}"></div>`;
    $(m + "On").addEventListener("change", ev => { S.dist[m].on = ev.target.checked; changed("dist"); });
    $(m + "Q").addEventListener("input", ev => { S.dist[m].q = +ev.target.value; changed("dist"); });
    $(m + "A0").addEventListener("input", ev => {
      const d = S.dist[m]; d.amin = +ev.target.value;
      if (d.amax < d.amin + 0.1) d.amax = Math.min(R.amax[1], d.amin + 0.1);
      changed("dist");
    });
    $(m + "A1").addEventListener("input", ev => {
      const d = S.dist[m]; d.amax = +ev.target.value;
      if (d.amin > d.amax - 0.1) d.amin = Math.max(R.amin[0], d.amax - 0.1);
      changed("dist");
    });
    $(m + "X").addEventListener("input", ev => { S.dist[m].abund = +ev.target.value; changed("dist"); });
  }
}
function syncControls() {
  document.body.className = "tab-" + S.tab;
  seg(["tGrain", "tDist"], S.tab, ["grain", "dist"]);
  seg(["mSil", "mGra"], S.mat, ["sil", "gra"]);
  seg(["sdMass", "sdArea"], S.sd, ["mass", "area"]);
  seg(["xH", "xV"], S.xn, ["H", "V"]);
  seg(["xsMat", "xsSize"], S.xs, ["mat", "size"]);
  seg(["spMat", "spSize"], S.sp, ["mat", "size"]);
  $("sU").value = S.U; $("vU").innerHTML = uHTML(S.U);
  $("sA").value = S.a; $("vA").innerHTML = aHTML(Math.pow(10, S.a));
  $("eLog").checked = S.elog; $("pLin").checked = S.plin;
  for (const m of ["sil", "gra"]) {
    const d = S.dist[m];
    $(m + "On").checked = d.on;
    $(m === "sil" ? "dSil" : "dGra").classList.toggle("off", !d.on);
    $(m + "Q").value = d.q; $(m + "Qv").textContent = d.q.toFixed(2);
    $(m + "A0").value = d.amin; $(m + "A0v").innerHTML = aHTML(Math.pow(10, d.amin));
    $(m + "A1").value = d.amax; $(m + "A1v").innerHTML = aHTML(Math.pow(10, d.amax));
    $(m + "X").value = d.abund; $(m + "Xv").textContent = d.abund.toFixed(m === "sil" ? 1 : 0);
  }
  const isMRN = ["sil", "gra"].every(m => { const d = S.dist[m], r = DEF.dist[m];
    return d.on && Math.abs(d.q - r.q) < 1e-9 && Math.abs(d.amin - r.amin) < 1e-9 && Math.abs(d.amax - r.amax) < 1e-9 && Math.abs(d.abund - r.abund) < 1e-9; });
  $("mrn").setAttribute("aria-pressed", String(isMRN));
}

// ================================================================== redraw orchestration
let raf = 0, dirty = new Set();
function changed(what) {
  syncControls();
  if (what === "grain" || what === "U" || what === "all") if (S.tab === "grain") { computeGrain(); dirty.add("grain"); }
  if (what === "dist" || what === "all") { if (S.tab === "dist") { computeDistOptics(); dirty.add("dist"); } }
  if (what === "U" && S.tab === "dist") dirty.add("dist");
  if (what === "view") dirty.add(S.tab);
  requestNodes();
  if (S.tab === "dist") computeDistEm();
  if (!raf) raf = requestAnimationFrame(flush);
  saveHash();
}
function flush() {
  raf = 0;
  if (S.tab === "grain") {
    if (!G1) computeGrain();
    drawQ(); drawW(); drawE(); drawP(); drawT(); readGrain();
    legend($("lgQ"), [{ color: C.ink, label: "<i>Q</i><sub>ext</sub>" }, { color: C.abs, label: "<i>Q</i><sub>abs</sub>" },
      { color: C.sca, dash: true, label: "<i>Q</i><sub>sca</sub>" }, { color: C.sca, label: "albedo" }, { color: C.ink, dot: true, label: "<i>g</i>" }]);
  } else {
    if (!DO) computeDistOptics();
    drawS(); drawX(); drawA(); drawI(); readDist();
  }
  dirty.clear();
}

// ================================================================== URL state
let hashTimer = 0;
function saveHash() {
  clearTimeout(hashTimer);
  hashTimer = setTimeout(() => {
    const d = m => { const x = S.dist[m]; return [x.on ? 1 : 0, x.q.toFixed(2), x.amin.toFixed(3), x.amax.toFixed(3), x.abund].join(","); };
    const h = `t=${S.tab}&m=${S.mat}&a=${S.a.toFixed(3)}&U=${S.U.toFixed(2)}&s=${d("sil")}&g=${d("gra")}` +
              `&el=${S.elog ? 1 : 0}&pl=${S.plin ? 1 : 0}&sd=${S.sd}&xn=${S.xn}&xs=${S.xs}&sp=${S.sp}`;
    try { history.replaceState(null, "", "#" + h); } catch (e) { /* sandboxed */ }
  }, 250);
}
function readHash() {
  let q;
  try { q = new URLSearchParams(location.hash.slice(1)); } catch (e) { return; }
  const num = (k, lo, hi) => { const v = parseFloat(q.get(k)); return isFinite(v) ? clamp(v, lo, hi) : null; };
  if (q.get("t") === "dist" || q.get("t") === "grain") S.tab = q.get("t");
  if (q.get("m") === "sil" || q.get("m") === "gra") S.mat = q.get("m");
  const a = num("a", R.a[0], R.a[1]); if (a !== null) S.a = a;
  const u = num("U", R.U[0], R.U[1]); if (u !== null) S.U = u;
  for (const [k, m] of [["s", "sil"], ["g", "gra"]]) {
    const v = (q.get(k) || "").split(",").map(parseFloat);
    if (v.length === 5 && v.every(isFinite)) {
      S.dist[m] = { on: v[0] > 0, q: clamp(v[1], R.q[0], R.q[1]), amin: clamp(v[2], R.amin[0], R.amin[1]),
                    amax: clamp(v[3], R.amax[0], R.amax[1]), abund: clamp(v[4], R.ab[m][0], R.ab[m][1]) };
    }
  }
  S.elog = q.get("el") === "1"; S.plin = q.get("pl") === "1";
  if (q.get("sd") === "area") S.sd = "area";
  if (q.get("xn") === "V") S.xn = "V";
  if (q.get("xs") === "size") S.xs = "size";
  if (q.get("sp") === "size") S.sp = "size";
}

// ================================================================== wiring
function wire() {
  $("tGrain").addEventListener("click", () => { S.tab = "grain"; changed("all"); });
  $("tDist").addEventListener("click", () => { S.tab = "dist"; changed("all"); });
  $("reset").addEventListener("click", () => { S = clone(DEF); HOV = null; changed("all"); });
  $("mrn").addEventListener("click", () => { S.dist = clone(DEF.dist); changed("dist"); });
  $("sU").addEventListener("input", ev => { S.U = +ev.target.value; changed("U"); });
  $("sA").addEventListener("input", ev => { S.a = +ev.target.value; changed("grain"); });
  $("mSil").addEventListener("click", () => { S.mat = "sil"; changed("grain"); });
  $("mGra").addEventListener("click", () => { S.mat = "gra"; changed("grain"); });
  $("eLog").addEventListener("change", ev => { S.elog = ev.target.checked; changed("view"); });
  $("pLin").addEventListener("change", ev => { S.plin = ev.target.checked; changed("view"); });
  $("sdMass").addEventListener("click", () => { S.sd = "mass"; changed("view"); });
  $("sdArea").addEventListener("click", () => { S.sd = "area"; changed("view"); });
  $("xH").addEventListener("click", () => { S.xn = "H"; changed("view"); });
  $("xV").addEventListener("click", () => { S.xn = "V"; changed("view"); });
  $("xsMat").addEventListener("click", () => { S.xs = "mat"; changed("view"); });
  $("xsSize").addEventListener("click", () => { S.xs = "size"; changed("view"); });
  $("spMat").addEventListener("click", () => { S.sp = "mat"; changed("view"); });
  $("spSize").addEventListener("click", () => { S.sp = "size"; changed("view"); });
  buildDistControls();

  $("cvQ").addEventListener("pointermove", ev => hoverQ(ev, AQ, $("cvQ")));
  $("cvW").addEventListener("pointermove", ev => hoverQ(ev, AW, $("cvW")));
  $("cvE").addEventListener("pointermove", hoverE);
  $("cvP").addEventListener("pointermove", hoverP);
  const cvT = $("cvT");
  cvT.addEventListener("pointerdown", ev => { tDrag = true; cvT.setPointerCapture(ev.pointerId); setAfromT(ev); });
  cvT.addEventListener("pointermove", ev => { if (tDrag) setAfromT(ev); hoverT(ev); });
  cvT.addEventListener("pointerup", () => { tDrag = false; });
  $("cvS").addEventListener("pointermove", hoverS);
  $("cvX").addEventListener("pointermove", ev => hoverX(ev, $("cvX"), AX));
  $("cvA").addEventListener("pointermove", ev => hoverX(ev, $("cvA"), AA));
  $("cvI").addEventListener("pointermove", hoverI);
  const out = grp => () => { hideTip(); if (HOV && HOV.grp === grp) { HOV = null; flush(); } };
  [["cvQ", "q"], ["cvW", "q"], ["cvE", "e"], ["cvP", "p"], ["cvS", "s"], ["cvX", "x"], ["cvA", "x"], ["cvI", "i"]]
    .forEach(([id, g]) => $(id).addEventListener("pointerleave", out(g)));
  $("cvT").addEventListener("pointerleave", hideTip);
  window.addEventListener("resize", () => { if (!raf) raf = requestAnimationFrame(flush); });
}

async function start() {
  const [meta, buf] = await Promise.all([fetch("data/dust.json").then(r => r.json()),
                                         fetch("data/dust.bin").then(r => r.arrayBuffer())]);
  P.init(meta, buf);
  LAMLOG = Float64Array.from(P.LAM, l => Math.log10(l));
  DEF = defaults();
  S = clone(DEF);
  readHash();
  wire();
  changed("all");
}
start();
})();
