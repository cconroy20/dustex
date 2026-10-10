// Run the browser physics under node for validate.py:  node testjs.js cases.json > out.json
// Each case is {kind: "mie" | "optics" | "pahopt" | "grain" | "dist" | "ccm", ...}; see validate.py.
"use strict";
const fs = require("fs"), path = require("path");
const Phys = require("./web/physics.js");
const meta = JSON.parse(fs.readFileSync(path.join(__dirname, "web/data/dust.json")));
const b = fs.readFileSync(path.join(__dirname, "web/data/dust.bin"));
Phys.init(meta, b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength));
const arr = a => Array.from(a);
const add = (acc, xs) => acc ? acc.map((v, k) => v.map((y, j) => y + xs[k][j])) : xs.map(arr);
const cases = JSON.parse(fs.readFileSync(process.argv[2]));
const out = cases.map(c => {
  const t0 = Date.now();
  if (c.kind === "ccm") return { y: c.x.map(x => Phys.ccm89(x, c.Rv)) };
  if (c.kind === "mie") return { q: Phys.mie(c.x, c.mr, c.mi) };
  if (c.kind === "optics") { const o = Phys.optics(c.mat, c.a); return { Qabs: arr(o.Qabs), Qsca: arr(o.Qsca), g: arr(o.g) }; }
  if (c.kind === "pahopt") { const o = Phys.pahOptics(c.a, c.charge); return { Qabs: arr(o.Qabs), Qsca: arr(o.Qsca), g: arr(o.g) }; }
  if (c.kind === "grain") {
    const opt = c.table !== undefined ? Phys.tableOptics(c.mat, c.table) : null;
    const g = new Phys.Grain(c.mat, c.a, c.U, opt, c.charge);
    const r = g.emission();
    return { Teq: g.Teq, Pabs: g.Pabs, sigmaT: r.sigmaT, stochastic: r.stochastic, em: arr(r.em),
             T: arr(r.T), P: arr(r.P), ms: Date.now() - t0 };
  }
  if (c.kind === "dist") {
    const res = { ext: {}, sca: {}, em: {}, mass: {}, pabs: {}, bins: null, extB: null, scaB: null, ms: 0 };
    for (const mat of Object.keys(c.params)) {
      const p = c.params[mat];
      let o, e;
      if (mat === "pah") {
        o = Phys.distOpticsPAH(p);
        const [i0, i1] = Phys.pahNodeRange(p.a0, p.sigma), nodes = { neu: [], ion: [] };
        for (let i = i0; i <= i1; i++) for (const ch of ["neu", "ion"]) nodes[ch][i] = Phys.nodeEmission("pah", i, c.U, ch);
        e = Phys.distEmissionPAH(p, nodes);
      } else {
        o = Phys.distOptics(mat, p);
        const [i0, i1] = Phys.nodeRange(p.amin, p.amax), nodes = [];
        for (let i = i0; i <= i1; i++) nodes[i] = Phys.nodeEmission(mat, i, c.U);
        e = Phys.distEmission(mat, p, nodes);
      }
      res.ext[mat] = arr(o.ext); res.sca[mat] = arr(o.sca); res.mass[mat] = o.mass;
      res.em[mat] = arr(e.em); res.pabs[mat] = e.pabs;
      res.bins = add(res.bins, e.bins); res.extB = add(res.extB, o.extB); res.scaB = add(res.scaB, o.scaB);
    }
    res.ms = Date.now() - t0;
    return res;
  }
  throw new Error("unknown case " + c.kind);
});
process.stdout.write(JSON.stringify(out));
