// Drive the page in headless Chrome and check that each interaction does what it should:
// sliders, tabs, toggles, hover tooltips on every panel, clicking the temperature panel,
// the background computation of the per-size emission, and the phone-width layout.
// Reading the code does not catch a dead handler; driving it does.
//   python3 web/serve.py 8781 &   then   node uitest.js [url] [outdir]
"use strict";
const { spawn } = require("child_process");
const fs = require("fs");
const path = require("path");

const URL0 = process.argv[2] || "http://localhost:8781/";
const OUT = process.argv[3] || path.join(__dirname, "cache", "shots");
const CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const PORT = 9335;
const sleep = ms => new Promise(r => setTimeout(r, ms));
let fails = 0;
function check(name, ok, detail) {
  if (!ok) fails++;
  console.log(`[${ok ? "PASS" : "FAIL"}] ${name}${detail !== undefined ? "  --  " + detail : ""}`);
}

(async () => {
  fs.mkdirSync(OUT, { recursive: true });
  const prof = fs.mkdtempSync("/tmp/dustex-chrome-");
  const chrome = spawn(CHROME, ["--headless=new", `--remote-debugging-port=${PORT}`,
    `--user-data-dir=${prof}`, "--window-size=1440,1500", "--hide-scrollbars", "about:blank"], { stdio: "ignore" });
  let tabs = null;
  for (let i = 0; i < 50 && !tabs; i++) {
    await sleep(200);
    try { tabs = await (await fetch(`http://127.0.0.1:${PORT}/json`)).json(); } catch (e) { /* not up yet */ }
  }
  const page = tabs.find(t => t.type === "page");
  const ws = new WebSocket(page.webSocketDebuggerUrl);
  await new Promise(r => ws.addEventListener("open", r));
  let id = 0; const pend = new Map(); const errors = [];
  ws.addEventListener("message", ev => {
    const m = JSON.parse(ev.data);
    if (m.id && pend.has(m.id)) { pend.get(m.id)(m); pend.delete(m.id); }
    if (m.method === "Runtime.exceptionThrown") errors.push(m.params.exceptionDetails.exception?.description || m.params.exceptionDetails.text);
    if (m.method === "Runtime.consoleAPICalled" && m.params.type === "error") errors.push(JSON.stringify(m.params.args));
    if (m.method === "Log.entryAdded" && m.params.entry.level === "error") errors.push(m.params.entry.text + " " + (m.params.entry.url || ""));
  });
  const send = (method, params = {}) => new Promise(r => { const i = ++id; pend.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
  const ev = async expr => (await send("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true })).result.result.value;
  const shot = async name => {
    const r = await send("Page.captureScreenshot", { format: "png", captureBeyondViewport: true });
    fs.writeFileSync(path.join(OUT, name), Buffer.from(r.result.data, "base64"));
  };
  const mouse = (type, x, y) => send("Input.dispatchMouseEvent", { type, x, y, button: "left", buttons: type === "mouseReleased" ? 0 : 1, clickCount: 1 });
  const move = (x, y) => send("Input.dispatchMouseEvent", { type: "mouseMoved", x, y });
  const rect = async sel => ev(`(r => ({x: r.left, y: r.top, w: r.width, h: r.height}))(document.querySelector(${JSON.stringify(sel)}).getBoundingClientRect())`);
  const scrollTo = async sel => { await ev(`document.querySelector(${JSON.stringify(sel)}).scrollIntoView({block: "center"})`); await sleep(80); };
  const setRange = async (sel, v) => { await ev(`(e => { e.value = ${v}; e.dispatchEvent(new Event("input", {bubbles: true})); })(document.querySelector(${JSON.stringify(sel)}))`); await sleep(60); };
  const clickSel = async sel => { await ev(`document.querySelector(${JSON.stringify(sel)}).click()`); await sleep(80); };
  const idle = async (ms = 15000) => {      // wait until the node queue has drained
    const t0 = Date.now();
    while (Date.now() - t0 < ms) {
      const busy = await ev(`document.getElementById("tSub").textContent + document.getElementById("iSub").textContent`);
      if (!busy) return Date.now() - t0;
      await sleep(100);
    }
    return -1;
  };
  // sweep the pointer across a canvas and count the distinct tooltip contents
  const sweep = async sel => {
    await scrollTo(sel);
    const r = await rect(sel), seen = new Set();
    for (let k = 1; k < 24; k++) {
      await move(r.x + r.w * (0.08 + 0.86 * k / 24), r.y + r.h * 0.5);
      await sleep(25);
      const t = await ev(`document.getElementById("tip").hidden ? "" : document.getElementById("tip").textContent`);
      if (t) seen.add(t);
    }
    await move(5, 5);
    return seen.size;
  };
  await send("Runtime.enable"); await send("Log.enable"); await send("Page.enable");
  await send("Page.navigate", { url: URL0 });
  await sleep(1500);

  // ---------------------------------------------------------------- single grain
  const t1 = await idle();
  check("page loads, single-grain nodes finish computing", t1 >= 0, `${t1} ms`);
  check("no errors on load", errors.length === 0, errors.join(" | "));
  const read = async () => ev(`[...document.querySelectorAll("#readG .rd")].map(r => r.textContent).join(" | ")`);
  const r0 = await read();
  check("readouts filled", /Equilibrium/.test(r0) && /K/.test(r0), r0.slice(0, 120));
  await shot("grain_default.png");
  for (const id of ["#cvQ", "#cvW", "#cvE", "#cvP", "#cvT"]) {
    const n = await sweep(id);
    check(`hover tooltip changes across ${id}`, n >= (id === "#cvP" ? 0 : 8), `${n} distinct`);
  }
  // radius slider: a small grain becomes stochastic
  await setRange("#sA", Math.log10(0.004));
  const r1 = await read();
  check("small grain: temperature spread reported", /Middle 90% of T[\d.]+ K – [\d.]+ K/.test(r1.replace(/\s+/g, " ")), r1.split(" | ")[4]);
  check("small grain: stochastic legend", /P\(T\)/.test(await ev(`document.getElementById("lgP").textContent`)));
  const np = await sweep("#cvP");
  check("hover tooltip on P(T) for a stochastic grain", np >= 6, `${np} distinct`);
  await shot("grain_small.png");
  // graphite, larger field, log energy axis
  await clickSel("#mGra");
  await setRange("#sU", 3);
  await ev(`(e => { e.checked = true; e.dispatchEvent(new Event("change")); })(document.getElementById("eLog"))`);
  const t2 = await idle();
  check("graphite at U = 1000 computes", t2 >= 0, `${t2} ms`);
  await shot("grain_gra_U1000.png");
  // clicking the temperature panel sets the radius
  await scrollTo("#cvT");
  const rt = await rect("#cvT");
  const before = await ev(`document.getElementById("sA").value`);
  await mouse("mousePressed", rt.x + rt.w * 0.75, rt.y + rt.h * 0.5);
  await mouse("mouseReleased", rt.x + rt.w * 0.75, rt.y + rt.h * 0.5);
  await sleep(150);
  const after = await ev(`document.getElementById("sA").value`);
  check("click on the temperature panel moves the radius", Math.abs(after - before) > 0.2, `${before} -> ${after}`);
  check("hash records state", /t=grain&m=gra/.test(await ev("location.hash")), await ev("location.hash"));

  // ---------------------------------------------------------------- distribution
  await clickSel("#reset");
  await clickSel("#tDist");
  const t3 = await idle();
  check("distribution tab computes its emission", t3 >= 0, `${t3} ms`);
  const rd = await ev(`[...document.querySelectorAll("#readD .rd")].map(r => r.textContent).join(" | ")`);
  check("distribution readouts filled", /R/.test(rd) && /%/.test(rd), rd);
  await shot("dist_mrn.png");
  for (const id of ["#cvS", "#cvX", "#cvA", "#cvI"]) {
    const n = await sweep(id);
    check(`hover tooltip changes across ${id}`, n >= 8, `${n} distinct`);
  }
  await setRange("#silQ", 3.0);
  await setRange("#graA1", Math.log10(1.0));
  await clickSel("#spSize"); await clickSel("#xV"); await clickSel("#sdArea");
  await clickSel("#xsSize");
  const lgx = await ev(`document.getElementById("lgX").textContent`);
  check("extinction split by size: legend lists the size bins", /0\.01–0\.1/.test(lgx) && !/Silicate/.test(lgx), lgx);
  const nx = await sweep("#cvX");
  check("hover tooltip on the size-split extinction", nx >= 8, `${nx} distinct`);
  await setRange("#sU", 2);
  const t4 = await idle();
  check("distribution recomputes after U and slope changes", t4 >= 0, `${t4} ms`);
  const rd2 = await ev(`[...document.querySelectorAll("#readD .rd")].map(r => r.textContent).join(" | ")`);
  check("readouts change with the distribution", rd2 !== rd);
  await shot("dist_changed.png");
  await ev(`(e => { e.checked = false; e.dispatchEvent(new Event("change")); })(document.getElementById("graOn"))`);
  await idle(); await sleep(150);
  check("graphite off: legend lists silicate only", !/Graphite/.test(await ev(`document.getElementById("lgX").textContent`)));
  await clickSel("#mrn");
  check("MRN button restores the preset", (await ev(`document.getElementById("mrn").getAttribute("aria-pressed")`)) === "true");

  // ---------------------------------------------------------------- reload from the hash
  await sleep(500);
  const h = await ev("location.hash");
  await send("Page.navigate", { url: "about:blank" }); await sleep(300);
  await send("Page.navigate", { url: URL0.replace(/#.*/, "") + h });
  await sleep(1200); await idle(); await sleep(500);
  check("state survives a reload from the hash", (await ev("location.hash")) === h, h);

  // ---------------------------------------------------------------- phone width
  await send("Emulation.setDeviceMetricsOverride", { width: 390, height: 900, deviceScaleFactor: 2, mobile: true });
  await sleep(400);
  const sw = await ev("document.documentElement.scrollWidth");
  check("no horizontal scroll at 390 px", sw <= 391, `scrollWidth ${sw}`);
  await shot("phone.png");
  await send("Emulation.clearDeviceMetricsOverride");

  check("no errors at the end", errors.length === 0, errors.join(" | "));
  console.log(fails ? `${fails} FAILED` : "all passed");
  ws.close(); chrome.kill();
  process.exit(fails ? 1 : 0);
})();
