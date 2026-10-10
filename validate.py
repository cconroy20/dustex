"""Checks of the Dust Explorer physics against its sources, and of the browser port against
the Python.  python3 validate.py   (about a minute; needs node for the port checks)

Reference values read from published figures were digitized from the 300 dpi renderings
of the arXiv figure files (rawdata/dl01, rawdata/ld01); each check says which.
"""
import gzip
import json
import os
import subprocess
import sys

import numpy as np

import dustlib as D

HERE = os.path.dirname(os.path.abspath(__file__))
results = []


def check(name, ok, detail=""):
    results.append((name, bool(ok)))
    print("%s  %-62s %s" % ("PASS" if ok else "FAIL", name, detail))


def readQ(path):
    """Draine's Q tables: {radius_um: rows of (lambda, Q_abs, Q_sca, g), lambda ascending}."""
    L = gzip.open(path, "rt").read().splitlines()
    nw = int(L[4].split()[0])
    out, i = {}, 5
    while i < len(L):
        if "radius" in L[i]:
            a = float(L[i].split()[0])
            rows = np.array([[float(v) for v in r.split()[:4]] for r in L[i + 2:i + 2 + nw]])
            out[a] = rows[::-1]
            i += 2 + nw
        else:
            i += 1
    return out


def bands_of(em, edges=((3, 8), (8, 20), (20, 50), (50, 150), (150, 1000))):
    return np.array([np.sum((D.WLN * D.LAM * em)[(D.LAM >= a) & (D.LAM < b)]) for a, b in edges])


# ------------------------------------------------------------------ Mie theory
q = D.mie(2 * np.pi * 0.525 / 0.6328, 1.55 + 0j)
check("BHMIE example (Bohren & Huffman App. A): Q_ext = Q_sca = 3.10543",
      abs(q[0][0] - 3.10543) < 5e-5 and abs(q[1][0] - 3.10543) < 5e-5, "%.5f %.5f" % (q[0][0], q[1][0]))
x = np.array([1e-3, 1e-2])
m = 1.7 + 0.03j
al = (m * m - 1) / (m * m + 2)
qe, qs, _ = D.mie(x, m)
check("dipole limits at small x (Q_abs = 4x Im alpha, Q_sca = 8/3 x^4 |alpha|^2)",
      np.all(abs((qe - qs) / (4 * x * al.imag) - 1) < 1e-3) and np.all(abs(qs / (8 / 3 * x ** 4 * abs(al) ** 2) - 1) < 1e-3))
check("extinction paradox: Q_ext(x = 300) -> 2", abs(D.mie(300.0, m)[0][0] - 2) < 0.06, "%.4f" % D.mie(300.0, m)[0][0])
try:
    import miepython
    worst = 0
    for mm in [1.7 + 0.03j, 1.5 + 0.001j, 3 + 4j, 30 + 40j, 400 + 400j]:
        for xx in [1e-3, 0.1, 1, 10, 100]:
            a = D.mie(xx, mm)
            b = miepython.efficiencies_mx(mm, xx)
            worst = max(worst, abs(a[0][0] / b[0] - 1), abs(a[1][0] / b[1] - 1))
    check("Q_ext, Q_sca vs miepython %s (25 cases, 1e-3 < x < 100)" % miepython.__version__, worst < 1e-5, "worst %.1e" % worst)
except ImportError:
    print("SKIP  miepython not installed")

# ------------------------------------------------------------------ optics vs Draine's own tables
for mat, fn, tol_bands in [("sil", "suvSil_81.gz", [(0.0912, 1000, 0.05)]),
                           ("gra", "Gra_81.gz", [(0.0912, 30, 0.03)])]:
    T = readQ(os.path.join(D.RAW, "draine", fn))
    worst = 0
    for a in [0.001, 0.01, 0.1, 1.0]:
        k = min(T, key=lambda r: abs(np.log(r / a)))
        w, qa, qsc, g = T[k].T
        for lo, hi, tol in tol_bands:
            sel = (w >= lo) & (w <= hi)
            Qa, Qs, G = D.optics(mat, k, w[sel])
            worst = max(worst, abs(np.median(Qa / qa[sel]) - 1))
    check("%s Q_abs vs Draine %s, median ratio over %g-%g um, 4 radii" % (mat, fn, lo, hi), worst < tol, "worst %.3f" % worst)
T = readQ(os.path.join(D.RAW, "draine", "Gra_81.gz"))
w, qa, _, _ = T[min(T, key=lambda r: abs(np.log(r / 0.1)))].T
sel = (w >= 100) & (w <= 1000)
r = np.median(D.optics("gra", 0.1, w[sel])[0] / qa[sel])
print("INFO  graphite 0.1 um, 100-1000 um: Q_abs(D03) / Q_abs(Gra_81, 1993 constants) = %.2f" % r)

# ------------------------------------------------------------------ radiation field
tl, tj = D._mmp83_table()
ratio = D.J_bb(tl) / tj
long = tl >= 0.3
check("MMP83 Table A3 vs Table A1 blackbodies, 0.35-8 um (incl. corrected 3.4 um)",
      np.all(abs(ratio[long] - 1) < 0.015), "range %.3f-%.3f" % (ratio[long].min(), ratio[long].max()))
res = tj - D.J_bb(tl)
i21 = int(np.argmin(abs(tl - 0.21)))
check("MMP83 UV component smooth through the corrected 0.21 um entry",
      res[i21 - 1] > res[i21] > res[i21 + 1], "%.2e > %.2e > %.2e" % (res[i21 - 1], res[i21], res[i21 + 1]))
lg = np.logspace(np.log10(0.0912), np.log10(8), 20000)
I = np.trapz(D.J_uv(lg) + D.J_bb(lg), lg)
check("int 4 pi J_lambda, 0.0912-8 um = 2.17e-2 (MMP83 Eq. A4)", abs(I / 2.17e-2 - 1) < 0.02, "%.4e" % I)
ucmb = D.integ_lam(D.u_cmb())
check("CMB energy density on the grid = a T^4", abs(ucmb / (7.5657e-15 * D.T_CMB ** 4) - 1) < 0.01,
      "%.4e" % ucmb)

# ------------------------------------------------------------------ heat capacities
fig2_sil = {10: 0.01158, 12: 0.01666, 15: 0.0260, 20: 0.04639, 25: 0.07268, 100: 0.9199, 200: 1.735, 500: 2.581}
c = D.heatcap_per_atom("sil", list(fig2_sil))
dev = np.max(abs(c / np.array(list(fig2_sil.values())) - 1))
check("silicate C/Nk vs DL01 Fig. 2 model curve (8 temperatures)", dev < 0.01, "max dev %.4f" % dev)
lab = {20: 0.0099, 30: 0.0223}           # deSorbo & Tyler graphite, as plotted in DL01 Fig. 2
c = D.heatcap_per_atom("gra", list(lab))
dev = np.max(abs(c / np.array(list(lab.values())) - 1))
check("graphite C/Nk vs the lab points in DL01 Fig. 2 (20, 30 K)", dev < 0.06, "max dev %.3f" % dev)
for mat in ("sil", "gra"):
    cl = D.heatcap_per_atom(mat, [4900.0])[0]
    check("%s C/Nk -> 3 at high T (equipartition; needs f_n prefactor n)" % mat, cl > 2.9, "%.3f at 4900 K" % cl)

# ------------------------------------------------------------------ equilibrium temperatures
# LD01 Fig. 3 (chi_MMP = 1) was computed with Draine's tabulated cross sections; with those
# tables the solver should reproduce it.  T in K, digitized.
fig3 = {("sil", 0.0102): 17.79, ("sil", 0.0486): 16.63, ("sil", 0.1948): 14.96, ("sil", 0.4868): 14.30,
        ("gra", 0.0105): 20.44, ("gra", 0.0486): 20.24, ("gra", 0.1948): 18.26, ("gra", 0.4868): 15.38}
tabs93 = {"sil": readQ(os.path.join(D.RAW, "draine", "suvSil_81.gz")),
          "gra": readQ(os.path.join(D.RAW, "draine", "Gra_81.gz"))}
worst = 0
for (mat, a), Tf in fig3.items():
    T = tabs93[mat]
    k = min(T, key=lambda r: abs(np.log(r / a)))
    w, qa, _, _ = T[k].T
    Qa = np.exp(np.interp(np.log(D.LAM), np.log(w), np.log(qa)))
    Qa = np.where(D.LAM > w.max(), qa[-1] * (w.max() / D.LAM) ** 2, Qa)
    g = D.Grain(mat, k, 1.0, opt=(Qa, 0 * Qa, 0 * Qa))
    worst = max(worst, abs(g.Teq - Tf))
check("T_eq with Draine's tables vs LD01 Fig. 3 (8 grains, U = 1)", worst < 0.35, "max |dT| %.2f K" % worst)

# ------------------------------------------------------------------ stochastic heating
worst = 0
for mat in ("sil", "gra"):
    for a in (0.003, 0.01, 0.03):
        for U in (1.0, 1e3):
            g = D.Grain(mat, a, U)
            T, P, info = g.stochastic()
            worst = max(worst, abs(info["Pemit"] / g.Pabs - 1))
check("stochastic steady state conserves energy (12 grains)", worst < 1e-3, "worst %.1e" % worst)
# LD01 Fig. 5: silicate, chi_MMP = 1; peak of dP/dlnE and its position E/hc (cm^-1)
fig5 = {0.005: (0.397, 767), 0.0075: (0.483, 6940), 0.010: (0.594, 2.80e4), 0.015: (1.154, 1.50e5),
        0.020: (1.863, 3.72e5), 0.025: (2.61, 7.20e5), 0.030: (3.354, 1.22e6)}
hd, pd = [], []
for a, (hf, ef) in fig5.items():
    g = D.Grain("sil", a, 1.0)
    T, P, info = g.stochastic()
    E = info["E"]
    dens = P / np.gradient(np.log(E))
    k = int(np.argmax(dens))
    hd.append(dens[k] / hf - 1)
    pd.append(E[k] / (D.H * D.C) / ef - 1)
check("P(E) peak heights vs LD01 Fig. 5 (silicate, 50-300 A)", max(abs(v) for v in hd) < 0.12,
      "dev " + " ".join("%+.2f" % v for v in hd))
check("P(E) peak positions vs LD01 Fig. 5 (see README: ~12% offset)", max(abs(v) for v in pd) < 0.2,
      "dev " + " ".join("%+.2f" % v for v in pd))
# convergence: default binning vs a much finer reference, bands holding > 0.1% of the power
worst = 0
for mat in ("sil", "gra"):
    for a in (0.003, 0.005, 0.01, 0.02):
        for U in (1.0, 100.0):
            g = D.Grain(mat, a, U)
            if g.sigma_T()[0] < g.SIGMA_EQ:
                continue
            e1 = bands_of(g.emission()[0])
            g.M, g.DELTA, g.M_MAX = 1000, 0.008, 12000
            e2 = bands_of(g.emission()[0])
            f = e2 / e2.sum()
            worst = max(worst, np.max(abs(e1 / e2 - 1)[f > 1e-3]))
check("emission converged in the bin number (< 3% in bands with > 0.1% of power)", worst < 0.03, "worst %.4f" % worst)
# the equilibrium switch: just below SIGMA_EQ, equilibrium and stochastic spectra agree
worst = 0
for mat, a, U in [("sil", 0.08, 1.0), ("gra", 0.07, 1.0), ("sil", 0.03, 100.0), ("sil", 0.12, 0.1)]:
    g = D.Grain(mat, a, U)
    s = g.sigma_T()[0]
    g.SIGMA_EQ = 0
    e1 = bands_of(g.emission()[0])
    e2 = bands_of(4 * np.pi * g.Cabs * D.B_lam(D.LAM, g.Teq))
    f = e1 / e1.sum()
    worst = max(worst, np.max(abs(e2 / e1 - 1)[f > 1e-3]))
check("equilibrium within 1% of stochastic near the switch (sigma_T/T ~ 0.01)", worst < 0.01, "worst %.4f" % worst)

# ------------------------------------------------------------------ size distributions
d = np.load(os.path.join(HERE, "cache", "optics.npz"))
tabs32 = {m: d[m].astype("f4").astype("f8") for m in ("sil", "gra")}
mrn = {m: dict(q=3.5, amin=0.005, amax=0.25, abund=D.mrn_abundance(m)) for m in ("sil", "gra")}
r = D.distribution(mrn, 1.0, tabs32)
pem = 4 * np.pi * np.sum(D.WLN * sum(r["em"].values()))
pab = sum(r["pabs"].values())
check("size distribution: emitted = absorbed power per H (MRN, U = 1)", abs(pem / pab - 1) < 1e-3, "%.5f" % (pem / pab))
check("MRN (DL84 normalization) locks up 34 ppm Si and 298 ppm C",
      abs(D.mrn_abundance("sil") - 34.16) < 0.05 and abs(D.mrn_abundance("gra") - 297.9) < 0.3)
bins = sum(bands_of(e) for e in r["em_size"]).sum()
check("emission split by size sums to the total", abs(bins / bands_of(sum(r["em"].values())).sum() - 1) < 1e-9)

# ------------------------------------------------------------------ browser port
cases = []
for mat in ("sil", "gra"):
    cases.append(dict(kind="optics", mat=mat, a=0.0137))
    cases.append(dict(kind="optics", mat=mat, a=0.42))
for mat, a, U in [("sil", 0.003, 1.0), ("gra", 0.01, 1.0), ("sil", 0.05, 1.0), ("gra", 0.3, 1.0),
                  ("sil", 0.02, 1000.0), ("gra", 0.004, 0.1)]:
    cases.append(dict(kind="grain", mat=mat, a=a, U=U))
cases.append(dict(kind="grain", mat="sil", a=float(D.A_GRID[40]), U=3.0, table=40))
dist = {"sil": dict(q=3.3, amin=0.004, amax=0.3, abund=30.0), "gra": dict(q=3.6, amin=0.006, amax=0.2, abund=200.0)}
cases.append(dict(kind="dist", U=2.0, params=dist))
cpath = os.path.join(HERE, "cache", "cases.json")
json.dump(cases, open(cpath, "w"))
try:
    js = json.loads(subprocess.run(["node", os.path.join(HERE, "testjs.js"), cpath], check=True,
                                   capture_output=True, text=True).stdout)
except (OSError, subprocess.CalledProcessError) as e:
    js = None
    check("browser port ran under node", False, str(e))
if js is not None:
    w_opt = w_g = w_d = 0
    for c, j in zip(cases, js):
        if c["kind"] == "optics":
            Qa, Qs, G = D.optics(c["mat"], c["a"])
            w_opt = max(w_opt, np.max(abs(np.array(j["Qabs"]) / Qa - 1)), np.max(abs(np.array(j["Qsca"]) / Qs - 1)),
                        np.max(abs(np.array(j["g"]) - G)))
        elif c["kind"] == "grain":
            opt = tuple(tabs32[c["mat"]][:, c["table"]]) if "table" in c else None
            g = D.Grain(c["mat"], c["a"], c["U"], opt=opt)
            em = g.emission()[0]
            b1, b2 = bands_of(em), bands_of(np.array(j["em"]))
            f = b1 / b1.sum()
            w_g = max(w_g, abs(j["Teq"] / g.Teq - 1), abs(j["Pabs"] / g.Pabs - 1), np.max(abs(b2 / b1 - 1)[f > 1e-3]))
        else:
            rp = D.distribution(c["params"], c["U"], tabs32)
            for mat in c["params"]:
                w_d = max(w_d, np.max(abs(np.array(j["ext"][mat]) / rp["ext"][mat] - 1)),
                          np.max(abs(np.array(j["sca"][mat]) / rp["sca"][mat] - 1)), abs(j["mass"][mat] / rp["mass"][mat] - 1))
                b1, b2 = bands_of(rp["em"][mat]), bands_of(np.array(j["em"][mat]))
                f = b1 / b1.sum()
                w_d = max(w_d, np.max(abs(b2 / b1 - 1)[f > 1e-3]))
    check("port: Mie optics (2 materials x 2 radii, all wavelengths)", w_opt < 1e-9, "worst %.1e" % w_opt)
    check("port: single grains (T_eq, P_abs, band emission; 7 cases)", w_g < 1e-4, "worst %.1e" % w_g)
    check("port: size distribution (extinction, scattering, mass, emission)", w_d < 1e-4, "worst %.1e" % w_d)

n = sum(ok for _, ok in results)
print("\n%d passed, %d failed" % (n, len(results) - n))
sys.exit(0 if n == len(results) else 1)
