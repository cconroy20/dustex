"""Write the data bundle the page loads: web/data/dust.json (layout and constants) and
web/data/dust.bin (the arrays).

  lam            wavelength grid (um), float64
  ustar          MMP83 starlight energy density for U = 1 (erg cm^-3 um^-1), float64
  m_sil          astronomical silicate refractive index on lam (re, im interleaved), float64
  eps_{pa,pe}{0.01,0.10}  graphite dielectric functions on lam (re, im interleaved), float64
  q_{sil,gra}    Q_abs then Q_sca on (A_GRID, lam), float32
  lamf, ustarf   the fine wavelength grid (PAHs and distribution sums) and the field on it
  pah_grid       PAH size nodes (um)
  tgrid_e        temperature grid for the enthalpy tables (K), float64
  e_{sil,gra}    E(T) / (N - 2) on tgrid_e (erg), float64

The page runs the Mie series itself for the single grain; the Q tables serve the size
distributions.  Run after any change to dustlib's grids or data:  python3 exportweb.py
"""
import json
import os

import numpy as np

import dustlib as D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "web", "data")


def main():
    os.makedirs(OUT, exist_ok=True)
    cache = os.path.join(HERE, "cache", "optics.npz")
    if os.path.exists(cache):
        d = np.load(cache)
        ok = d["lam"].shape == D.LAM.shape and np.allclose(d["lam"], D.LAM) and np.allclose(d["a"], D.A_GRID)
    else:
        ok = False
    if ok:
        tabs = {m: d[m] for m in ("sil", "gra")}
    else:
        tabs = {m: D.optics_table(m) for m in ("sil", "gra")}
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        np.savez(cache, sil=tabs["sil"], gra=tabs["gra"], lam=D.LAM, a=D.A_GRID)
    E = D._load_eps()

    def cplx(z):
        return np.column_stack([z.real, z.imag]).ravel()

    arrays = [
        ("lam", D.LAM, "f8"),
        ("ustar", D.u_star(), "f8"),
        ("m_sil", cplx(np.sqrt(E["sil"])), "f8"),
        ("eps_pa0.01", cplx(E["pa0.01"]), "f8"),
        ("eps_pa0.10", cplx(E["pa0.10"]), "f8"),
        ("eps_pe0.01", cplx(E["pe0.01"]), "f8"),
        ("eps_pe0.10", cplx(E["pe0.10"]), "f8"),
        ("q_sil", tabs["sil"][:2].astype("f4"), "f4"),
        ("q_gra", tabs["gra"][:2].astype("f4"), "f4"),
        ("lamf", D.LAMF, "f8"),
        ("ustarf", D.u_star(D.LAMF), "f8"),
        ("pah_grid", D.PAH_GRID, "f8"),
        ("tgrid_e", D.TGRID_E, "f8"),
        ("e_sil", D.enthalpy_table("sil"), "f8"),
        ("e_gra", D.enthalpy_table("gra"), "f8"),
    ]
    layout, blobs, off = [], [], 0
    for name, arr, dt in arrays:
        b = np.ascontiguousarray(arr, dtype=dt).tobytes()
        pad = (-off) % 8
        blobs.append(b"\0" * pad)
        off += pad
        layout.append(dict(name=name, dtype=dt, offset=off, length=int(np.asarray(arr).size),
                           shape=list(np.shape(arr))))
        blobs.append(b)
        off += len(b)
    meta = dict(
        arrays=layout,
        a_grid=dict(lo=-3.0, hi=1.0, n=len(D.A_GRID)),
        nlam=D.NLAM,
        T_CMB=D.T_CMB,
        mat={m: dict(rho=D.MAT[m]["rho"], atoms_per_cm3=D.MAT[m]["atoms_per_cm3"],
                     per_key=7.0 if m == "sil" else 1.0) for m in ("sil", "gra")},
        mrn={m: D.mrn_abundance(m) for m in ("sil", "gra")},
    )
    with open(os.path.join(OUT, "dust.bin"), "wb") as f:
        f.write(b"".join(blobs))
    with open(os.path.join(OUT, "dust.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print("wrote %s (%.0f kB)" % (os.path.join(OUT, "dust.bin"), off / 1e3))


if __name__ == "__main__":
    main()
