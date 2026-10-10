"""dustlib: the physics of the Dust Explorer.

This is the reference implementation.  web/physics.js is a port of the parts the page
computes live (Mie series, radiation field, equilibrium and stochastic temperatures,
emission, size-distribution sums), and validate.py checks the two against each other.

Grains are homogeneous spheres of two materials:
  astronomical silicate  Draine (2003) dielectric function (rawdata/draine/callindex.out_silD03),
                         rho = 3.5 g cm^-3, MgFeSiO4 (Li & Draine 2001, Sect. 3), so
                         N_atom = rho V / (172.23/7 amu)
  graphite               Draine (2003) dielectric functions for E parallel and perpendicular to
                         the c axis, free-electron damping for a = 0.01 and 0.10 um
                         (callindex.out_C{pa,pe}D03_{0.01,0.10}); rho = 2.24 g cm^-3;
                         N_C = 468 (a / 1 nm)^3 (Draine's dust.diel page); randomly oriented
                         spheres in the "1/3-2/3" approximation (Draine & Malhotra 1993):
                         Q = Q(eps_par)/3 + 2 Q(eps_perp)/3
Between a = 0.01 and 0.10 um the graphite dielectric function is interpolated linearly in
log a; outside that range the nearer file is used.

Optics: Mie theory, the BHMIE algorithm of Bohren & Huffman (1983, Appendix A),
vectorized over wavelength.

Radiation field (Li & Draine 2001, Eq. 6 form):
  u_lambda = U [u_UV + sum_i W_i (4 pi / c) B_lambda(T_i)] + (4 pi / c) B_lambda(T_CMB)
with (T_i, W_i) = (7500 K, 1e-14), (4000 K, 1e-13), (3000 K, 4e-13) from MMP83 Table A1 and
the UV component taken as MMP83 Table A3 (D_G = 10 kpc) minus the three blackbodies,
which is how MMP83 define it.  T_CMB = 2.7255 K (Fixsen 2009).

Heat capacities: Draine & Li (2001), Eqs. (eq:E_PAH_Debye) with N_H = 0 for graphite and
(eq:E_sil_Debye) for silicate,
  E_gra(T) = (N_C - 2) [k Theta_op f_2(T/863 K) + 2 k Theta_ip f_2(T/2504 K)]
  E_sil(T) = (N_a - 2) [2 k Theta_2 f_2(T/500 K) + k Theta_3 f_3(T/1500 K)]
  f_n(x)   = n int_0^1 y^n dy / [exp(y/x) - 1].
⚠ DL01 (and LD01) print f_n with a prefactor 1/n.  The n-dimensional Debye model with the
mode counts DL01 give (N_m = N_C - 2 out of plane, 2(N_C - 2) in plane) requires n, and only
n gives equipartition, E -> 3 (N - 2) k T at high T, and reproduces the model curves in their
Fig. 2 (graphite C/Nk = 1.0 at 300 K).  validate.py checks both.

Stochastic heating: the "thermal-continuous" approximation of DL01 Sect. 7.2, solved by the
Guhathakurta & Draine (1989) recursion (DL01 Sect. 9.2), with the finite-bin-width upward
rates of DL01 Eqs. (eq:T_ul)-(eq:Tu,u-1) and the thermal emission of Eq. (eq:nuFnu_thermal).
Grains large enough that the temperature fluctuations are below SIGMA_EQ of T_eq are
treated as being at T_eq.

Units: wavelength in um; u_lambda in erg cm^-3 um^-1; powers in erg s^-1; cross sections in
cm^2; energies in erg unless a name says eV.
"""
import gzip
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "rawdata")
DATA = os.path.join(HERE, "data")

# ------------------------------------------------------------------ constants (cgs, CODATA 2018)
H = 6.62607015e-27
C = 2.99792458e10
KB = 1.380649e-16
EV = 1.602176634e-12
AMU = 1.66053906660e-24
UM = 1e-4
T_CMB = 2.7255                       # Fixsen (2009)
LYMAN_UM = 0.0912                    # 13.6 eV; no starlight shortward of this

# ------------------------------------------------------------------ grids
NLAM = 404                           # 80 per decade from the Lyman limit to 1 cm
LAM = 10 ** np.linspace(np.log10(LYMAN_UM), 4.0, NLAM)
DLNLAM = np.log(LAM[1] / LAM[0])
A_GRID = 10 ** np.linspace(-3.0, 1.0, 101)    # 0.001 - 10 um, 25 per decade

# trapezoid weights in ln(lambda): int F dlambda = sum(WLN * LAM * F)
WLN = np.full(NLAM, DLNLAM)
WLN[0] = WLN[-1] = 0.5 * DLNLAM


def integ_lam(F):
    """int F dlambda over the wavelength grid (F per um)."""
    return np.sum(WLN * LAM * F, axis=-1)


# ------------------------------------------------------------------ materials
MAT = {
    "sil": dict(name="astronomical silicate", rho=3.5,
                atoms_per_cm3=3.5 / (172.232 / 7 * AMU),
                debye=[(2, 500.0, 2.0), (3, 1500.0, 1.0)]),     # (n, Theta, weight)
    "gra": dict(name="graphite", rho=2.24,
                atoms_per_cm3=468.0 / (4 / 3 * np.pi * 1e-21),  # N_C = 468 (a/nm)^3
                debye=[(2, 863.0, 1.0), (2, 2504.0, 2.0)]),
}


def n_atoms(mat, a_um):
    return MAT[mat]["atoms_per_cm3"] * 4 / 3 * np.pi * (a_um * UM) ** 3


def grain_mass(mat, a_um):
    return MAT[mat]["rho"] * 4 / 3 * np.pi * (a_um * UM) ** 3


# ------------------------------------------------------------------ dielectric functions
def read_callindex(path):
    """Return lambda (um, ascending) and complex m from one of Draine's callindex files."""
    with open(path) as f:
        lines = f.readlines()
    nw = int(lines[3].split()[0])
    d = np.array([[float(v) for v in ln.split()[:5]] for ln in lines[5:5 + nw]])
    d = d[np.argsort(d[:, 0])]
    return d[:, 0], (1 + d[:, 3]) + 1j * d[:, 4]


def m_on_grid(path, lam=LAM):
    """Refractive index on the grid: Re(m) - 1 linear and Im(m) logarithmic in log lambda."""
    lw, m = read_callindex(path)
    x = np.log(lw)
    re = np.interp(np.log(lam), x, m.real - 1)
    im = np.exp(np.interp(np.log(lam), x, np.log(np.maximum(m.imag, 1e-30))))
    return (1 + re) + 1j * im


_EPS = {}


def _load_eps(lam=LAM):
    """Dielectric functions on a wavelength grid (cached for the standard grid)."""
    std = lam is LAM
    if std and _EPS:
        return _EPS
    D = os.path.join(RAW, "draine")
    out = {"sil": m_on_grid(os.path.join(D, "callindex.out_silD03"), lam) ** 2}
    for o in ("pa", "pe"):
        for s in ("0.01", "0.10"):
            out[o + s] = m_on_grid(os.path.join(D, "callindex.out_C%sD03_%s" % (o, s)), lam) ** 2
    if std:
        _EPS.update(out)
    return out


def refr_index(mat, a_um, lam=LAM):
    """Refractive index on lam; for graphite a pair (E parallel, E perpendicular to c)."""
    E = _load_eps(lam)
    if mat == "sil":
        return np.sqrt(E["sil"])
    w = np.clip(np.log10(a_um / 0.01), 0.0, 1.0)
    return tuple(np.sqrt((1 - w) * E[o + "0.01"] + w * E[o + "0.10"]) for o in ("pa", "pe"))


# ------------------------------------------------------------------ Mie theory
def mie(x, m):
    """Q_ext, Q_sca and g for spheres: BHMIE vectorized over the elements of x and m."""
    x = np.atleast_1d(np.asarray(x, float))
    m = np.broadcast_to(np.asarray(m, complex), x.shape).copy()
    nstop = (x + 4 * x ** (1 / 3) + 2).astype(int)
    NS = int(nstop.max())
    y = m * x
    nmx = int(max(NS, np.abs(y).max())) + 15
    D = np.zeros((NS + 1,) + x.shape, complex)
    Dn = np.zeros(x.shape, complex)
    for n in range(nmx, 0, -1):            # logarithmic derivative, downward
        Dn = n / y - 1 / (Dn + n / y)      # this is D_{n-1}
        if n - 1 <= NS:
            D[n - 1] = Dn
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = psi1 - 1j * chi1
    qs = np.zeros(x.shape)
    qe = np.zeros(x.shape)
    gs = np.zeros(x.shape)
    an1 = bn1 = np.zeros(x.shape, complex)
    with np.errstate(all="ignore"):
        for n in range(1, NS + 1):
            psi = (2 * n - 1) * psi1 / x - psi0
            chi = (2 * n - 1) * chi1 / x - chi0
            xi = psi - 1j * chi
            ta = D[n] / m + n / x
            tb = D[n] * m + n / x
            an = (ta * psi - psi1) / (ta * xi - xi1)
            bn = (tb * psi - psi1) / (tb * xi - xi1)
            on = n <= nstop                 # beyond nstop the upward recursion is unstable
            an = np.where(on, an, 0)
            bn = np.where(on, bn, 0)
            qs += (2 * n + 1) * (np.abs(an) ** 2 + np.abs(bn) ** 2)
            qe += (2 * n + 1) * (an + bn).real
            gs += (2 * n + 1) / (n * (n + 1)) * (an * np.conj(bn)).real
            if n > 1:
                gs += (n - 1) * (n + 1) / n * (an1 * np.conj(an) + bn1 * np.conj(bn)).real
            psi0, psi1, chi0, chi1 = psi1, psi, chi1, chi
            xi1 = psi1 - 1j * chi1
            an1, bn1 = an, bn
    g = np.where(x < 1e-3, 0.0, 2 * gs / qs)   # true g ~ x^2 there; the series loses it
    return 2 / x ** 2 * qe, 2 / x ** 2 * qs, g


def optics(mat, a_um, lam=LAM):
    """Q_abs, Q_sca, g on the wavelength grid for one grain."""
    x = 2 * np.pi * a_um / lam
    m = refr_index(mat, a_um, lam)
    if mat == "sil":
        qe, qs, g = mie(x, m)
        return qe - qs, qs, g
    qe1, qs1, g1 = mie(x, m[0])
    qe2, qs2, g2 = mie(x, m[1])
    qe = qe1 / 3 + 2 * qe2 / 3
    qs = qs1 / 3 + 2 * qs2 / 3
    g = (qs1 * g1 / 3 + 2 * qs2 * g2 / 3) / qs
    return qe - qs, qs, g


# ------------------------------------------------------------------ radiation field
XMAX = 80.0       # h nu / k T beyond which a Planck term is dropped (< e^-70 of the peak)


def B_lam(lam_um, T):
    """Planck function per um: erg s^-1 cm^-2 sr^-1 um^-1."""
    l = lam_um * UM
    x = H * C / (l * KB * T)
    with np.errstate(over="ignore"):
        return np.where(x < XMAX, 2 * H * C ** 2 / l ** 5 / np.expm1(np.minimum(x, XMAX)), 0.0) * UM


STARS = [(7500.0, 1.0e-14), (4000.0, 1.0e-13), (3000.0, 4.0e-13)]   # MMP83 Table A1


def _mmp83_table():
    d = np.loadtxt(os.path.join(DATA, "mmp83_tableA3.txt"))
    return d[:, 0], d[:, 1]


def J_bb(lam_um):
    """4 pi J_lambda of the three diluted blackbodies (erg cm^-2 s^-1 um^-1)."""
    return 4 * np.pi * sum(W * B_lam(lam_um, T) for T, W in STARS)


def J_uv(lam_um):
    """MMP83's UV component: Table A3 minus the blackbodies, log-log interpolated over
    0.0912-0.25 um (where it is positive) and zero elsewhere."""
    tl, tj = _mmp83_table()
    sel = tl <= 0.25
    res = tj[sel] - J_bb(tl[sel])
    lam_um = np.asarray(lam_um, float)
    out = np.exp(np.interp(np.log(lam_um), np.log(tl[sel]), np.log(res)))
    return np.where((lam_um >= LYMAN_UM) & (lam_um <= 0.25), out, 0.0)


def u_star(lam_um=LAM):
    """Starlight energy density for U = 1, erg cm^-3 um^-1."""
    lam_um = np.asarray(lam_um, float)
    return np.where(lam_um >= LYMAN_UM, (J_uv(lam_um) + J_bb(lam_um)) / C, 0.0)


def u_cmb(lam_um=LAM):
    return 4 * np.pi * B_lam(lam_um, T_CMB) / C


def u_field(U, lam_um=LAM):
    return U * u_star(lam_um) + u_cmb(lam_um)


# ------------------------------------------------------------------ heat capacity
def f_debye(n, x):
    """f_n(x) = n int_0^1 y^n / (exp(y/x) - 1) dy (see the module docstring on the prefactor)."""
    x = np.atleast_1d(np.asarray(x, float))
    out = np.empty_like(x)
    # low-x limit: n Gamma(n+1) zeta(n+1) x^(n+1); the neglected upper tail is ~exp(-1/x)
    zeta = {2: 1.2020569031595942, 3: 1.0823232337111382}[n]
    gam = {2: 2.0, 3: 6.0}[n]
    lo = x < 0.03
    out[lo] = n * gam * zeta * x[lo] ** (n + 1)
    if np.any(~lo):
        yg, wg = np.polynomial.legendre.leggauss(200)
        yy = 0.5 * (yg + 1)
        ww = 0.5 * wg
        xx = x[~lo][:, None]
        out[~lo] = n * np.sum(ww * yy ** n / np.expm1(yy / xx), axis=1)
    return out


TGRID_E = 10 ** np.linspace(0.0, np.log10(5000.0), 600)   # 1 - 5000 K, for the E(T) table


def e_per_mode_set(mat, T):
    """E(T) / (N - 2) in erg."""
    return sum(w * KB * th * f_debye(n, np.asarray(T) / th) for n, th, w in MAT[mat]["debye"])


def heatcap_per_atom(mat, T, dT=1e-4):
    """C / (N k) in the large-N limit, by differencing E(T)."""
    T = np.asarray(T, float)
    return (e_per_mode_set(mat, T * (1 + dT)) - e_per_mode_set(mat, T * (1 - dT))) / (2 * dT * T) / KB


class Enthalpy:
    """E(T) and T(E) for one grain, log-log interpolated on TGRID_E."""

    def __init__(self, mat, a_um):
        self.N = n_atoms(mat, a_um)
        self.lT = np.log(TGRID_E)
        self.lE = np.log((self.N - 2) * e_per_mode_set(mat, TGRID_E))

    def E(self, T):
        return np.exp(np.interp(np.log(T), self.lT, self.lE))

    def T(self, E):
        return np.exp(np.interp(np.log(E), self.lE, self.lT))

    def C(self, T, d=1e-4):
        return (self.E(T * (1 + d)) - self.E(T * (1 - d))) / (2 * d * T)


# ------------------------------------------------------------------ one grain in a field
class Grain:
    """A grain of material mat and radius a_um in the field U (optics computed or given)."""

    SIGMA_EQ = 0.01       # treat as equilibrium below this fractional temperature spread
    M = 400               # enthalpy bins, first pass and minimum for the second
    K_COOL = 30.0         # cool for this many mean photon intervals to set the bottom bin

    def __init__(self, mat, a_um, U, opt=None):
        self.mat, self.a, self.U = mat, a_um, U
        self.Qabs, self.Qsca, self.g = optics(mat, a_um) if opt is None else opt
        self.Cabs = np.pi * (a_um * UM) ** 2 * self.Qabs
        self.mass = grain_mass(mat, a_um)
        self.H = Enthalpy(mat, a_um)
        self.u = u_field(U)
        self.us = U * u_star()
        self.pabs_lam = C * self.Cabs * self.u          # absorbed power per um
        self.Pabs = integ_lam(self.pabs_lam)
        self.Teq = self.equilibrium_T()

    # emitted power at T, all wavelengths, or only photons below Emax (erg)
    def L(self, T, Emax=None):
        F = 4 * np.pi * self.Cabs * B_lam(LAM, T)
        if Emax is not None:
            F = np.where(H * C / (LAM * UM) < Emax, F, 0.0)
        return integ_lam(F)

    def equilibrium_T(self):
        lo, hi = np.log(1.0), np.log(5000.0)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if self.L(np.exp(mid)) > self.Pabs:
                hi = mid
            else:
                lo = mid
        return np.exp(0.5 * (lo + hi))

    # photon statistics of the starlight alone
    def photon_stats(self):
        Eph = H * C / (LAM * UM)
        p = C * self.Cabs * self.us
        R = integ_lam(p / Eph)
        S1 = integ_lam(p)
        S2 = integ_lam(p * Eph)
        return R, S1, S2

    def sigma_T(self):
        """Fractional temperature spread for steady fluctuations about T_eq: an
        Ornstein-Uhlenbeck estimate, var(E) = tau_relax/2 * int (h nu)^2 dR."""
        R, S1, S2 = self.photon_stats()
        T = self.Teq
        dL = (self.L(T * 1.001) - self.L(T * 0.999)) / (0.002 * T)
        Cv = self.H.C(T)
        tau = Cv / dL
        sE = np.sqrt(0.5 * tau * S2)
        return sE / Cv / T, sE

    # ---------------------------------------------------------------- stochastic solution
    def _photon_cumulants(self):
        """Cumulative int f dE and int E f dE with f(E) = C_abs u_E, on the ascending photon
        energy grid; trapezoid in ln E."""
        E = (H * C / (LAM * UM))[::-1]
        f = (self.Cabs * self.u * LAM / (H * C / (LAM * UM)))[::-1]    # u_E = u_lambda * lambda / E
        g0 = f * E
        g1 = f * E * E
        I0 = np.concatenate([[0.0], np.cumsum(0.5 * (g0[1:] + g0[:-1]) * DLNLAM)])
        I1 = np.concatenate([[0.0], np.cumsum(0.5 * (g1[1:] + g1[:-1]) * DLNLAM)])
        return E, I0, I1

    def stochastic(self, M=None, Ehi_scale=1.0, verbose=False):
        """Return (T_j, P_j, info) for the thermal-continuous steady state."""
        M = M or self.M
        Hh = self.H
        Teq, Eeq = self.Teq, self.H.E(self.Teq)
        _, sE = self.sigma_T()
        R = self.photon_stats()[0]
        Lcmb = self.L(T_CMB)
        # bottom: cool from T_eq for K_COOL photon intervals, or 10 sigma, whichever is lower
        Ts = np.exp(np.linspace(np.log(Teq), np.log(T_CMB * 1.0005), 200))
        Ln = np.array([self.L(t) for t in Ts]) - Lcmb
        dt = -np.diff(np.log(Ts)) * 0.5 * ((Hh.C(Ts) * Ts / np.maximum(Ln, 1e-300))[1:] +
                                           (Hh.C(Ts) * Ts / np.maximum(Ln, 1e-300))[:-1])
        tcum = np.concatenate([[0], np.cumsum(dt)])
        k = np.searchsorted(tcum, self.K_COOL / R)
        Tcool = Ts[min(k, len(Ts) - 1)]
        Elo = max(Hh.E(T_CMB), min(Eeq - 10 * sE, Hh.E(Tcool)))
        Ehi = Eeq + 13.6 * EV + 10 * sE
        Ehi = Elo + (Ehi - Elo) * Ehi_scale
        for it in range(8):
            Tlo, Thi = Hh.T(Elo), min(Hh.T(Ehi), 4900.0)
            T, P, info = self._solve(np.exp(np.linspace(np.log(Tlo), np.log(Thi), M + 1)))
            if P[-1] < 1e-10 or Thi >= 4900.0:
                break
            Ehi = Elo + 1.5 * (Ehi - Elo)
        info["iterations"] = it + 1
        # Second pass on the range that matters.  The first range is deliberately generous,
        # and for a grain whose distribution is narrow most of its bins are empty; the
        # continuous-cooling steps then add numerical diffusion of order one bin width and
        # broaden P(E).  Keep P above TRIM in both tails, and the upper tail while it still
        # carries more than TRIM_POW of the emitted power.
        jl, ju = self._trim(P, info["down"])
        T2 = self._second_grid(T, P, info["E"], jl, ju, M)
        T, P, info2 = self._solve(T2)
        info2["iterations"] = info["iterations"] + 1
        return T, P, info2

    DELTA = 0.03          # second-pass bin width / mean absorbed photon energy, near the peak
    M_MAX = 1500

    def _second_grid(self, T, P, E, jl, ju, M):
        """Log-T bins over the trimmed range, fine enough that a bin near where the grain
        spends its time is DELTA of the mean absorbed photon energy: the continuous-cooling
        steps add a numerical diffusion whose variance, relative to that of the photon
        heating, is about (bin width) / (mean photon energy)."""
        _, S1, S2 = self.photon_stats()
        eps = S2 / S1
        dlnE = np.gradient(np.log(E))
        k = int(np.argmax(P / dlnE))
        Eref = max(E[k], eps)
        Tref = self.H.T(Eref)
        slope = np.log(self.H.E(Tref * 1.001) / self.H.E(Tref / 1.001)) / np.log(1.001 ** 2)
        dlnT = self.DELTA * eps / (Eref * slope)
        lo, hi = np.log(T[jl]), np.log(T[ju])
        n = int(np.clip(np.ceil((hi - lo) / dlnT), M, self.M_MAX))
        return np.exp(np.linspace(lo, hi, n + 1))

    TRIM, TRIM_POW = 1e-14, 1e-10

    def _trim(self, P, down):
        M = len(P) - 1
        lo = np.cumsum(P)
        hi = np.cumsum(P[::-1])[::-1]
        pw = P * down
        hipow = np.cumsum(pw[::-1])[::-1] / pw.sum()
        jl = max(0, int(np.argmax(lo > self.TRIM)) - 1)
        keep = (hi > self.TRIM) | (hipow > self.TRIM_POW)
        ju = min(M, int(np.nonzero(keep)[0].max()) + 1)
        return jl, ju

    def _solve(self, T):
        Hh = self.H
        M = len(T) - 1
        E = Hh.E(T)
        edge = np.empty(M + 2)
        edge[1:-1] = 0.5 * (E[1:] + E[:-1])
        edge[0] = max(0.0, E[0] - 0.5 * (E[1] - E[0]))
        edge[-1] = E[-1] + 0.5 * (E[-1] - E[-2])
        Emin, Emax = edge[:-1], edge[1:]
        dE = Emax - Emin
        Eg, I0, I1 = self._photon_cumulants()
        Etop = Eg[-1]

        def cum(Ix, w):                 # linear in E; zero below the grid
            return np.interp(w, Eg, Ix)

        def ramp_up(a, b, w0):      # int_a^b (E - w0) f dE
            return (cum(I1, b) - cum(I1, a)) - w0 * (cum(I0, b) - cum(I0, a))

        def ramp_dn(a, b, w4):      # int_a^b (w4 - E) f dE
            return w4 * (cum(I0, b) - cum(I0, a)) - (cum(I1, b) - cum(I1, a))

        A = np.zeros((M + 1, M + 1))    # A[u, l]: rate l -> u
        for l in range(M):
            u = np.arange(l + 1, M)
            W1 = Emin[u] - Emax[l]
            W4 = Emax[u] - Emin[l]
            ok = W1 < Etop
            u, W1, W4 = u[ok], W1[ok], W4[ok]
            if u.size:
                a1 = Emin[u] - Emin[l]
                a2 = Emax[u] - Emax[l]
                W2 = np.minimum(a1, a2)
                W3 = np.maximum(a1, a2)
                W1c = np.maximum(W1, 0.0)
                G = (ramp_up(W1c, W2, W1) + np.minimum(dE[l], dE[u]) * (cum(I0, W3) - cum(I0, W2))
                     + ramp_dn(W3, W4, W4)) / (dE[u] * dE[l])
                A[u, l] = C * dE[u] / (E[u] - E[l]) * G
            # top bin takes everything that would go beyond it
            W1 = Emin[M] - Emax[l]
            Wc = Emin[M] - Emin[l]
            W1c = max(W1, 0.0)
            A[M, l] = C / (E[M] - E[l]) * (ramp_up(W1c, Wc, W1) / (Wc - W1) + (cum(I0, 1e300) - cum(I0, Wc)))
            # intrabin absorptions in bin l promote to l + 1
            A[l + 1, l] += C / (E[l + 1] - E[l]) * ramp_dn(0.0, dE[l], dE[l]) / dE[l]
        # continuous cooling u -> u - 1
        down = np.array([self.L(T[u], Emax=E[u]) for u in range(M + 1)])
        Dn = np.zeros(M + 1)
        Dn[1:] = down[1:] / (E[1:] - E[:-1])
        # Guhathakurta & Draine recursion
        B = np.cumsum(A[::-1], axis=0)[::-1]        # B[j, i] = sum_{u >= j} A[u, i]
        X = np.zeros(M + 1)
        X[0] = 1.0
        for j in range(1, M + 1):
            X[j] = np.dot(B[j, :j], X[:j]) / Dn[j]
            if X[j] > 1e250:
                X[:j + 1] *= 1e-250
        P = X / X.sum()
        info = dict(E=E, down=down, Pemit=np.dot(P, down))
        return T, P, info

    # ---------------------------------------------------------------- emission
    def emission(self, stochastic=True):
        """Emitted power per um per grain, the temperatures, P(T) and a flag."""
        sT, _ = self.sigma_T()
        if not stochastic or sT < self.SIGMA_EQ:
            return 4 * np.pi * self.Cabs * B_lam(LAM, self.Teq), np.array([self.Teq]), np.array([1.0]), False
        T, P, info = self.stochastic()
        E = info["E"]
        Eph = H * C / (LAM * UM)
        em = np.zeros(NLAM)
        for Tj, Pj, Ej in zip(T, P, E):
            if Pj < 1e-30:
                continue
            em += Pj * np.where(Eph < Ej, B_lam(LAM, Tj), 0.0)
        return 4 * np.pi * self.Cabs * em, T, P, True


# ------------------------------------------------------------------ tables
def optics_table(mat, a_grid=A_GRID):
    """Q_abs, Q_sca, g on (a_grid, LAM)."""
    out = np.zeros((3, len(a_grid), NLAM))
    for i, a in enumerate(a_grid):
        out[:, i] = optics(mat, a)
    return out


def enthalpy_table(mat):
    """E(T)/(N - 2) on TGRID_E (erg)."""
    return e_per_mode_set(mat, TGRID_E)


# ------------------------------------------------------------------ size distributions
def powerlaw_weights(mat, q, amin, amax, abund_ppm, a_grid=A_GRID, nsub=8):
    """Nodes and weights for int dn/da F(a) da with dn/da = A a^-q on [amin, amax],
    normalized so the grains hold abund_ppm (Si or C atoms per 1e6 H) of their key element.
    Returns (a_nodes, w) with sum(w * F(a_nodes)) = (1/n_H) int dn/da F da, F per grain;
    nodes are a fine log grid (nsub per table interval) on which F is interpolated."""
    la = np.log(a_grid)
    n = int(np.ceil((np.log(amax) - np.log(amin)) / (la[1] - la[0]) * nsub)) + 1
    lx = np.linspace(np.log(amin), np.log(amax), max(n, 2))
    ax = np.exp(lx)
    tw = np.full(lx.size, lx[1] - lx[0])
    tw[0] = tw[-1] = 0.5 * (lx[1] - lx[0])
    shape = ax ** (1 - q) * tw                      # a dn/da d ln a, up to A
    # key-element atoms per grain: silicate one Si per formula unit (7 atoms), graphite all C
    per = n_atoms(mat, ax) / (7.0 if mat == "sil" else 1.0)
    Anorm = abund_ppm * 1e-6 / np.sum(shape * per)
    return ax, Anorm * shape


def interp_in_a(table, a_nodes, a_grid=A_GRID):
    """Interpolate a (len(a_grid), ...) table to a_nodes, linearly in ln a."""
    la = np.log(a_grid)
    x = np.log(a_nodes)
    i = np.clip(np.searchsorted(la, x) - 1, 0, len(la) - 2)
    t = ((x - la[i]) / (la[i + 1] - la[i])).reshape((-1,) + (1,) * (table.ndim - 1))
    return (1 - t) * table[i] + t * table[i + 1]


# ------------------------------------------------------------------ a whole size distribution
MAG = 2.5 / np.log(10)              # A = 1.0857 tau

# MRN (Mathis, Rumpl & Nordsieck 1977) as normalized by Draine & Lee (1984) and quoted by
# WD01: (1/n_H) dn/da = C a^-3.5 for 50 A < a < 0.25 um, C = 10^-25.11 (silicate) and
# 10^-25.13 (graphite) cm^2.5.  The abundances below are the key-element atoms that
# normalization puts in grains (computed by mrn_abundance()).
MRN_C = {"sil": 10 ** -25.11, "gra": 10 ** -25.13}


def mrn_abundance(mat, q=3.5, amin=0.005, amax=0.25):
    """ppm of Si (silicate) or C (graphite) per H for the DL84 MRN normalization."""
    ax, w = powerlaw_weights(mat, q, amin, amax, 1.0)        # weights for 1 ppm
    per = n_atoms(mat, ax) / (7.0 if mat == "sil" else 1.0)
    # with dn/da = C a^-q (a in cm): atoms per H = C int a^-q N(a) da
    lx = np.log(ax)
    tw = np.full(lx.size, lx[1] - lx[0])
    tw[0] = tw[-1] = 0.5 * (lx[1] - lx[0])
    acm = ax * UM
    return 1e6 * MRN_C[mat] * np.sum(acm ** (1 - q) * tw * per)


SIZE_BINS = [(0.0, 0.01), (0.01, 0.1), (0.1, np.inf)]   # um, for splitting the emission


def distribution(params, U, tables, emis_cache=None):
    """Extinction, scattering, emission and mass per H for power-law size distributions.

    params: {mat: dict(q=, amin=, amax=, abund=)} for the materials that are switched on.
    tables: {mat: optics_table(mat)} on A_GRID.
    Returns per-H quantities on LAM: 'ext', 'sca' (cm^2/H), 'em' (lambda I_lambda / N_H,
    erg s^-1 sr^-1 H^-1), each a dict over materials, plus 'em_size' (same, split by
    SIZE_BINS), 'mass' (g/H) and 'pabs' (erg s^-1 H^-1)."""
    out = dict(ext={}, sca={}, em={}, mass={}, pabs={}, em_size=[np.zeros(NLAM) for _ in SIZE_BINS])
    for mat, p in params.items():
        ax, w = powerlaw_weights(mat, p["q"], p["amin"], p["amax"], p["abund"])
        Q = interp_in_a(np.moveaxis(tables[mat], 0, 1), ax)          # (n, 3, NLAM)
        geo = np.pi * (ax * UM) ** 2
        out["ext"][mat] = np.sum((w * geo)[:, None] * (Q[:, 0] + Q[:, 1]), axis=0)
        out["sca"][mat] = np.sum((w * geo)[:, None] * Q[:, 1], axis=0)
        out["mass"][mat] = np.sum(w * grain_mass(mat, ax))
        # emission per grain at the table nodes that bracket [amin, amax], then
        # interpolated in ln a as (power per grain) / a^3
        i0 = max(0, np.searchsorted(A_GRID, p["amin"]) - 1)
        i1 = min(len(A_GRID) - 1, np.searchsorted(A_GRID, p["amax"]))
        node = np.zeros((len(A_GRID), NLAM))
        pab = np.zeros(len(A_GRID))
        for i in range(i0, i1 + 1):
            key = (mat, i, U)
            if emis_cache is not None and key in emis_cache:
                node[i], pab[i] = emis_cache[key]
                continue
            g = Grain(mat, A_GRID[i], U, opt=tuple(tables[mat][:, i]))
            em, _, _, _ = g.emission()
            node[i], pab[i] = em, g.Pabs
            if emis_cache is not None:
                emis_cache[key] = (em, g.Pabs)
        per_a3 = interp_in_a(node / A_GRID[:, None] ** 3, ax) * (ax ** 3)[:, None]
        pab_a = interp_in_a((pab / A_GRID ** 3)[:, None], ax)[:, 0] * ax ** 3
        out["em"][mat] = LAM * np.sum(w[:, None] * per_a3, axis=0) / (4 * np.pi)
        out["pabs"][mat] = np.sum(w * pab_a)
        for k, (lo, hi) in enumerate(SIZE_BINS):
            sel = (ax >= lo) & (ax < hi)
            out["em_size"][k] += LAM * np.sum((w * sel)[:, None] * per_a3, axis=0) / (4 * np.pi)
    return out
