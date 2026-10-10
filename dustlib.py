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


class Grid:
    """A wavelength grid (um, ascending, not necessarily uniform) with trapezoid weights in
    ln(lambda): int F dlambda = sum(wln * lam * F)."""

    def __init__(self, lam):
        self.lam = np.asarray(lam, float)
        self.n = len(self.lam)
        d = np.diff(np.log(self.lam))
        self.wln = np.zeros(self.n)
        self.wln[:-1] += 0.5 * d
        self.wln[1:] += 0.5 * d
        self.eph = H * C / (self.lam * UM)

    def integ(self, F):
        return np.sum(self.wln * self.lam * F, axis=-1)


def _fine_lambda():
    """The PAH grid: 80 per decade below 1 um and above 20 um, 200 per decade over 1-3 um
    for the near-infrared bands of PAH cations, 500 per decade over 3-20 um, where the
    narrowest bands have FWHM ~1.2% (trapezoid error on their integrals < 1e-3)."""
    seg = lambda lo, hi, n: 10 ** np.arange(np.log10(lo), np.log10(hi) - 1e-12, 1.0 / n)
    return np.concatenate([seg(LYMAN_UM, 1.0, 80), seg(1.0, 3.0, 200), seg(3.0, 20.0, 500),
                           seg(20.0, 1e4, 80), [1e4]])


COARSE = Grid(LAM)              # silicate and graphite
FINE = Grid(_fine_lambda())     # PAHs, and the distribution sums
LAMF = FINE.lam


def to_fine(arr):
    """Coarse-grid array onto the fine grid, linear in ln(lambda)."""
    return np.interp(np.log(LAMF), np.log(LAM), arr)


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


def e_first_mode(mat, N):
    """Energy of the lowest vibrational mode (erg), from the model mode spectra of DL01:
    hbar omega_1 = k Theta [(1 - beta)(1 - delta_1)/N_m + beta]^(1/n), delta_1 = 1/2.
    Carbon (PAH, graphite): out-of-plane C-C modes, n = 2, Theta = 863 K, N_m = N_C - 2, beta
    from DL01 Eq. (beta_graph).  Silicate: the lower of the n = 2 (Theta = 500 K, N_m = 2(N - 2))
    and n = 3 (Theta = 1500 K, N_m = N - 2) sets, beta from Eq. (beta_sil).  Below this energy a
    grain is in its vibrational ground state."""
    if mat == "sil":
        out = []
        for n, th, Nm in ((2, 500.0, 2 * (N - 2)), (3, 1500.0, N - 2)):
            beta = (Nm ** (1 - n / 3) - 1) / (2 * Nm - 1)
            out.append(KB * th * ((1 - beta) * 0.5 / Nm + beta) ** (1.0 / n))
        return min(out)
    NC, Nm = N, N - 2
    if NC <= 54:
        beta = 0.0
    elif NC <= 102:
        beta = (NC - 54) / 52 / (2 * Nm - 1)
    else:
        beta = ((NC - 2) / 52 * (102 / NC) ** (2 / 3) - 1) / (2 * Nm - 1)
    return KB * 863.0 * ((1 - beta) * 0.5 / Nm + beta) ** 0.5


class Enthalpy:
    """E(T) and T(E) for one grain, log-log interpolated on TGRID_E."""

    def __init__(self, mat, a_um):
        self.N = n_atoms(mat, a_um)
        self.E1 = e_first_mode(mat, self.N)
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

    def __init__(self, mat, a_um, U, opt=None, charge="neu"):
        self.mat, self.a, self.U, self.charge = mat, a_um, U, charge
        self.G = FINE if mat == "pah" else COARSE
        if opt is None:
            opt = pah_optics(a_um, charge) if mat == "pah" else optics(mat, a_um)
        self.Qabs, self.Qsca, self.g = opt
        self.Cabs = np.pi * (a_um * UM) ** 2 * self.Qabs
        if mat == "pah":
            self.H = PAHEnthalpy(a_um)
            self.mass = self.H.mass
        else:
            self.H = Enthalpy(mat, a_um)
            self.mass = grain_mass(mat, a_um)
        self.u = u_field(U, self.G.lam)
        self.us = U * u_star(self.G.lam)
        self.pabs_lam = C * self.Cabs * self.u          # absorbed power per um
        self.Pabs = self.G.integ(self.pabs_lam)
        self.Teq = self.equilibrium_T()

    # emitted power at T, all wavelengths, or only photons below Emax (erg)
    def L(self, T, Emax=None):
        F = 4 * np.pi * self.Cabs * B_lam(self.G.lam, T)
        if Emax is not None:
            F = np.where(self.G.eph < Emax, F, 0.0)
        return self.G.integ(F)

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
        Eph = self.G.eph
        p = C * self.Cabs * self.us
        R = self.G.integ(p / Eph)
        S1 = self.G.integ(p)
        S2 = self.G.integ(p * Eph)
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
        lam = self.G.lam
        E = self.G.eph[::-1]
        f = (self.Cabs * self.u * lam / self.G.eph)[::-1]    # u_E = u_lambda * lambda / E
        g0 = f * E
        g1 = f * E * E
        dl = np.diff(np.log(E))
        I0 = np.concatenate([[0.0], np.cumsum(0.5 * (g0[1:] + g0[:-1]) * dl)])
        I1 = np.concatenate([[0.0], np.cumsum(0.5 * (g1[1:] + g1[:-1]) * dl)])
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
        # never below the lowest vibrational mode (nor below two quanta of the longest
        # wavelength on the grid, where no radiative transition is possible): energies under
        # E1 belong to the vibrational ground state, which the bottom bin then stands for
        Erad = 2 * self.G.eph.min()
        Elo = max(Hh.E(T_CMB), Erad, Hh.E1, min(Eeq - 10 * sE, Hh.E(Tcool)))
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
        Dn[1:] = np.maximum(down[1:], 1e-300) / (E[1:] - E[:-1])
        # Guhathakurta & Draine recursion
        B = np.cumsum(A[::-1], axis=0)[::-1]        # B[j, i] = sum_{u >= j} A[u, i]
        X = np.zeros(M + 1)
        X[0] = 1.0
        for j in range(1, M + 1):
            X[j] = np.dot(B[j, :j], X[:j]) / Dn[j]
            if X[j] > 1e250:
                X[:j + 1] *= 1e-250
        P = X / X.sum()
        # the bottom bin has no downward transition: it is the ground state (or empty), so it
        # radiates nothing, and the power balance is sum_{j>=1} P_j L_j
        info = dict(E=E, down=down, Pemit=np.dot(P[1:], down[1:]), ground=P[0] if E[0] <= 1.0001 * self.H.E1 else 0.0)
        return T, P, info

    # ---------------------------------------------------------------- emission
    def emission(self, stochastic=True):
        """Emitted power per um per grain, the temperatures, P(T) and a flag."""
        lam = self.G.lam
        sT, _ = self.sigma_T()
        if not stochastic or sT < self.SIGMA_EQ:
            return 4 * np.pi * self.Cabs * B_lam(lam, self.Teq), np.array([self.Teq]), np.array([1.0]), False
        T, P, info = self.stochastic()
        E = info["E"]
        em = np.zeros(self.G.n)
        for j, (Tj, Pj, Ej) in enumerate(zip(T, P, E)):
            if j == 0 or Pj < 1e-30:                 # the bottom bin radiates nothing
                continue
            em += Pj * np.where(self.G.eph < Ej, B_lam(lam, Tj), 0.0)
        return 4 * np.pi * self.Cabs * em, T, P, True


# ------------------------------------------------------------------ PAHs
# Absorption per C atom: Li & Draine (2001, LD01) Eqs. (5)-(11) and their Drude profiles,
# with the band parameters and changes of Draine & Li (2007, DL07, Sect. 2, Table 1),
# including the near-infrared continuum and bands of PAH cations.  Each band is
# (lambda_j [um], gamma_j, sigma_int neutral, sigma_int ionized [1e-20 cm/C], times H/C?).
# ⚠ DL07 Table 1 row 29 (18.92 um) reads gamma = 0.10, sigma = 0.10 (neutral), 0.17
# (ionized) in the arXiv source, after a commented-out older row; used as printed.
PAH_UV = [(0.0722, 0.195, 7.97e7), (0.2175, 0.217, 1.23e7)]          # S_1, S_2
PAH_BANDS = {
    "DL07": [(1.050, 0.055, 0.0, 2.0e4, False), (1.260, 0.11, 0.0, 0.078, False),
             (1.905, 0.09, 0.0, -146.5, False), (3.300, 0.012, 394.0, 89.4, True),
             (5.270, 0.034, 2.5, 20.0, False), (5.700, 0.035, 4.0, 32.0, False),
             (6.220, 0.030, 29.4, 235.0, False), (6.690, 0.070, 7.35, 59.0, False),
             (7.417, 0.126, 20.8, 181.0, False), (7.598, 0.044, 18.1, 163.0, False),
             (7.850, 0.053, 21.9, 197.0, False), (8.330, 0.052, 6.94, 48.0, True),
             (8.610, 0.039, 27.8, 194.0, True), (10.68, 0.020, 0.3, 0.3, True),
             (11.23, 0.012, 18.9, 17.7, True), (11.33, 0.032, 52.0, 49.0, True),
             (11.99, 0.045, 24.2, 20.5, True), (12.62, 0.042, 35.0, 31.0, True),
             (12.69, 0.013, 1.3, 1.3, True), (13.48, 0.040, 8.0, 8.0, True),
             (14.19, 0.025, 0.45, 0.45, False), (15.90, 0.020, 0.04, 0.04, False),
             (16.45, 0.014, 0.5, 0.5, False), (17.04, 0.065, 2.22, 2.22, False),
             (17.375, 0.012, 0.11, 0.11, False), (17.87, 0.016, 0.067, 0.067, False),
             (18.92, 0.10, 0.10, 0.17, False), (15.0, 0.8, 50.0, 50.0, False)],
    # LD01 Table 1 with E_6.2 = 3, E_7.7 = E_8.6 = 2 (their diffuse-ISM values); used only
    # to check this code against Draine's tabulated LD01 cross sections
    "LD01": [(3.3, 0.012, 197.0, 44.7, True), (6.2, 0.032, 3 * 19.6, 3 * 157.0, False),
             (7.7, 0.091, 2 * 60.9, 2 * 548.0, False), (8.6, 0.047, 2 * 34.7, 2 * 242.0, True),
             (11.3, 0.018, 427.0 / 3, 400.0 / 3, True), (11.9, 0.025, 72.7 / 3, 61.4 / 3, True),
             (12.7, 0.024, 167.0 / 3, 149.0 / 3, True), (16.4, 0.010, 5.52, 5.52, False),
             (18.3, 0.036, 6.04, 6.04, False), (21.2, 0.038, 10.8, 10.8, False),
             (23.1, 0.046, 2.78, 2.78, False), (26.0, 0.69, 15.2, 15.2, False)],
}
NC_PER_NM3 = {"DL07": 460.0, "LD01": 470.2}     # N_C = k (a / 1 nm)^3 (DL07 Eq. 2; LD01 note 2)
PAH_AMIN = 3.55e-4                              # um; N_C = 20, the smallest PAH (DL07)
HC_OVER_K = 1.4387769                           # cm K
CH_MODES = [886.0, 1161.0, 3030.0]              # cm^-1, C-H bending/stretching (DL01)
M_C, M_H = 12.011 * AMU, 1.008 * AMU


def nc_pah(a_um, model="DL07"):
    return NC_PER_NM3[model] * (a_um / 1e-3) ** 3


def hc_ratio(NC):
    """H/C for compact PAHs (LD01 Eq. 4, DL07 Eq. 3)."""
    return np.where(NC <= 25, 0.5, np.where(NC <= 100, 0.5 * np.sqrt(25.0 / np.maximum(NC, 1)), 0.25))


def drude(lam_um, lj, gj, sig):
    """S_j per C atom (cm^2) for sigma_int in 1e-20 cm/C."""
    return 2 / np.pi * gj * lj * UM * sig * 1e-20 / ((lam_um / lj - lj / lam_um) ** 2 + gj ** 2)


def pah_cabs_per_c(lam_um, NC, charge, model="DL07"):
    """C_abs / N_C (cm^2) of a PAH with N_C carbon atoms, neutral ("neu") or ionized ("ion")."""
    lam = np.asarray(lam_um, float)
    x = 1.0 / lam
    ion = charge == "ion"
    hc = float(hc_ratio(NC))
    S1 = drude(lam, *PAH_UV[0])
    S2 = drude(lam, *PAH_UV[1])
    out = np.zeros_like(lam)
    s = (x > 10) & (x < 15)
    out[s] = S1[s] + (-3.0 + 1.35 * x[s]) * 1e-18
    s = (x > 7.7) & (x <= 10)
    out[s] = (66.302 - 24.367 * x[s] + 2.950 * x[s] ** 2 - 0.1057 * x[s] ** 3) * 1e-18
    s = (x > 5.9) & (x <= 7.7)
    out[s] = S2[s] + (1.8687 + 0.1905 * x[s] + 0.4175 * (x[s] - 5.9) ** 2 + 0.04370 * (x[s] - 5.9) ** 3) * 1e-18
    s = (x > 3.3) & (x <= 5.9)
    out[s] = S2[s] + (1.8687 + 0.1905 * x[s]) * 1e-18
    s = x <= 3.3
    M = (0.3 if NC < 40 else 0.4) * NC                  # fused benzenoid rings (LD01 App. A)
    lc = 1.0 / (2.282 / np.sqrt(M) + 0.889) if ion else 1.0 / (3.804 / np.sqrt(M) + 1.052)
    y = lc / lam[s]
    cut = np.arctan(1e3 * (y - 1) ** 3 / y) / np.pi + 0.5
    ir = 34.58e-18 * 10 ** (-3.431 / x[s]) * cut
    for lj, gj, sn, si, byhc in PAH_BANDS[model]:
        ir = ir + drude(lam[s], lj, gj, (si if ion else sn) * (hc if byhc else 1.0))
    out[s] = ir
    if ion and model == "DL07":                         # DL07 Eq. (1), Mattioda et al. (2005)
        out = out + 3.5e-19 * 10 ** (-1.45 / x) * np.exp(-0.1 * x ** 2)
    return out


def xi_gra(a_um):
    """Graphite weight in the PAH-graphite blend (DL07 Eqs. 4-5; LD01 Eq. 3)."""
    return np.where(a_um <= 0.005, 0.01, 0.01 + 0.99 * (1 - (0.005 / a_um) ** 3))


def graphite_on_fine(a_um):
    """Graphite C_abs, C_sca (cm^2) and g on the fine grid for a sphere of radius a: Mie theory
    on the coarse grid at max(a, 0.001 um), scaled to a in the small-particle limit
    (C_abs ~ a^3, C_sca ~ a^6), then interpolated in ln(lambda)."""
    ae = max(a_um, 1e-3)
    qa, qs, g = optics("gra", ae)
    geo = np.pi * (ae * UM) ** 2
    s = a_um / ae
    return to_fine(geo * qa) * s ** 3, to_fine(geo * qs) * s ** 6, to_fine(g)


def pah_optics(a_um, charge, model="DL07"):
    """Q_abs, Q_sca, g on the fine grid for a PAH-graphite grain (DL07 Eq. 4).  Scattering is
    that of a graphite sphere (negligible at these sizes)."""
    NC = nc_pah(a_um, model)
    ca_gra, cs_gra, g = graphite_on_fine(a_um)
    cabs = (1 - xi_gra(a_um)) * NC * pah_cabs_per_c(LAMF, NC, charge, model) + xi_gra(a_um) * ca_gra
    # ⚠ For cations the negative 1.905 um "resonance" of DL07 outweighs everything else and
    # leaves C_abs < 0 near 1.9 um (-1.6e-22 cm^2/C); DL07 intend it only to suppress
    # absorption there, so the total is floored at zero.
    cabs = np.maximum(cabs, 0.0)
    geo = np.pi * (a_um * UM) ** 2
    return cabs / geo, cs_gra / geo, g


class PAHEnthalpy(Enthalpy):
    """E(T) for a PAH: N_H Einstein oscillators for each C-H mode plus the 2-D Debye C-C modes
    (DL01 Eq. eq:E_PAH_Debye), with H/C from hc_ratio."""

    def __init__(self, a_um, model="DL07"):
        self.NC = nc_pah(a_um, model)
        self.NH = float(hc_ratio(self.NC)) * self.NC
        self.N = self.NC + self.NH
        self.E1 = e_first_mode("pah", self.NC)
        self.mass = self.NC * M_C + self.NH * M_H
        T = TGRID_E
        ech = sum(KB * HC_OVER_K * w / np.expm1(np.minimum(HC_OVER_K * w / T, 700)) for w in CH_MODES)
        self.lT = np.log(T)
        self.lE = np.log(self.NH * ech + (self.NC - 2) * e_per_mode_set("gra", T))


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

# PAHs: a log-normal in radius (DL07 Eq. 6, one component), dn/da = (n0/a) exp[-ln^2(a/a0)/2 sigma^2]
# for a >= PAH_AMIN, normalized to abund ppm of C per H.  Emission and optics are computed at
# PAH_GRID nodes (15 per decade) and interpolated as (per-grain quantity)/a^3, as for the others.
PAH_GRID = 10 ** np.linspace(np.log10(PAH_AMIN), np.log10(0.05), 34)


def pah_amax(a0, sigma):
    return float(np.clip(a0 * np.exp(4 * sigma), 2 * a0, 0.05))


def lognormal_weights(a0, sigma, abund, nper=120):
    """Nodes and weights for (1/n_H) int dn/da F da over [PAH_AMIN, a0 e^{4 sigma}]."""
    lo, hi = np.log(PAH_AMIN), np.log(pah_amax(a0, sigma))
    n = max(int(np.ceil((hi - lo) / np.log(10) * nper)) + 1, 2)
    lx = np.linspace(lo, hi, n)
    ax = np.exp(lx)
    tw = np.full(n, lx[1] - lx[0])
    tw[0] = tw[-1] = 0.5 * (lx[1] - lx[0])
    shape = np.exp(-0.5 * ((lx - np.log(a0)) / sigma) ** 2) * tw      # dn/dln a dln a, up to n0
    A = abund * 1e-6 / np.sum(shape * nc_pah(ax))
    return ax, A * shape


def interp_nodes(nodes, vals, a):
    """Interpolate vals (len(nodes), ...) to radii a, linearly in ln a (clamped)."""
    la = np.log(nodes)
    x = np.clip(np.log(a), la[0], la[-1])
    i = np.clip(np.searchsorted(la, x) - 1, 0, len(la) - 2)
    t = ((x - la[i]) / (la[i + 1] - la[i])).reshape((-1,) + (1,) * (vals.ndim - 1))
    return (1 - t) * vals[i] + t * vals[i + 1]


def pah_node_data(U, cache=None):
    """Per-node optics (C_abs neutral/ion, C_sca) and emission (neutral/ion) for the PAHs."""
    key = ("pah", U)
    if cache is not None and key in cache:
        return cache[key]
    n = len(PAH_GRID)
    d = dict(cabs={c: np.zeros((n, FINE.n)) for c in ("neu", "ion")}, csca=np.zeros((n, FINE.n)),
             em={c: np.zeros((n, FINE.n)) for c in ("neu", "ion")}, pabs={c: np.zeros(n) for c in ("neu", "ion")})
    for i, a in enumerate(PAH_GRID):
        geo = np.pi * (a * UM) ** 2
        for c in ("neu", "ion"):
            opt = pah_optics(a, c)
            g = Grain("pah", a, U, opt=opt, charge=c)
            d["cabs"][c][i] = geo * opt[0]
            d["csca"][i] = geo * opt[1]
            d["em"][c][i] = g.emission()[0]
            d["pabs"][c][i] = g.Pabs
    if cache is not None:
        cache[key] = d
    return d


def distribution(params, U, tables, emis_cache=None):
    """Extinction, scattering, emission and mass per H for the size distributions in params.

    params: {mat: dict(...)} for the populations switched on: "sil", "gra" with q, amin, amax,
    abund (power laws), and "pah" with a0, sigma, abund, fion (log-normal; fion = ionized
    fraction).  tables: {mat: optics_table(mat)} on A_GRID for silicate and graphite.
    Returns per-H quantities on the fine grid LAMF: 'ext', 'sca' (cm^2/H), 'em' (lambda I_lambda
    / N_H, erg s^-1 sr^-1 H^-1), each a dict over materials; 'ext_size', 'sca_size', 'em_size'
    (summed over materials, split by SIZE_BINS); 'mass' (g/H) and 'pabs' (erg s^-1 H^-1).
    Silicate and graphite are computed on the coarse grid and interpolated onto the fine one."""
    nF = FINE.n
    out = dict(ext={}, sca={}, em={}, mass={}, pabs={}, em_size=[np.zeros(nF) for _ in SIZE_BINS],
               ext_size=[np.zeros(nF) for _ in SIZE_BINS], sca_size=[np.zeros(nF) for _ in SIZE_BINS])
    for mat, p in params.items():
        if mat == "pah":
            ax, w = lognormal_weights(p["a0"], p["sigma"], p["abund"])
            d = pah_node_data(U, emis_cache)
            f = p["fion"]
            a3 = (ax ** 3)[:, None]
            n3 = (PAH_GRID ** 3)[:, None]
            cabs = interp_nodes(PAH_GRID, ((1 - f) * d["cabs"]["neu"] + f * d["cabs"]["ion"]) / n3, ax) * a3
            csca = interp_nodes(PAH_GRID, d["csca"] / n3, ax) * a3
            em = interp_nodes(PAH_GRID, ((1 - f) * d["em"]["neu"] + f * d["em"]["ion"]) / n3, ax) * a3
            pab = interp_nodes(PAH_GRID, (((1 - f) * d["pabs"]["neu"] + f * d["pabs"]["ion"]) / PAH_GRID ** 3)[:, None], ax)[:, 0] * ax ** 3
            ext_k, sca_k, em_k = w[:, None] * (cabs + csca), w[:, None] * csca, LAMF * w[:, None] * em / (4 * np.pi)
            NC = nc_pah(ax)
            out["mass"][mat] = np.sum(w * (NC * M_C + hc_ratio(NC) * NC * M_H))
            out["pabs"][mat] = np.sum(w * pab)
        else:
            ax, w = powerlaw_weights(mat, p["q"], p["amin"], p["amax"], p["abund"])
            Q = interp_in_a(np.moveaxis(tables[mat], 0, 1), ax)          # (n, 3, NLAM)
            geo = np.pi * (ax * UM) ** 2
            ext_k = np.array([to_fine(v) for v in (w * geo)[:, None] * (Q[:, 0] + Q[:, 1])])
            sca_k = np.array([to_fine(v) for v in (w * geo)[:, None] * Q[:, 1]])
            out["mass"][mat] = np.sum(w * grain_mass(mat, ax))
            # emission per grain at the table nodes that bracket [amin, amax], interpolated in
            # ln a as (power per grain) / a^3
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
            em_k = np.array([to_fine(v) for v in LAM * w[:, None] * per_a3 / (4 * np.pi)])
            out["pabs"][mat] = np.sum(w * pab_a)
        out["ext"][mat] = ext_k.sum(axis=0)
        out["sca"][mat] = sca_k.sum(axis=0)
        out["em"][mat] = em_k.sum(axis=0)
        for k, (lo, hi) in enumerate(SIZE_BINS):
            sel = (ax >= lo) & (ax < hi)
            out["ext_size"][k] += ext_k[sel].sum(axis=0)
            out["sca_size"][k] += sca_k[sel].sum(axis=0)
            out["em_size"][k] += em_k[sel].sum(axis=0)
    return out


# ------------------------------------------------------------------ reference extinction curve
# Cardelli, Clayton & Mathis (1989, ApJ 345, 245), Eqs. (1)-(5), read from the ADS scan
# (rawdata/ccm89): A_lambda / A_V = a(x) + b(x) / R_V for 0.3 <= x = 1/lambda <= 10 um^-1.
# Per H it is scaled by A_V / N_H = R_V / (N_H / E(B-V)) with N_H / E(B-V) = 5.8e21 cm^-2 mag^-1
# (Bohlin, Savage & Drake 1978, abstract).
NH_EBV = 5.8e21


def ccm89(x, Rv=3.1):
    """A_lambda / A_V; NaN outside 0.3-10 um^-1."""
    x = np.atleast_1d(np.asarray(x, float))
    a = np.full(x.shape, np.nan)
    b = np.full(x.shape, np.nan)
    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61
    op = (x >= 1.1) & (x < 3.3)
    y = x[op] - 1.82
    a[op] = np.polyval([0.32999, -0.77530, 0.01979, 0.72085, -0.02427, -0.50447, 0.17699, 1.0], y)
    b[op] = np.polyval([-2.09002, 5.30260, -0.62251, -5.38434, 1.07233, 2.28305, 1.41338, 0.0], y)
    uv = (x >= 3.3) & (x <= 8)
    xu = x[uv]
    fa = np.where(xu >= 5.9, -0.04473 * (xu - 5.9) ** 2 - 0.009779 * (xu - 5.9) ** 3, 0.0)
    fb = np.where(xu >= 5.9, 0.2130 * (xu - 5.9) ** 2 + 0.1207 * (xu - 5.9) ** 3, 0.0)
    a[uv] = 1.752 - 0.316 * xu - 0.104 / ((xu - 4.67) ** 2 + 0.341) + fa
    b[uv] = -3.090 + 1.825 * xu + 1.206 / ((xu - 4.62) ** 2 + 0.263) + fb
    fu = (x > 8) & (x <= 10)
    z = x[fu] - 8
    a[fu] = -1.073 - 0.628 * z + 0.137 * z ** 2 - 0.070 * z ** 3
    b[fu] = 13.670 + 4.257 * z - 0.420 * z ** 2 + 0.374 * z ** 3
    return a + b / Rv
