# Dust Explorer

An interactive page on interstellar dust grains: change the grain size, or the size distribution, and see the
extinction, the albedo, the grain temperature and the infrared emission that result. Two tabs:

- **Single grain**: efficiencies Q_ext, Q_abs, Q_sca with the albedo and asymmetry parameter g; the power
  absorbed and emitted per gram (equal areas on a log-wavelength axis); the temperature distribution P(T);
  and temperature against radius, which shows where single photons start to matter.
- **Size distribution**: power laws dn/da ∝ a^-q between a_min and a_max for silicate and graphite,
  normalized by the Si or C atoms locked in grains; the extinction curve (per H or per A_V, against 1/λ or
  λ) and albedo, with the Cardelli, Clayton & Mathis (1989) R_V = 3.1 curve as a reference, and the
  infrared spectrum λI_λ/N_H, each split by material or by size; and summary numbers per H nucleus.

The page runs the whole calculation in the browser, including the Mie series for the single grain and the
stochastic heating of every size in the distribution; only the Mie efficiencies on the size grid are
precomputed. PAHs are included as a third material (single grain, neutral or ionized) and as a log-normal
population in the size distribution, with an ionized-fraction slider and an "MRN + PAHs" preset.

## Physics and assumptions

- **Materials**: homogeneous spheres of astronomical silicate (MgFeSiO4, 3.5 g cm^-3) and graphite
  (2.24 g cm^-3, N_C = 468 (a/1 nm)^3), with the Draine (2003) dielectric functions from Draine's web site.
  Graphite uses the 1/3–2/3 approximation (Draine & Malhotra 1993); its free-electron damping depends on
  radius, and Draine supplies a = 0.01 and 0.10 µm files, which are interpolated linearly in log a and held
  fixed outside that range.
- **Optics**: Mie theory, the BHMIE algorithm of Bohren & Huffman (1983), on 404 wavelengths from the Lyman
  limit to 1 cm (80 per decade); sizes 0.001–10 µm at 25 per decade for the distribution tables.
- **Radiation field**: u_λ = U [u_UV + Σ W_i (4π/c) B_λ(T_i)] + (4π/c) B_λ(T_CMB) (the form of Li & Draine
  2001, Eq. 6), with the Mathis, Mezger & Panagia (1983) blackbodies (7500, 4000, 3000 K; W = 1e-14, 1e-13,
  4e-13; their Table A1) and the UV component taken as their Table A3 (D_G = 10 kpc) minus those blackbodies,
  which is how MMP83 define it. T_CMB = 2.7255 K. The CMB is not scaled with U.
- **Heat capacities**: the Debye models of Draine & Li (2001): graphite (N_C − 2)[kΘ f_2(T/863 K) +
  2kΘ f_2(T/2504 K)], silicate (N_a − 2)[2kΘ f_2(T/500 K) + kΘ f_3(T/1500 K)].
- **Temperatures**: equilibrium from ∫C_abs c u_λ dλ = ∫C_abs 4π B_λ(T) dλ. A grain whose estimated
  temperature spread (an Ornstein–Uhlenbeck estimate from photon statistics) exceeds 1% of T_eq is treated
  stochastically: the thermal-continuous approximation of Draine & Li (2001, Sect. 7.2) with their
  finite-bin-width upward rates, solved by the Guhathakurta & Draine (1989) recursion, emission from the
  thermal approximation with photons limited to the grain's energy.
- **PAHs** (Draine & Li 2007, DL07): cross sections per C atom from Li & Draine (2001) Eqs. 5–11 with the
  DL07 Table 1 bands and their near-infrared continuum and bands for cations; N_C = 460 (a/1 nm)³; H/C
  from DL07 Eq. 3; a 1% graphite admixture below 50 Å, rising to pure graphite above (DL07 Eqs. 4–5); heat
  capacity with N_H Einstein oscillators for the C–H modes plus the Debye C–C modes (DL01); the smallest
  PAH has 20 C atoms (3.55 Å). Computed on a finer wavelength grid (500 per decade over 3–20 µm), which
  also carries the size-distribution sums. The size distribution is one log-normal in radius (center,
  width, C atoms in PAHs); neutral and ionized PAHs are computed separately and mixed by the ionized
  fraction. The preset takes DL07's first log-normal component (4 Å, σ = 0.4) with 45 ppm of C moved out
  of graphite; it gives q_PAH = 5.8% (DL07's Milky Way model: 4.6%).
- **Ground state**: the energy grid never goes below the grain's lowest vibrational mode (from DL01's model
  mode spectra); that bottom bin is the ground state, radiates nothing, and is reported as the fraction
  of time spent there. The temperature panels show the temperatures that radiate 90% of the energy.
- **Reference extinction curve**: CCM89 Eqs. (1)–(5) for R_V = 3.1, read from the ADS scan and identical to
  `dust_extinction`'s implementation to 3e-15; defined for 0.3 ≤ 1/λ ≤ 10 µm⁻¹. Per H it is scaled by
  A_V/N_H = 3.1 / 5.8e21 mag cm², the gas-to-reddening ratio of Bohlin, Savage & Drake (1978, abstract).
- **Left out**: non-spherical grains, ices, coatings and porosity, grain charging (the PAH ionized fraction
  is set by hand), sublimation (temperatures are capped at 4900 K), and the infrared part of the field.

## Things worth knowing (found while building this)

1. **DL01's f_n has the wrong prefactor as printed.** Draine & Li (2001) and Li & Draine (2001) both write
   f_n(x) = (1/n)∫_0^1 y^n dy/(e^{y/x} − 1). The n-dimensional Debye model with their mode counts needs a
   prefactor **n**, and only n gives equipartition (C → 3Nk at high T) and reproduces their Fig. 2.
   validate.py checks the silicate curve of that figure to 0.3%.
2. **DL01 Fig. 2 draws graphite with Θ_op ≈ 950 K** (the Krumhansl & Brooks value quoted in their text),
   while their equations and the LD01 emission models use 863 K. 863 K is used here; it matches the lab
   points in the same figure at 20–30 K to 3%.
3. **Two misprints in MMP83 Table A3** (10 kpc column), each a factor of 10: 0.21 µm reads 9.61E-2 and 3.4 µm
   reads 6.49E-3. Corrected in `data/mmp83_tableA3.txt`; the blackbody sum confirms both (ratio 0.991–0.995
   over 0.35–8 µm).
4. **Graphite runs ~9% warmer than in LD01** (22.3 vs 20.4 K at 0.01 µm, U = 1) because the Draine (2003)
   graphite absorbs about half as much in the far infrared as the 1993 constants LD01 used. With Draine's
   1993 tables the solver reproduces LD01 Fig. 3 to 0.14 K.
5. **P(E) peaks sit 7–18% higher in energy than LD01 Fig. 5**, with heights matching to 3%. That is what a
   silicate grain with ~12% more atoms per volume would show; LD01 state only "3.5 g cm^-3 for MgFeSiO4",
   which is what is used here. Not resolved; it moves the onset of stochastic heating by ~4% in radius.
6. **Numerical diffusion in the continuous-cooling scheme.** Each cooling step to the next bin adds variance
   of about one bin width times the cooling flux, which broadens P(E) whenever bins are not small compared
   with the photon energy. The solver therefore runs twice: a generous first range, then a trimmed range
   with bins of 3% of the mean absorbed photon energy near the peak (400–1500 bins). The transition matrix
   is banded (no photon exceeds 13.6 eV), so the fine binning is cheap in the browser.

7. **CCM89's Table 3 does not match its own equations at B, H and K** (B: 1.337 printed vs 1.323 from Eq. 3;
   H 0.190 vs 0.192; K 0.114 vs 0.116), although the footnote says the table comes from Eqs. (2)–(3). U, V,
   R, I and J agree. The equations are used here, as in `dust_extinction`.

8. **DL07's cation cross section goes negative near 1.9 µm.** Their negative 1.905 µm "resonance"
   (σ_int = −146.5e-20 cm/C, meant to suppress absorption at 1.8–2.0 µm) outweighs the continuum and the
   1.05 µm band wing even with the 1% graphite admixture, leaving C_abs ≈ −1.6e-22 cm²/C. The total is
   floored at zero; it shows as a small notch at 1.9 µm.
9. **DL07's text and Table 1 disagree on several band widths** (6.2 µm γ = 0.0284 in the text, 0.030 in the
   table; also 11.99 and 12.7 µm), and Table 1 row 29 (18.92 µm) reads γ = 0.10, σ = 0.10 (neutral), 0.17
   (ionized) after a commented-out older row. The table is used as printed. Against DL07 Fig. 3a the
   implementation is 1–24% high at band peaks; against their Fig. 6 (band emission fractions vs size, 205
   digitized points) the median ratio is 1.07, by band 1.00 (7.7, 8.6 µm) to 1.17 (11.3 µm).
10. **Draine's tabulated PAH cross sections (PAHneu_30, PAHion_30) follow LD01, but with broader 11.3 and
    12.7 µm bands than LD01 Table 1.** With LD01's own parameters the code matches the tables to 5% except
    in 10–13.5 µm, where the band integral agrees to 7% and only the profile widths differ.
11. **The smallest grains need a ground state.** Without it the thermal-continuous scheme lets a 5 Å PAH
    linger at 35–50 K between photons; with the energy floor at the lowest vibrational mode it spends 98.9%
    of the time in the ground state. For silicate this reproduces LD01's P_0 at 5 and 10 Å (1.000 vs 0.999,
    0.991 vs 0.994); at 25 Å it gives 0.09 against their exact 0.23.

## Files

| File | Purpose |
|---|---|
| `dustlib.py` | the physics, reference implementation (documents every source) |
| `web/physics.js` | the browser port of dustlib (checked against it to ~1e-14) |
| `web/app.js`, `web/index.html`, `web/style.css` | the page |
| `exportweb.py` | writes `web/data/dust.json` + `dust.bin` (Q tables, dielectric functions, field, enthalpies) |
| `validate.py` | 39 checks against the sources and of the port; `python3 validate.py` (~1 min) |
| `testjs.js` | runs physics.js under node for validate.py |
| `uitest.js` | drives the page in headless Chrome (sliders, tabs, hovers, clicks, hash, phone width) |
| `fetch.sh` | downloads every primary source into `rawdata/` |
| `data/mmp83_tableA3.txt` | the transcribed MMP83 table (with the two corrections) |
| `data/dl07_fig6.json` | DL07 Fig. 6 band emission fractions, digitized for validate.py |
| `deploy.sh` | publishes `web/` to the gh-pages branch |

## Running

```
./fetch.sh                      # once, into rawdata/
python3 exportweb.py            # builds web/data/
python3 web/serve.py 8781       # http://localhost:8781/ (stamps scripts so webviews cannot cache them)
python3 validate.py
node uitest.js                  # with the server running
```

## Sources

Bohren & Huffman (1983), *Absorption and Scattering of Light by Small Particles*, App. A ·
Draine & Lee (1984, ApJ 285, 89) · Draine & Malhotra (1993, ApJ 414, 632) · Draine (2003, ApJ 598, 1026) ·
Draine & Li (2001, ApJ 551, 807) · Li & Draine (2001, ApJ 554, 778) · Guhathakurta & Draine (1989, ApJ 345,
230) · Mathis, Mezger & Panagia (1983, A&A 128, 212) · Cardelli, Clayton & Mathis (1989, ApJ 345, 245) ·
Bohlin, Savage & Drake (1978, ApJ 224, 132) · Mathis, Rumpl & Nordsieck (1977, ApJ 217, 425) ·
Fixsen (2009, ApJ 707, 916) · Draine & Li (2007, ApJ 657, 810) · Mattioda, Allamandola & Hudgins (2005).
