# MAIA calorimeter timing-cut study

**Goal:** Understand how the timing cuts applied along the MAIA reconstruction chain
remove hits from calorimeter clusters, in the ECAL and the HCAL, and arrive at a
suggested set of timing cuts that does not strongly shape the cluster energy.

**Setup:** the tutorial-tagged Gaudi chain
([mcd-wiki tutorial](https://mcd-wiki.web.cern.ch/software/howto/gaudi/)):
container `ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v3.0`, steering from
[mucoll-benchmarks](https://github.com/MuonColliderSoft/mucoll-benchmarks) with the
[MAIAConfig](https://github.com/MuonColliderSoft/MAIAConfig) submodule (`main`),
algorithms from [k4Reco](https://github.com/MuonColliderSoft/k4Reco). Geometry `MAIA_v0`.
**No BIB overlay.**

All times in this study are *relative*:
`delta_t = t_arrival - r_cell/c`, the difference between the actual arrival time of a
deposit and the time a speed-of-light particle from the IP would need to reach that
detector cell. This is the convention the reconstruction code itself uses (see
inventory below).

**Samples** (per user spec): 100 events per sample, guns at theta = 90 deg exactly
(barrel-only study; endcap + theta-dependence deferred to a follow-up).

- Phase 1 (ECAL): photon guns, E = 10, 50, 200 GeV.
- Phase 2 (ECAL+HCAL combined): pi- guns, E = 10, 50, 200 GeV.

---

## 1. Where the chain cuts on time (verified against source)

The chain applies timing requirements at **two** places; a third stage (Pandora) applies
none. The digitization stage lumps contributions together, which is why this study works
at the simhit level (each Geant4 contribution carries its own time).

| # | Stage | Algorithm (instance names) | Cut variable | Window / value | Acts on |
|---|-------|----------------------------|--------------|----------------|---------|
| 1 | digi | `RealisticCaloDigiSilicon` (`ECALBarrelDigi`, `ECALEndcapDigi`), `RealisticCaloDigiScinPpd` (`HCALBarrelDigi`, `HCALEndcapDigi`) | per-contribution `delta_t = t_contrib - r_cell/c` | `timingWindowMin/Max = (-0.5, +10) ns` | each SimCalorimeterHit contribution |
| 1b | digi | same | energy threshold on windowed cell energy | ECAL `5e-5 GeV` (~0.32 MIP); HCAL `0.5 MIP` | each cell |
| 2 | digi (selection) | `CaloHitSelector` (`EcalBarrelSelector`, ..., `HcalEndcapSelector`) | stored digi hit time (see note) | `TimeWindowMin/Max = (-0.3, +0.3) ns`, exclusive | each digitized cell |
| 2b | digi (selection) | same | energy threshold | ECAL: per-(theta, layer) BIB mode map (`ECAL_Thresholds_10TeV.root:th_2dmode_sym`, `Nsigma=0`); HCAL: flat `5e-5 GeV` | each cell |
| 3 | reco | `DDPandoraPFANewAlgorithm` | — | **no timing cuts** | — |

Code references (k4Reco @ main, MAIAConfig @ main):

- `k4Reco/CaloDigi/src/RealisticCaloDigi.cc`, `StandardIntegration()`:
  `relativetime = t_i - sqrt(x^2+y^2+z^2)/c` (CLHEP `c_light`, mm/ns); contributions with
  `timingWindowMin < relativetime < timingWindowMax` are summed;
  **the stored digi hit time is the EARLIEST accepted relative time** (already
  propagation-corrected); `timingResolution = 0` (no smearing). Both calorimeters use
  the `Standard` integration method (the `ROC` fast/slow-shaper option is not used).
- `MAIAConfig/CaloDigi/calorimetry_EM.py`, `calorimetry_HAD.py`: the windows and
  thresholds quoted above; `timingCorrectForPropagation = 1`, `timingCut = 1`.
- `MAIAConfig/CaloDigi/calo_coning.py`: `CaloConer` (0.6 rad truth cone around MC
  particles — a no-op for single-gun no-BIB samples) then `CaloHitSelector`.
- `k4Reco/BIBUtils/src/CaloHitSelector.cpp`: threshold first, then
  `relativeTime = hit.getTime() - r/TMath::C()`, cut `(-0.3, +0.3)` exclusive.

### Important subtleties (the "used in conjunction" part)

1. **The selector cuts on the earliest-contribution time.** The digi hit time is the
   earliest accepted contribution's delta_t, while the digi hit *energy* integrates
   everything out to +10 ns. So the +/-0.3 ns selector removes *whole cells* whose first
   in-window deposit is late (typical for neutron-dominated or shower-tail cells), but
   keeps every late deposit that shares a cell with a prompt one. A tighter digi window
   and a tighter selector window therefore remove different energy.
2. **A would-be double time-of-flight correction.** The digi hit time is already
   relative, yet `CaloHitSelector` subtracts `r/TMath::C()` again. It is saved by a
   units accident: `r` is in mm and `TMath::C()` in m/s, so the second correction is
   ~1e-6 ns and numerically irrelevant. Functionally the selector cuts on the digi's
   already-corrected time — which is the sensible behavior — but the code only does so
   by accident. **Worth reporting upstream** (if someone "fixes" the units it becomes a
   real double correction).
3. **Digi threshold couples to the timing window.** The cell threshold (1b) is applied to
   the *windowed* energy, so tightening the digi window also kills cells via the
   threshold, not just by removing their late energy.
4. **Stochastic digitization** (affects validation only, not the physics conclusions):
   ECAL silicon applies a Poisson e-h-pair fluctuation (`silicon_pairEnergy = 3.6 eV`,
   sub-% effect); HCAL scintillator applies binomial SiPM pixel smearing
   (2000 px, 15 pe/MIP, ~26% at 1 MIP). Per-cell closure checks against the real chain
   are therefore tolerance-based; the emulation itself uses the deterministic mean.
5. **The ECAL selector threshold is large even without BIB.** The per-(theta, layer)
   mode map is ~6.5 MeV (rec scale, ~1 MIP) per cell at theta = 90 - a *BIB*-derived
   threshold that also bites on pure signal. This study keeps it in both numerator and
   denominator so the *timing* effect is isolated, but it is worth noting when
   interpreting absolute cluster energies.
6. **Found while validating (upstream k4Reco bugs, independent of timing):**
   `RealisticCaloReco.cc:72` does `int cellid = hit->getCellID()`, truncating the
   64-bit cellID to 32 bits: every Rec/Coned/Sel calorimeter hit loses its `x`/`y`
   fields (bits 32-63). Layer/stave/module survive (bits < 32), so calibration and the
   selector still work, but any downstream consumer of full Sel-hit cellIDs gets
   degenerate IDs. Together with the units accident in note 2, both are worth an
   upstream report.

### Environment quirks pinned down during setup (documented for reproducibility)

- MAIAConfig `main` requires a newer `k4ActsTracking` than the v3.0 container ships
  (`CKF_Chi2CutOffOutlier`), i.e. the tutorial's "compatibility issues" warning. We pin
  the MAIAConfig submodule to `d290d02` (2026-06-17, contemporary with the image;
  identical calo settings to main) in `config.sh` / `setup_benchmarks.sh`.
- The image entrypoint launches an interactive shell (swallowing commands), so
  `production/run_container.sh` bypasses it and sources `/opt/setup_mucoll.sh` itself.
- The selector threshold maps ship with `MyBIBUtils` in this image; the wrapper exports
  `MUCOLL_CALO_THRESHOLDS_DIR` accordingly.
- `steer_baseline.py` hardcodes an input file; `production/sim_gun_steer.py` wraps it
  and clears `SIM.inputFiles` so `--enableGun` works.

---

## 2. Method

**Primary: simhit-level emulation.** From the sim-level `SimCalorimeterHit`
contributions (saved by `steer_baseline.py` via `enableDetailedShowerMode = True`) we
recompute the full digi + selector chain in numpy, with every timing window a free
parameter:

1. per contribution: `delta_t = t - |cell position|/c`;
2. digi window scan -> windowed cell energy + earliest accepted time;
3. digi threshold (ECAL: GeV-deposit; HCAL: MIP with mean SiPM saturation);
4. reco calibration (MIP -> GeV factors from `calorimetry_*.py`);
5. selector threshold (extracted ECAL map / flat HCAL) + selector window scan;
6. per event: retained energy, at every (digi window x selector window) grid point.

Key computational note: with the digi window *minimum* fixed, the earliest accepted time
per cell does not depend on the digi window *maximum* — the max edge only changes the
cell's energy and survival. This makes the 2D scan cheap and the coupling exact.

**Validation: the real chain.** Re-run `k4run digi_steer.py` + `reco_steer.py` on the
same sim files with property overrides (e.g. `--EcalBarrelSelector.TimeWindowMax 0.15`)
at a few points, and compare (a) summed `*CollectionSel` energy against the emulation,
(b) the leading Pandora cluster energy against the Sel sum (quantifies what clustering
adds on top of hit selection).

### Scan grid

- digi `timingWindowMax`: {1, 2, 3, 5, 10*, 15, 25} ns; extended with {50, 100} for HCAL.
- digi `timingWindowMin`: {-0.25, -0.5*, -1} ns (cross-check).
- selector symmetric half-width w: {0.05, 0.1, 0.15, 0.2, 0.25, 0.3*, 0.4, 0.5, 0.75, 1.0} ns.
- selector asymmetric: min = -0.3 fixed, max in {0.1 ... 5} ns.

(* = current defaults.)

### Metrics

- Energy-weighted delta_t spectra, at contribution level and at cell-earliest level,
  per subdetector / gun energy / layer group -> where the late tail lives.
- Retained Sel-energy fraction vs cut value (curves per gun energy).
- Response (median retained energy / no-timing-cut energy) and resolution (68% IQR /
  2 x median) vs cut value.
- Recommendation criterion: the tightest windows keeping the median energy loss below
  ~0.5-1% at all gun energies, with no residual energy-dependence (shaping) beyond
  statistics.

---

## 3. Workflow

```
config.sh                        # every knob in one place (image, grid, paths)
production/setup_benchmarks.sh   # clone mucoll-benchmarks + MAIAConfig (pinned commits)
production/run_container.sh CMD  # docker (or apptainer) wrapper
production/run_chain.sh          # gun -> ddsim -> digi -> reco for one grid point
production/launch_grid.sh        # all grid points, parallel, writes manifest
production/extract_thresholds.sh # ECAL selector threshold map -> data/thresholds_ecal.npz
analysis/flatten.py              # (in container) EDM4hep -> npz arrays
analysis/sanity_check.py         # closure gates (must pass before scanning)
analysis/emulate.py              # the scan engine
analysis/plot_distributions.py   # delta_t spectra
analysis/plot_scans.py           # retained fraction / response / resolution vs cuts
production/run_digireco_variant.sh + analysis/summarize_reco.py + analysis/validate.py
```

Data files stay in `timing_studies/data/` (gitignored, reproducible from seeds via the
manifest); final plots in `timing_studies/plots/`; numbers + conclusions get written
back into this file.

### Verification gates (run before any scan is trusted)

1. flatten self-consistency: sum of contribution energies == SimCalorimeterHit energy;
2. prompt-photon delta_t peaks at ~0 (no unit/offset errors);
3. per-cell join of emulated vs actual Digi hits at default cuts (time exact,
   energy within stochastic tolerance, matching cell sets);
4. same vs Sel collection (validates threshold-map lookup + windows in one shot);
5. event-level emulated vs actual Sel sums;
6. one real-chain variant run reproduced by the emulation.

## 4. Status / decision log

- [x] Cut inventory verified against k4Reco/MAIAConfig source (this file, section 1).
- [x] Phase 0: scaffolding, container smoke test, threshold-map extraction.
- [x] Phase 1: photon samples (E = 10/50/200 GeV, 100 evts, theta = 90); ECAL results.
- [x] Phase 2: pion samples (same grid); ECAL+HCAL combined results.
- [x] Suggested timing cuts written here (ECAL: section 5.3; combined: section 5.5).

## 5. Results

### 5.1 Phase 1 closure (photons, theta = 90, no BIB)

All verification gates pass. Highlights (100 events/point):

- Per-cell emulated vs actual **Digi** hits (ECAL barrel): cell sets match to
  <0.05%, energies to 0.2% (the Poisson e-h smear), times to <1e-6 ns.
- Emulated vs actual **Sel** event sums (ECAL barrel): median deviation 0.01-0.04%.
- **Real-chain variant runs** (selector +/-0.15 ns, +/-1.0 ns; digi max 2 ns, at all
  three energies): emulation matches the actual Sel sums with median deviation
  < 0.01% and 68% width < 0.06% -> the scan below is equivalent to re-running the
  chain at every grid point.
- **Leading Pandora cluster / Sel-hit sum = 1.021-1.023** at every variant point,
  i.e. exactly Pandora's `ECalToEMGeVCalibration = 1.0237`: for single photons the
  cluster energy tracks the hit-selection energy perfectly, so hit-level retained
  fractions ARE cluster-level retained fractions.
- HCAL closure at this occupancy (photon leakage, ~2 cells/event near threshold) is
  limited by the un-emulated binomial SiPM smear, as expected; at >= 20 cells/event
  it closes to <1%. Properly exercised in phase 2.

### 5.2 Phase 1 ECAL results (plots in `plots/ecal_photons/`)

Energy-weighted delta_t spectra (`dt_*_EcalBarrel.png`): a prompt peak at
~0.02 ns containing ~99.5% of the energy within +/-0.5 ns, falling by three
orders of magnitude within ~1 ns, then a flat few-1e-5/ns tail out past 10 ns.

Median retained ECAL-barrel energy fraction (denominator: same thresholds, no
timing cuts; digi window at default (-0.5, 10) ns):

| selector window [ns] | 10 GeV | 50 GeV | 200 GeV |
|---------------------|--------|--------|---------|
| abs(t) < 0.2         | 0.9860 | 0.9860 | 0.9880 |
| abs(t) < 0.25        | 0.9918 | 0.9911 | 0.9917 |
| **abs(t) < 0.3 (current)** | **0.9942** | **0.9940** | **0.9941** |
| abs(t) < 0.4         | 0.9966 | 0.9968 | 0.9969 |
| abs(t) < 0.5         | 0.9983 | 0.9979 | 0.9982 |
| abs(t) < 1.0         | 1.0000 | 0.9990 | 0.9991 |
| -0.3 < t < 0.5       | 0.9983 | 0.9979 | 0.9982 |
| -0.3 < t < 1.0       | 1.0000 | 0.9990 | 0.9991 |

Event-to-event 68% spread of the retained fraction (a direct resolution
contribution): 0.46% / 0.21% / 0.13% (10/50/200 GeV) at +/-0.3 ns, halving at
+/-0.5 ns and halving again at (-0.3, +1.0) ns.

Key observations:

1. **The digi window is a non-issue for the ECAL**: even `timingWindowMax = 1 ns`
   keeps > 99.8% of the energy; at 2 ns the loss is < 0.1% at all energies
   (`retained_vs_digimax.png`). The current 10 ns is generously safe - or could be
   *tightened to ~2-3 ns* for free if BIB integration in the +10 ns tail proves
   harmful.
2. **All the timing loss is on the late side**: the asymmetric rows equal the
   symmetric rows with the same max. The early edge (-0.3 ns) removes nothing from
   photon clusters - backward smearing of the prompt peak is < 1e-4.
3. **The current +/-0.3 ns selector cut removes a flat 0.6%** of cluster energy,
   *independent of energy* between 10 and 200 GeV - a calibratable offset, not a
   shape, with a small resolution penalty (<= 0.5%). Strong energy-dependent shaping
   only starts below ~0.2 ns half-width.

### 5.3 Suggested ECAL timing cuts (photons, no BIB - to be revisited with BIB)

- The **selector late edge is the only knob that matters**. `TimeWindowMax = 0.3 ns`
  (current) is acceptable (0.6% flat loss); raising it to **0.5 ns** cuts the loss
  to < 0.25% and the resolution penalty to < 0.3%; **+1.0 ns** makes the cut
  effectively lossless (<= 0.1%).
- Keep the early edge tight: `TimeWindowMin = -0.3 ns` costs nothing.
- The digi window can stay at (-0.5, +10) ns; if BIB pileup in the integration
  window is a concern, it can be tightened to `timingWindowMax ~ 2-3 ns` with
  < 0.1% signal cost.

(How much BIB the looser selector edge re-admits is outside this no-BIB study and
must be checked against the BIB occupancy before adopting +0.5/+1.0 ns.)

### 5.4 Phase 2 results (pi-, theta = 90, no BIB; plots in `plots/hcal_pions/`)

Closure at hadronic occupancy: all gates pass. ECAL as in phase 1 (<0.1%); HCAL
event Sel sums close to 1.5% / 1.0% / 0.3% (10/50/200 GeV) - the residual is
entirely the un-emulated binomial SiPM smear, and the real-chain variant runs
agree with the emulation to 0.1-0.8% median (68% width 0.5-4%).

The delta_t spectra show why hadrons are different: only ~59% (HCAL) and ~78%
(ECAL) of the pion sim energy arrives within +/-0.5 ns; the neutron-dominated
tail extends over tens (HCAL: hundreds) of ns.

Median retained energy fraction, **ECAL+HCAL barrel summed** (68% event-to-event
spread in parentheses; denominator: same thresholds, no timing cuts):

| cut | 10 GeV | 50 GeV | 200 GeV |
|---|---|---|---|
| **current (digi 10 ns, sel +/-0.3 ns)** | **0.709 (0.129)** | **0.789 (0.066)** | **0.845 (0.046)** |
| sel +/-0.5 ns | 0.799 (0.116) | 0.852 (0.048) | 0.889 (0.036) |
| sel +/-1.0 ns | 0.883 (0.071) | 0.908 (0.034) | 0.929 (0.022) |
| sel +/-2.0 ns | 0.923 (0.051) | 0.938 (0.024) | 0.952 (0.017) |
| sel +/-5.0 ns | 0.965 (0.025) | 0.969 (0.013) | 0.975 (0.009) |
| no selector (digi 10 ns alone) | 0.983 (0.013) | 0.987 (0.007) | 0.987 (0.005) |
| digi 25 ns, no selector | 0.994 (0.007) | 0.996 (0.003) | 0.996 (0.002) |

HCAL-only at the current cut retains just 0.52 / 0.61 / 0.76 with 37% / 16% / 10%
event-to-event spread; ECAL-only for pions retains 0.79 / 0.87 / 0.92 (hadronic
interactions in the ECAL also produce late energy). Asymmetric windows again equal
symmetric ones: all loss is on the late side.

Key observations:

1. **The +/-0.3 ns selector cut strongly shapes hadronic cluster energy**: it
   removes 15-29% of the calorimeter energy with an ~14% response non-linearity
   between 10 and 200 GeV, and adds a 5-13% event-to-event spread - a major
   resolution term. Cross-checked in the real chain: the median 10 GeV pion
   leading-cluster energy moves from 5.35 GeV (selector +/-0.15 ns) to 8.90 GeV
   (+/-1.0 ns).
2. **No reasonable window recovers everything**: the late tail is intrinsic to
   hadronic showers, so the window choice is a genuine trade-off against BIB.
   Doubling the late edge roughly halves both the energy loss and the added spread
   (+1 ns -> ~90% retained; +2 ns -> ~94%; +5 ns -> ~97%).
3. **The digi window matters for hadrons too, but mildly**: 10 ns costs 1.3-1.7%
   (energy-independent); 25 ns recovers to <0.6%.
4. The interplay documented in section 1 is visible: with the selector open, energy
   integrated out to 10 ns survives via prompt-earliest cells; the selector then
   discards whole late-first cells, which dominate the hadronic loss.
5. The 10 GeV HCAL 68%-spread curve is non-monotonic (peaks ~0.40 at w ~ 0.4-0.75 ns)
   because the per-event retained-fraction distribution is bimodal: a population of
   events retains exactly f = 0 (17/96 at the current +/-0.3 ns, still 9 at +/-2 ns -
   no cell with an early-enough first deposit) while prompt-rich events sit near
   f = 1. The spread peaks where the median crosses mid-range (p(1-p)-type maximum)
   and collapses once everything converges to f ~ 1; at 50/200 GeV showers average
   over enough cells that the distribution is single-peaked and the spread monotonic.
   See `analysis/plot_fraction_dists.py` ->
   `plots/hcal_pions/fraction_dist_HcalBarrel.png`. Consequence: at 10 GeV and
   +/-0.3 ns the HCAL flips between seeing the shower and seeing nothing event by
   event - the "spread" there is population splitting, not a Gaussian resolution.

### 5.5 Per-subdetector selector-window tuning and suggested cuts (no BIB)

The selector acts per cell, so the ECAL and HCAL windows factorize exactly and can
be tuned independently. Median retained fraction per subdetector (68% spread), with
the pion response linearity ratio retained(200 GeV)/retained(10 GeV) as the shaping
metric:

**ECAL barrel window** (must serve photons AND the ECAL half of hadronic showers):

| window [ns] | gamma 10 | gamma 200 | pi 10 | pi 50 | pi 200 | pi 200/10 |
|---|---|---|---|---|---|---|
| +/-0.3 (current) | 0.994 | 0.994 | 0.794 (0.18) | 0.865 (0.08) | 0.916 (0.05) | 1.153 |
| +/-0.5 | 0.998 | 0.998 | 0.870 (0.13) | 0.911 (0.06) | 0.944 (0.04) | 1.084 |
| +/-1.0 | 1.000 | 0.999 | 0.919 (0.07) | 0.944 (0.04) | 0.964 (0.03) | 1.049 |
| +/-2.0 | 1.000 | 1.000 | 0.959 (0.04) | 0.965 (0.03) | 0.973 (0.02) | 1.016 |
| +/-5.0 | 1.000 | 1.000 | 0.979 (0.02) | 0.984 (0.01) | 0.987 (0.01) | 1.008 |
| open | 1.000 | 1.000 | 0.990 (0.01) | 0.991 (0.01) | 0.992 (0.01) | 1.002 |

**HCAL barrel window** (pions; photons contribute only ~0.05% leakage):

| window [ns] | pi 10 | pi 50 | pi 200 | pi 200/10 |
|---|---|---|---|---|
| +/-0.3 (current) | 0.515 (0.37) | 0.613 (0.16) | 0.762 (0.10) | 1.479 |
| +/-0.5 | 0.651 (0.38) | 0.730 (0.15) | 0.827 (0.08) | 1.270 |
| +/-1.0 | 0.785 (0.31) | 0.839 (0.11) | 0.886 (0.05) | 1.129 |
| +/-2.0 | 0.846 (0.19) | 0.886 (0.07) | 0.919 (0.04) | 1.087 |
| +/-5.0 | 0.934 (0.10) | 0.947 (0.03) | 0.960 (0.02) | 1.028 |
| open | 0.983 (0.04) | 0.978 (0.02) | 0.982 (0.01) | 1.000 |

Suggested cuts (tightest windows whose shaping is small; every value to be
re-weighed against BIB admission):

- **ECAL selector: (-0.3, +2.0) ns.** Photons stay exactly lossless (they already
  are from +0.5 ns); the ECAL part of hadronic showers recovers to 96-97% retained
  with response non-linearity down from 15% to 1.6% and spread <= 4%. If BIB forces
  a tighter edge, +1.0 ns is the compromise (5% non-linearity); for an
  *EM-only* trigger-style selection +0.5 ns suffices.
- **HCAL selector: (-0.3, +5.0) ns**, floor (-0.3, +2.0) ns. At +5 ns the hadronic
  response flattens (non-linearity 2.8%, retained 93-96%, spread <= 10%); at the
  +2 ns floor the non-linearity is 8.7%. The current +/-0.3 ns (48% non-linearity,
  37% spread at 10 GeV) should not be used for hadronic calorimetry without BIB
  pressure that demands it.
- **Early edges: -0.3 ns everywhere** - they remove nothing from signal at any
  energy, for either particle type.
- **Digi windows** (secondary): ECAL can be tightened 10 -> 2-3 ns for free; HCAL
  loosened 10 -> 25 ns recovers 1-2% of hadronic energy (relevant mainly if the
  selector is opened to +5 ns, where the digi edge becomes the limit).
- Any selector retune shifts the hadronic energy scale by 10-25%, so Pandora's
  hadronic calibration (`HCalToHadGeVCalibration`, software compensation weights)
  must be re-derived alongside.

The natural follow-up is repeating this with BIB overlay to place the other side of
the trade-off on the same axes (signal retained vs BIB admitted per window), and a
theta scan to extend the barrel-only tuning to the endcaps.
