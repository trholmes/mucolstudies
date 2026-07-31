# PR title

Fix CaloHitSelector pseudo TOF correction and RealisticCaloReco cellID truncation

# PR body (paste below this line)

This fixes two small but consequential issues found while auditing the calorimeter
timing cuts of the MAIA reconstruction chain (simhit-level timing study, validated
against chain re-runs in the `mucoll-sim-ubuntu24:v3.0` stack).

## 1. `CaloHitSelector`: remove the pseudo time-of-flight correction

`k4Reco/BIBUtils/src/CaloHitSelector.cpp` subtracted `r / TMath::C()` from the hit
time before applying the time window. Two problems compound there:

- The hit time stored by `RealisticCaloDigi` (with `timingCorrectForPropagation`)
  is **already** corrected for the time of flight from the IP — it is the earliest
  accepted contribution's corrected time — so subtracting a TOF again would
  double-count it.
- The subtracted term mixes units: `r` is in mm while `TMath::C()` is in m/s, so
  it actually evaluates to ~7×10⁻⁶ ns instead of the ~6 ns a real TOF would be —
  a numerical no-op that made the double correction invisible.

The two bugs cancel, which is exactly why this is worth fixing: if someone later
"fixes" the units alone, the selector silently starts double-correcting and the
time window moves by the full TOF. This PR removes the term and documents that the
input time is already TOF-corrected, cutting on `hit.getTime()` directly.

## 2. `RealisticCaloReco`: keep the full 64-bit cellID

`k4Reco/CaloDigi/src/RealisticCaloReco.cc` had `int cellid = hit->getCellID()`,
truncating the 64-bit cellID to 32 bits before writing it to the output hit. With
cellID encodings that use the upper word — e.g. the MAIA calorimeter encoding
`system:0:5,side:5:-2,module:7:8,stave:15:4,layer:19:9,submodule:28:4,x:32:-16,y:48:-16`
— every reconstructed (and downstream Coned/Sel) calorimeter hit loses its `x`/`y`
cell indices, and hits from different cells become indistinguishable by cellID.
Fixed by using `std::uint64_t`.

## Validation: outputs are unchanged (except the intended cellID repair)

Rebuilt this branch inside `ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v3.0` and
re-ran the full MAIAConfig digi+reco chain (no BIB, same seeds) on two 100-event
MAIA_v0 samples, comparing every calorimeter hit collection (Digi/Rec/Coned/Sel ×
ECAL/HCAL × barrel/endcap: cellID, energy, time, hit-by-hit) and all Pandora
cluster/PFO energies against the unpatched code:

- **50 GeV photon gun (θ = 90°):** bit-identical everywhere — same hits, energies
  and times to the last bit, identical Pandora clusters/PFOs. Rec/Coned/Sel
  cellIDs now carry the full 64 bits (their 32-bit truncation equals the old
  values, and they match the digitized hits' cellIDs).
- **10 GeV π⁻ gun (θ = 90°):** identical except **2 hits out of ~19,000 selected**,
  each sitting within the removed ~7×10⁻⁶ ns pseudo-term of the +0.3 ns selector
  edge (t = 0.300007 and 0.300003 ns). The old code kept them by accident; the
  fixed code correctly rejects them (t ≥ TimeWindowMax). Knock-on effect: one
  Pandora cluster shifts by 0.17 GeV in each of those two events; everything else
  is bit-identical.

So behavior is preserved to float precision by construction, and the only decision
changes are hits within ~10⁻⁵ ns of the window edge, where the new behavior is the
correct one.

---
🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01QfkMtYniJDh3tYADrTGibX
