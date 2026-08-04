"""Constants and cut logic of the MAIA calorimeter digitization chain.

Everything the emulation needs to reproduce the chain lives here, with the
source of each number cited. Chain (per calorimeter region):

  SimCalorimeterHit contributions
    -> RealisticCaloDigi{Silicon,ScinPpd}   (timing window on delta_t, threshold)
    -> RealisticCaloReco{Silicon,ScinPpd}   (calibration only, keeps time)
    -> CaloConer                            (truth cone, no-op for guns)
    -> CaloHitSelector                      (threshold + time window on stored time)
    -> Pandora                              (no timing cuts)

All parameter values from MAIAConfig/MAIAConfig/CaloDigi/{calorimetry_EM.py,
calorimetry_HAD.py, calo_coning.py} @ main; algorithm behavior from
MuonColliderSoft/k4Reco @ main (see ../PLAN.md section 1).
"""

import numpy as np

C_MM_NS = 299.792458  # speed of light [mm/ns] (CLHEP c_light)

# Region indices used across the flat arrays
REGIONS = ["EcalBarrel", "EcalEndcap", "HcalBarrel", "HcalEndcap"]
SIM_COLLECTIONS = {
    "EcalBarrel": "ECalBarrelCollection",
    "EcalEndcap": "ECalEndcapCollection",
    "HcalBarrel": "HCalBarrelCollection",
    "HcalEndcap": "HCalEndcapCollection",
}
ECAL = (0, 1)
HCAL = (2, 3)

# --- digitizer settings (calorimetry_EM.py / calorimetry_HAD.py) ------------
DIGI_WINDOW_MIN = -0.5   # ns, timingWindowMin (all regions)
DIGI_WINDOW_MAX = 10.0   # ns, timingWindowMax (all regions)

CALIB_MIP = {            # calibration_mip [GeV deposit per MIP]
    "EcalBarrel": 0.0001575,
    "EcalEndcap": 0.0001575,
    "HcalBarrel": 0.0004925,
    "HcalEndcap": 0.0004725,
}
# digi threshold: ECAL "5e-5 GeV", HCAL "0.5 MIP" -> expressed in MIPs
DIGI_THRESHOLD_MIP = {
    "EcalBarrel": 5e-5 / CALIB_MIP["EcalBarrel"],   # ~0.317 MIP
    "EcalEndcap": 5e-5 / CALIB_MIP["EcalEndcap"],
    "HcalBarrel": 0.5,
    "HcalEndcap": 0.5,
}
# HCAL SiPM (RealisticCaloDigiScinPpd): ppd_mipPe, ppd_npix
PPD_PE_PER_MIP = 15.0
PPD_N_PIXELS = 2000.0

# --- reco calibration (RealisticCaloReco*, single layer group) --------------
CALIB_MIP_TO_GEV = {     # calibration_factorsMipGev
    "EcalBarrel": 0.00641222630095,
    "EcalEndcap": 0.00641222630095,
    "HcalBarrel": 0.0287783798145,
    "HcalEndcap": 0.0285819096797,
}

# --- CaloHitSelector settings (calo_coning.py) -------------------------------
SEL_WINDOW_MIN = -0.3    # ns, exclusive
SEL_WINDOW_MAX = 0.3     # ns, exclusive
SEL_FLAT_THRESHOLD_HCAL = 5e-5  # GeV, on the calibrated (rec) hit energy
# ECAL uses the per-(theta, layer) mode map from ECAL_Thresholds_10TeV.root
# (th_2dmode_sym, Nsigma = 0), extracted by production/extract_thresholds.sh.


def digi_cell_response_mip(e_gev_windowed, region):
    """Windowed deposited energy [GeV] -> digitized cell energy in MIPs
    (deterministic mean of RealisticCaloDigi{Silicon,ScinPpd}; the silicon
    Poisson e-h and scintillator binomial pixel smearing are stochastic and
    intentionally omitted -> validate with tolerances)."""
    e_mip = np.asarray(e_gev_windowed, dtype=np.float64) / CALIB_MIP[region]
    if region in ("HcalBarrel", "HcalEndcap"):
        # mean SiPM saturation in npe, expressed back in (saturated) MIPs
        npe = e_mip * PPD_PE_PER_MIP
        npe = PPD_N_PIXELS * (1.0 - np.exp(-npe / PPD_N_PIXELS))
        e_mip = npe / PPD_PE_PER_MIP
    return e_mip


def digi_threshold_pass(e_mip_digitized, region):
    """RealisticCaloDigi threshold: energyDig > threshold (converted to the
    digitizer's internal unit, i.e. MIP for silicon / NPE for scint; both
    reduce to a MIP comparison here). Applied to the saturated energy."""
    return e_mip_digitized > DIGI_THRESHOLD_MIP[region]


def rec_energy_gev(e_mip_digitized, region):
    """Digitized (saturated) cell energy -> calibrated rec energy [GeV].
    RealisticCaloRecoScinPpd de-saturates with the exact inverse of the mean
    saturation, so for the deterministic emulation the HCAL round trip is the
    identity below the linearization point (0.95 * npix ~ 133 MIP)."""
    e_mip = np.asarray(e_mip_digitized, dtype=np.float64)
    if region in ("HcalBarrel", "HcalEndcap"):
        npe = e_mip * PPD_PE_PER_MIP
        r = 0.95
        lin = npe < r * PPD_N_PIXELS
        npe = np.where(
            lin,
            -PPD_N_PIXELS * np.log(np.clip(1.0 - npe / PPD_N_PIXELS, 1e-12, None)),
            (npe - r * PPD_N_PIXELS) / (1.0 - r) - PPD_N_PIXELS * np.log(1.0 - r),
        )
        e_mip = npe / PPD_PE_PER_MIP
    return e_mip * CALIB_MIP_TO_GEV[region]


class SelectorThresholdMap:
    """Replicates CaloHitSelector's per-(theta, layer) threshold lookup.

    threshold = mode(theta_sym, layer) + Nsigma * stddev  (Nsigma = 0 here);
    theta is folded around pi/2 (th_2dmode_sym binning); lookup mirrors
    TH2::FindBin including under/overflow clamping.
    """

    def __init__(self, npz_path, det="ECAL"):
        d = np.load(npz_path)
        self.xedges = d[f"{det}_th_2dmode_sym_xedges"]
        self.yedges = d[f"{det}_th_2dmode_sym_yedges"]
        self.values = d[f"{det}_th_2dmode_sym_values"]  # incl. under/overflow

    def threshold(self, theta, layer):
        theta = np.asarray(theta, dtype=np.float64)
        theta = np.where(theta > np.pi / 2, np.pi - theta, theta)
        ix = np.searchsorted(self.xedges, theta, side="right")  # 0..nx+1 = ROOT bin
        iy = np.searchsorted(self.yedges, np.asarray(layer, dtype=np.float64),
                             side="right")
        return self.values[ix, iy]


def selector_pass(e_rec_gev, t_cell, theta, layer, region,
                  window_min=SEL_WINDOW_MIN, window_max=SEL_WINDOW_MAX,
                  ecal_map=None):
    """Full CaloHitSelector decision for already-digitized cells.

    t_cell is the stored digi hit time = earliest accepted contribution
    delta_t. (The selector subtracts another r/TMath::C(), but with r in mm
    and C in m/s that is ~1e-6 ns; treated as exactly zero here - see PLAN.md.)
    Both cuts exclusive, threshold first, mirroring CaloHitSelector.cpp.
    """
    if region in ("HcalBarrel", "HcalEndcap"):
        thr = SEL_FLAT_THRESHOLD_HCAL
    else:
        if ecal_map is None:
            raise ValueError("ECAL selector needs the extracted threshold map")
        thr = ecal_map.threshold(theta, layer)
    ok_e = np.asarray(e_rec_gev) > thr
    ok_t = (t_cell > window_min) & (t_cell < window_max)
    return ok_e & ok_t


class BitFieldCoder:
    """Minimal pure-python DD4hep BitFieldCoder (decode only).

    Parses encodings like "system:5,side:-2,layer:9,module:8,..." where a
    leading '-' on the width marks a signed field and an optional explicit
    start bit can be given as name:start:width.
    """

    def __init__(self, encoding):
        self.fields = {}
        pos = 0
        for tok in encoding.split(","):
            parts = tok.split(":")
            name = parts[0]
            if len(parts) == 2:
                start, width = pos, int(parts[1])
            else:
                start, width = int(parts[1]), int(parts[2])
            signed = width < 0
            width = abs(width)
            self.fields[name] = (start, width, signed)
            pos = start + width

    def get(self, cellid, name):
        start, width, signed = self.fields[name]
        val = (np.asarray(cellid).astype(np.uint64) >> np.uint64(start)) \
            & np.uint64((1 << width) - 1)
        val = val.astype(np.int64)
        if signed:
            sign_bit = np.int64(1 << (width - 1))
            val = np.where(val >= sign_bit, val - np.int64(1 << width), val)
        return val
