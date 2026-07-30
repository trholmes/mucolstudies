// Generates MAIA_timing_cuts_study.pptx — slides for the calorimeter
// timing-cut study (context, method, results, suggested cuts).
// Run: node make_slides.js   (from timing_studies/slides/)
const pptxgen = require("pptxgenjs");

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5 in

// palette
const NAVY = "1E2761";
const ICE = "CADCFC";
const ICE_LIGHT = "EAF1FD";
const ACCENT = "C0392B"; // used for "current cut" highlights
const GRAY = "5A6472";
const WHITE = "FFFFFF";

const W = 13.33, H = 7.5, M = 0.6;
const PLOTS = "../plots";

function titleBar(slide, title, subtitle) {
  slide.addText(title, {
    x: M, y: 0.32, w: W - 2 * M, h: 0.62, fontSize: 30, bold: true,
    color: NAVY, fontFace: "Calibri", margin: 0,
  });
  if (subtitle) {
    slide.addText(subtitle, {
      x: M, y: 0.92, w: W - 2 * M, h: 0.34, fontSize: 13, color: GRAY,
      fontFace: "Calibri", margin: 0,
    });
  }
}

// config chip: states the FULL timing configuration behind a plot/table
function configChip(slide, text, y) {
  slide.addShape("roundRect", {
    x: M, y: y, w: W - 2 * M, h: 0.6, fill: { color: ICE_LIGHT },
    line: { color: ICE, width: 1 }, rectRadius: 0.06,
  });
  slide.addText([{ text: "Cut configuration:  ", options: { bold: true, color: NAVY } },
                 { text: text, options: { color: "222222" } }], {
    x: M + 0.15, y: y, w: W - 2 * M - 0.3, h: 0.6, fontSize: 11.5,
    fontFace: "Calibri", margin: 0, valign: "middle",
  });
}

function bullets(slide, items, opts) {
  const runs = items.map((t, i) => {
    const o = { bullet: { code: "2022", indent: 12 }, paraSpaceAfter: 8, breakLine: true };
    if (typeof t === "string") return { text: t, options: o };
    return { text: t.text, options: Object.assign(o, t.options || {}) };
  });
  slide.addText(runs, Object.assign({
    fontSize: 14, color: "222222", fontFace: "Calibri", valign: "top", margin: 0,
  }, opts));
}

function footer(slide, n) {
  slide.addText(`MAIA calorimeter timing-cut study  ·  ${n}`, {
    x: W - 4.2, y: H - 0.38, w: 4.2 - M + 0.3, h: 0.3, fontSize: 9,
    color: GRAY, align: "right", fontFace: "Calibri", margin: 0,
  });
}

const TBL_BASE = {
  fontFace: "Calibri", fontSize: 12, color: "222222",
  border: { type: "solid", color: "D8DEE8", pt: 0.75 }, valign: "middle",
};
function headerRow(cells) {
  return cells.map((c) => ({ text: c, options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "center", fontSize: 12 } }));
}
function row(cells, highlight) {
  return cells.map((c, i) => ({
    text: c,
    options: {
      align: i === 0 ? "left" : "center",
      fill: { color: highlight ? "F9E5E2" : WHITE },
      bold: i === 0 && highlight,
      color: highlight && i === 0 ? ACCENT : "222222",
    },
  }));
}

// ---------------------------------------------------------------- 1 · title
{
  const s = pres.addSlide();
  s.background = { color: NAVY };
  s.addText("Calorimeter timing cuts in the MAIA reconstruction chain", {
    x: M, y: 2.1, w: W - 2 * M, h: 1.5, fontSize: 40, bold: true, color: WHITE,
    fontFace: "Calibri", margin: 0,
  });
  s.addText("How the ECAL and HCAL timing windows remove hits from clusters,\nand what cut values avoid shaping the cluster energy", {
    x: M, y: 3.6, w: W - 2 * M, h: 1.0, fontSize: 18, color: ICE, fontFace: "Calibri", margin: 0,
  });
  s.addText("MAIA_v0  ·  mucoll-sim-ubuntu24:v3.0 Gaudi chain  ·  no BIB overlay  ·  guns at θ = 90°\n" +
            "branch claude/maia-timing-cuts-study-wbtv1n  ·  full write-up in timing_studies/PLAN.md", {
    x: M, y: 5.6, w: W - 2 * M, h: 0.9, fontSize: 13, color: "8FA8D8", fontFace: "Calibri", margin: 0,
  });
}

// ------------------------------------------------------------- 2 · context
{
  const s = pres.addSlide();
  titleBar(s, "Context and goal", "Why timing cuts, and what “not shaping the energy” means here");
  bullets(s, [
    { text: "Beam-induced background (BIB) arrives out of time, so the reconstruction applies aggressive timing cuts to calorimeter hits.", options: {} },
    { text: "These cuts are enforced at more than one point in the chain, and their combined effect on signal clusters had not been quantified.", options: {} },
    { text: "Goal: measure how each timing cut removes hits/energy from ECAL and HCAL clusters, and suggest windows that do not strongly shape the cluster energy.", options: {} },
    { text: "This round is deliberately BIB-free: it isolates the signal side of the trade-off. The BIB side (occupancy re-admitted per window) is the natural follow-up.", options: {} },
    { text: "All times are relative:  Δt = t_arrival − r_cell/c  (time-of-flight from the IP at the speed of light already subtracted — the convention the code itself uses).", options: { bold: true } },
  ], { x: M, y: 1.55, w: 6.6, h: 5.0 });

  // simple chain diagram, right side
  const bx = 7.7, bw = 4.9, bh = 0.78;
  const steps = [
    ["Geant4 SimCalorimeterHits", "every energy deposit keeps its own time (detailed shower mode)"],
    ["Digitization  ·  RealisticCaloDigi", "per-contribution window Δt ∈ (−0.5, +10) ns  +  cell threshold"],
    ["Selection  ·  CaloHitSelector", "cell time window (−0.3, +0.3) ns  +  energy threshold → *Sel"],
    ["PandoraPFA clustering", "no timing cuts — clusters inherit whatever survives"],
  ];
  steps.forEach((st, i) => {
    const y = 1.7 + i * 1.32;
    s.addShape("roundRect", { x: bx, y: y, w: bw, h: bh + 0.28, fill: { color: i === 3 ? ICE_LIGHT : ICE }, line: { color: NAVY, width: 1 }, rectRadius: 0.07 });
    s.addText([{ text: st[0] + "\n", options: { bold: true, fontSize: 13, color: NAVY } },
               { text: st[1], options: { fontSize: 11, color: "333333" } }], {
      x: bx + 0.18, y: y + 0.06, w: bw - 0.36, h: bh + 0.16, fontFace: "Calibri", margin: 0, valign: "middle",
    });
    if (i < 3) s.addShape("downArrow", { x: bx + bw / 2 - 0.11, y: y + bh + 0.30, w: 0.22, h: 0.24, fill: { color: NAVY } });
  });
  footer(s, 2);
}

// ------------------------------------------------------------ 3 · software
{
  const s = pres.addSlide();
  titleBar(s, "Setup and samples", "Tutorial-tagged Gaudi chain (mcd-wiki “Running with Gaudi”), fully containerized and reproducible");
  bullets(s, [
    "Container ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v3.0; steering from mucoll-benchmarks + MAIAConfig; algorithms from k4Reco; geometry MAIA_v0.",
    "MAIAConfig pinned to d290d02 (2026-06-17): main requires a newer k4ActsTracking than the v3.0 container ships (the tutorial's known compatibility warning). Calo settings identical to main.",
    "Chain per sample: particle gun → ddsim (FTFP_BERT_HP, detailed shower mode: contributions keep individual times) → k4run digi → k4run reco (Pandora).",
    "Every quirk needed to make this reproducible (entrypoint, threshold-map location, EvtMax, gun steering) is scripted in timing_studies/production/ and documented in PLAN.md.",
  ], { x: M, y: 1.5, w: 7.2, h: 3.4 });
  const tbl = [
    headerRow(["samples", "energies", "events", "θ", "BIB"]),
    row(["photon gun  (ECAL phase)", "10 / 50 / 200 GeV", "100 per point", "90° (barrel)", "none"]),
    row(["π⁻ gun  (ECAL+HCAL phase)", "10 / 50 / 200 GeV", "100 per point", "90° (barrel)", "none"]),
  ];
  s.addTable(tbl, Object.assign({ x: M, y: 5.3, w: W - 2 * M, colW: [4.13, 2.6, 2.2, 1.8, 1.4] }, TBL_BASE));
  s.addText("Barrel-only by construction (θ = 90°); endcaps and θ dependence deferred to a follow-up θ scan.", {
    x: M, y: 6.55, w: W - 2 * M, h: 0.3, fontSize: 11, italic: true, color: GRAY, fontFace: "Calibri", margin: 0 });
  footer(s, 3);
}

// ------------------------------------------------- 4 · where the chain cuts
{
  const s = pres.addSlide();
  titleBar(s, "Where the chain cuts on time", "Verified against k4Reco / MAIAConfig source — two timing-cut levels act in sequence; Pandora adds none");
  const tbl = [
    headerRow(["level", "algorithm (instances)", "variable the cut acts on", "window / value"]),
    row(["1 · digitization", "RealisticCaloDigiSilicon / ScinPpd  (ECAL/HCAL, barrel+endcap)",
         "each Geant4 contribution:  Δt = t − r_cell/c", "(−0.5, +10) ns"], false),
    row(["1b · digi threshold", "same", "cell energy summed over IN-WINDOW contributions", "ECAL 5×10⁻⁵ GeV (≈0.32 MIP)  ·  HCAL 0.5 MIP"], false),
    row(["2 · selection", "CaloHitSelector  (→ *CollectionSel, what Pandora sees)",
         "stored digi hit time = EARLIEST accepted contribution Δt", "(−0.3, +0.3) ns, exclusive"], true),
    row(["2b · sel. threshold", "same", "calibrated cell energy", "ECAL: BIB mode map ≈ 1 MIP  ·  HCAL: 5×10⁻⁵ GeV"], false),
    row(["3 · clustering", "DDPandoraPFANewAlgorithm", "—", "no timing cuts"], false),
  ];
  s.addTable(tbl, Object.assign({ x: M, y: 1.6, w: W - 2 * M, colW: [1.85, 3.55, 3.63, 3.1], fontSize: 11.5 }, TBL_BASE));
  bullets(s, [
    { text: "The key interplay: a cell's ENERGY integrates everything out to +10 ns, but its TIME is the earliest accepted contribution. The ±0.3 ns selector therefore drops whole cells whose first deposit is late, while keeping every late deposit that shares a cell with a prompt one.", options: { bold: true } },
    "Tightening the digi window also kills cells through the threshold (1b acts on the windowed energy) — the two levels are coupled, which is why both are scanned together below.",
  ], { x: M, y: 4.9, w: W - 2 * M, h: 1.9, fontSize: 13 });
  footer(s, 4);
}

// --------------------------------------------------------- 5 · code findings
{
  const s = pres.addSlide();
  titleBar(s, "Things found while auditing the code", "Worth knowing independently of the tuning (details + line refs in PLAN.md §1)");
  const cards = [
    ["Double TOF correction — saved by a units accident",
     "The digi hit time is already propagation-corrected, yet CaloHitSelector subtracts r/TMath::C() again. With r in mm and C in m/s the second term is ~10⁻⁶ ns, so the cut works on the intended variable — but only by accident. If someone “fixes” the units it becomes a real double correction. → report upstream."],
    ["cellID truncation in RealisticCaloReco (k4Reco)",
     "int cellid = hit->getCellID() truncates 64 → 32 bits: every Rec/Coned/Sel calorimeter hit loses its x/y fields. Layer/stave/module survive, so calibration and the selector still work, but any consumer of full Sel-hit cellIDs gets degenerate IDs. → report upstream."],
    ["BIB-derived thresholds bite without BIB",
     "The ECAL selector threshold map (≈ 6.5 MeV ≈ 1 MIP per cell at θ = 90°) is derived from BIB occupancy but applied always. This study keeps it in numerator AND denominator so the timing effect is isolated — but it matters when interpreting absolute cluster energies."],
    ["Stochastic digitization (validation nuance)",
     "ECAL silicon: Poisson e–h smear (sub-%); HCAL scintillator: binomial SiPM pixel smear (~26% at 1 MIP), neither emulated. Sets the closure tolerances; timing behaviour itself is deterministic and reproduced exactly."],
  ];
  cards.forEach((c, i) => {
    const cx = M + (i % 2) * 6.2, cy = 1.55 + Math.floor(i / 2) * 2.7;
    s.addShape("roundRect", { x: cx, y: cy, w: 5.95, h: 2.5, fill: { color: i < 2 ? ICE_LIGHT : WHITE }, line: { color: ICE, width: 1.25 }, rectRadius: 0.08 });
    s.addText([{ text: c[0] + "\n", options: { bold: true, fontSize: 13.5, color: NAVY } },
               { text: c[1], options: { fontSize: 11.5, color: "222222" } }], {
      x: cx + 0.2, y: cy + 0.12, w: 5.55, h: 2.26, fontFace: "Calibri", margin: 0, valign: "top" });
  });
  footer(s, 5);
}

// -------------------------------------------------------------- 6 · method
{
  const s = pres.addSlide();
  titleBar(s, "Method: simhit-level emulation of the full chain", "Why simhit level: digitization lumps contributions together — only the sim record keeps each deposit's time");
  const steps = [
    ["1", "Flatten", "read SimCalorimeterHit contributions from the EDM4hep reco file (podio); per contribution: Δt = t − |cell position|/c and energy; per cell: position, layer, plus the ACTUAL Digi/Sel hits and Pandora clusters for closure."],
    ["2", "Emulate digitization", "for any window (min, max): cell energy = Σ in-window contributions; cell time = earliest accepted Δt; mean detector response (HCAL SiPM saturation) + thresholds exactly as in k4Reco."],
    ["3", "Emulate selection", "reco calibration → selector energy threshold (extracted ECAL map / flat HCAL) → selector time window on the cell time. Both windows are free parameters."],
    ["4", "Scan", "full 2D grid: digi max × selector window (symmetric and asymmetric), vectorized numpy → per-event retained energy at every grid point in seconds."],
  ];
  steps.forEach((st, i) => {
    const y = 1.6 + i * 1.15;
    s.addShape("ellipse", { x: M, y: y + 0.12, w: 0.52, h: 0.52, fill: { color: NAVY } });
    s.addText(st[0], { x: M, y: y + 0.12, w: 0.52, h: 0.52, align: "center", valign: "middle", color: WHITE, bold: true, fontSize: 16, fontFace: "Calibri", margin: 0 });
    s.addText([{ text: st[1] + " — ", options: { bold: true, color: NAVY, fontSize: 13.5 } },
               { text: st[2], options: { fontSize: 12, color: "222222" } }], {
      x: M + 0.75, y: y, w: 8.0, h: 1.1, fontFace: "Calibri", margin: 0, valign: "top" });
  });
  s.addShape("roundRect", { x: 9.7, y: 1.6, w: 3.0, h: 4.5, fill: { color: NAVY }, rectRadius: 0.08 });
  s.addText([
    { text: "Structural trick\n", options: { bold: true, fontSize: 14, color: WHITE } },
    { text: "With the digi window MIN fixed, the cell's earliest-time variable does not depend on the digi MAX — the max edge only changes cell energy and survival. The 2D scan is exact and cheap.\n\n", options: { fontSize: 11.5, color: ICE } },
    { text: "Emulation, not re-running\n", options: { bold: true, fontSize: 14, color: WHITE } },
    { text: "One sim pass per sample; every cut point afterwards is a numpy reduction, validated against genuine chain re-runs (next slide).", options: { fontSize: 11.5, color: ICE } },
  ], { x: 9.9, y: 1.75, w: 2.6, h: 4.2, fontFace: "Calibri", margin: 0, valign: "top" });
  footer(s, 6);
}

// ---------------------------------------------------------- 7 · validation
{
  const s = pres.addSlide();
  titleBar(s, "How the emulation was validated", "Closure gates run per sample before any scan is trusted; then genuine re-runs of the chain at shifted cut settings");
  const stats = [
    ["< 0.1%", "ECAL closure, per event", "emulated vs actual Sel-collection energy sums, at default cuts AND every variant point"],
    ["0.3 – 1.5%", "HCAL closure, per event", "entirely the un-emulated binomial SiPM smear; → 0.3% at high occupancy (200 GeV π)"],
    ["× 1.0237", "leading cluster / Sel sum", "exactly Pandora's EM calibration: cluster energy tracks hit selection — hit-level fractions ARE cluster-level fractions"],
  ];
  stats.forEach((st, i) => {
    const cx = M + i * 4.15;
    s.addShape("roundRect", { x: cx, y: 1.5, w: 3.9, h: 2.15, fill: { color: ICE_LIGHT }, line: { color: ICE, width: 1 }, rectRadius: 0.08 });
    s.addText(st[0], { x: cx, y: 1.6, w: 3.9, h: 0.7, align: "center", fontSize: 30, bold: true, color: NAVY, fontFace: "Calibri", margin: 0 });
    s.addText([{ text: st[1] + "\n", options: { bold: true, fontSize: 12, color: "222222" } },
               { text: st[2], options: { fontSize: 10, color: GRAY } }], {
      x: cx + 0.18, y: 2.32, w: 3.54, h: 1.28, align: "center", fontFace: "Calibri", margin: 0 });
  });
  // left: the gate ladder, run per sample before any scan is trusted
  s.addText([{ text: "Closure-gate ladder (run per sample, must pass before scanning)", options: { bold: true, fontSize: 13.5, color: NAVY, breakLine: true, paraSpaceAfter: 6 } }].concat([
    "1 · flatten self-consistency: Σ contribution energies = SimCalorimeterHit energy (exact).",
    "2 · prompt Δt peak at ≈ 0 with sub-ns width → no unit or r/c-offset errors.",
    "3 · per-cell join, emulated vs ACTUAL Digi hits (full 64-bit cellIDs): identical cell sets to <0.05%, energies to the stochastic tolerance, times to 10⁻⁶ ns.",
    "4 · emulated vs ACTUAL Sel hits per (event, layer): validates threshold-map lookup, MIP thresholds, SiPM saturation and the time window in one comparison.",
    "5 · event-level Sel energy sums (numbers above).",
  ].map((t) => ({ text: t, options: { bullet: { code: "2022", indent: 12 }, paraSpaceAfter: 5, breakLine: true, fontSize: 11.5, color: "222222" } }))), {
    x: M, y: 3.95, w: 6.05, h: 3.1, fontFace: "Calibri", margin: 0, valign: "top" });
  // right: the real-chain variant matrix
  s.addText([{ text: "6 · Real-chain variant re-runs (the decisive check)", options: { bold: true, fontSize: 13.5, color: NAVY, breakLine: true, paraSpaceAfter: 6 } }].concat([
    "k4run digi+reco genuinely re-run on the same sim files with property overrides on all four Digi / all four Selector instances coherently.",
    "Matrix: {selector ±0.15 ns, selector ±1.0 ns, digi max 2 ns} × {γ, π⁻} × {10, 50, 200 GeV} — 18 variant reconstructions.",
    "Emulation matches every variant's Sel sums per event to <0.1% (ECAL) / ~1% median (HCAL), and the Pandora cluster follows the Sel sum at each point.",
  ].map((t) => ({ text: t, options: { bullet: { code: "2022", indent: 12 }, paraSpaceAfter: 5, breakLine: true, fontSize: 11.5, color: "222222" } })).concat([
    { text: "⇒ one sim pass + emulation is equivalent to re-running digitization and selection at every scan point.", options: { bold: true, fontSize: 12, color: NAVY, breakLine: false } },
  ])), {
    x: 7.05, y: 3.95, w: 5.65, h: 3.1, fontFace: "Calibri", margin: 0, valign: "top" });
  footer(s, 7);
}

// ----------------------------------------------- 8 · photon dt distributions
{
  const s = pres.addSlide();
  titleBar(s, "Photons: where the ECAL energy is in time");
  configChip(s, "spectra built at the DEFAULT digi window (−0.5, +10) ns; no selector cut applied — guide lines mark the current selector ±0.3 ns and digi max 10 ns edges. Cell time = earliest accepted contribution.", 1.35);
  s.addImage({ path: `${PLOTS}/ecal_photons/dt_cellearliest_EcalBarrel.png`, x: M, y: 2.0, w: 6.55, h: 4.68 });
  bullets(s, [
    "Energy-weighted Δt of ECAL barrel cells: a prompt peak at Δt ≈ +0.02 ns carrying ≈99.5% of the energy within ±0.5 ns.",
    "The spectrum falls three orders of magnitude within ~1 ns, then a flat tail at the 10⁻⁵/ns level extends past 10 ns.",
    "Identical shape at 10 / 50 / 200 GeV → whatever a window cuts, it cuts nearly the same fraction at every energy (little scope for shaping — quantified next).",
    "Same conclusion at contribution level (no cuts at all): plot in plots/ecal_photons/dt_contrib_EcalBarrel.png.",
  ], { x: 7.5, y: 2.2, w: 5.2, h: 4.4, fontSize: 13.5 });
  footer(s, 8);
}

// --------------------------------------------------- 9 · ECAL selector scan
{
  const s = pres.addSlide();
  titleBar(s, "ECAL / photons — selector window scan");
  configChip(s, "digi window FIXED at default (−0.5, +10) ns; selector window varied, symmetric ±w. Denominator everywhere: same energy thresholds, NO timing cuts at either level (digi max 100 ns, selector open).", 1.35);
  s.addImage({ path: `${PLOTS}/ecal_photons/retained_vs_selwidth.png`, x: M, y: 2.0, w: 4.85, h: 4.85 });
  const tbl = [
    headerRow(["selector window", "10 GeV", "50 GeV", "200 GeV"]),
    row(["±0.3 ns  (current)", "0.9942", "0.9940", "0.9941"], true),
    row(["±0.5 ns", "0.9983", "0.9979", "0.9982"], false),
    row(["±1.0 ns", "1.0000", "0.9990", "0.9991"], false),
    row(["−0.3 < t < +0.5 ns", "0.9983", "0.9979", "0.9982"], false),
    row(["−0.3 < t < +1.0 ns", "1.0000", "0.9990", "0.9991"], false),
  ];
  s.addTable(tbl, Object.assign({ x: 5.75, y: 2.0, w: 7.0, colW: [2.8, 1.4, 1.4, 1.4], fontSize: 11.5 }, TBL_BASE));
  bullets(s, [
    { text: "Current ±0.3 ns: a flat 0.6% loss, identical at all three energies — an offset, not a shape.", options: { bold: true } },
    "Asymmetric rows equal symmetric ones → ALL loss is on the late side; the early edge (−0.3 ns) removes nothing.",
    "Event-to-event spread of the retained fraction (a resolution term): 0.46 / 0.21 / 0.13% at ±0.3 ns, halving at ±0.5 ns.",
  ], { x: 5.75, y: 4.55, w: 7.0, h: 2.2, fontSize: 13 });
  footer(s, 9);
}

// ------------------------------------------------------- 10 · ECAL digi scan
{
  const s = pres.addSlide();
  titleBar(s, "ECAL / photons — digi window scan");
  configChip(s, "digi max varied with digi min FIXED at −0.5 ns. Two curve families: selector OPEN (solid — isolates the digi cut) and selector at the DEFAULT ±0.3 ns (dashed — the two cuts combined). Denominator as before: no timing cuts.", 1.35);
  s.addImage({ path: `${PLOTS}/ecal_photons/retained_vs_digimax.png`, x: M, y: 2.05, w: 6.55, h: 4.68 });
  bullets(s, [
    { text: "The digi window is a non-issue for the ECAL: with the selector open, even max = 1 ns keeps > 99.8%; at 2 ns the loss is < 0.1% at all energies.", options: { bold: true } },
    "With the selector at ±0.3 ns (dashed), the curves saturate at the selector's own 99.4% — beyond ~2 ns the digi edge is irrelevant because the selector dominates.",
    "Implication: the current 10 ns is generously safe — or can be tightened to ~2–3 ns for free if BIB energy integrated in the (0.3, 10) ns tail of prompt cells proves harmful.",
    "Digi min scan (−0.25 / −0.5 / −1.0 ns): no effect on photons.",
  ], { x: 7.5, y: 2.2, w: 5.2, h: 4.4, fontSize: 13.5 });
  footer(s, 10);
}

// --------------------------------------------------- 11 · pion distributions
{
  const s = pres.addSlide();
  titleBar(s, "Pions: hadronic showers are late — by construction");
  configChip(s, "spectra built at the DEFAULT digi window (−0.5, +10) ns; no selector cut applied — guide lines mark the current selector ±0.3 ns and digi max 10 ns edges.", 1.35);
  s.addImage({ path: `${PLOTS}/hcal_pions/dt_cellearliest_HcalBarrel.png`, x: M, y: 2.05, w: 6.55, h: 4.68 });
  bullets(s, [
    { text: "Only ~59% of the HCAL and ~78% of the ECAL pion energy arrives within ±0.5 ns.", options: { bold: true } },
    "The neutron-dominated late component spreads over tens of ns in cell earliest-time (and hundreds of ns at contribution level) — it is intrinsic to hadronic showers, not an artifact.",
    "Many late deposits land in cells whose FIRST deposit is prompt: those survive the selector via the +10 ns digi integration. Cells that light up late are dropped whole — the dominant hadronic loss mechanism.",
    "So for hadrons the window choice is a genuine trade-off; no reasonable window recovers everything.",
  ], { x: 7.5, y: 2.2, w: 5.2, h: 4.4, fontSize: 13.5 });
  footer(s, 11);
}

// ------------------------------------------------- 12 · pion combined result
{
  const s = pres.addSlide();
  titleBar(s, "Pions: what the current cuts do to cluster energy");
  configChip(s, "each row states BOTH cut levels explicitly. Retained fraction = median over events, ECAL+HCAL barrel summed; 68% event spread in parentheses. Denominator: same thresholds, no timing cuts.", 1.35);
  const tbl = [
    headerRow(["digi window  +  selector window", "10 GeV", "50 GeV", "200 GeV"]),
    row(["(−0.5, +10) ns  +  ±0.3 ns   ← current defaults", "0.709 (0.129)", "0.789 (0.066)", "0.845 (0.046)"], true),
    row(["(−0.5, +10) ns  +  ±1.0 ns", "0.883 (0.071)", "0.908 (0.034)", "0.929 (0.022)"], false),
    row(["(−0.5, +10) ns  +  ±2.0 ns", "0.923 (0.051)", "0.938 (0.024)", "0.952 (0.017)"], false),
    row(["(−0.5, +10) ns  +  ±5.0 ns", "0.965 (0.025)", "0.969 (0.013)", "0.975 (0.009)"], false),
    row(["(−0.5, +10) ns  +  selector OPEN  (digi cut alone)", "0.983 (0.013)", "0.987 (0.007)", "0.987 (0.005)"], false),
    row(["(−0.5, +25) ns  +  selector OPEN", "0.994 (0.007)", "0.996 (0.003)", "0.996 (0.002)"], false),
  ];
  s.addTable(tbl, Object.assign({ x: M, y: 2.0, w: W - 2 * M, colW: [5.53, 2.2, 2.2, 2.2], fontSize: 11.5 }, TBL_BASE));
  bullets(s, [
    { text: "The current defaults remove 15–29% of the hadronic calorimeter energy, with ~14% response non-linearity between 10 and 200 GeV and a 5–13% event-to-event spread — genuine cluster-energy shaping.", options: { bold: true, color: ACCENT } },
    "Cross-checked in the REAL chain (not emulation): the median 10 GeV π leading Pandora cluster moves 5.35 → 8.90 GeV between selector ±0.15 ns and ±1.0 ns (digi window at default in both).",
  ], { x: M, y: 4.95, w: W - 2 * M, h: 1.8, fontSize: 13 });
  footer(s, 12);
}

// -------------------------------------------- 13 · ECAL tuning (with pions)
{
  const s = pres.addSlide();
  titleBar(s, "Tuning the ECAL window alone (\u03b3 + \u03c0)");
  configChip(s, "ECAL barrel only; digi window FIXED at default (−0.5, +10) ns; ECAL selector window varied. Selector windows factorize per cell, so ECAL and HCAL tune independently. Denominator: no timing cuts.", 1.35);
  s.addImage({ path: `${PLOTS}/hcal_pions/retained_vs_selwidth_ecal.png`, x: M, y: 2.15, w: 4.4, h: 4.4 });
  const tbl = [
    headerRow(["ECAL window", "γ 10", "γ 200", "π 10", "π 50", "π 200", "π 200/10"]),
    row(["±0.3 ns (cur.)", "0.994", "0.994", "0.794", "0.865", "0.916", "1.153"], true),
    row(["±0.5 ns", "0.998", "0.998", "0.870", "0.911", "0.944", "1.084"], false),
    row(["±1.0 ns", "1.000", "0.999", "0.919", "0.944", "0.964", "1.049"], false),
    row(["±2.0 ns", "1.000", "1.000", "0.959", "0.965", "0.973", "1.016"], false),
    row(["±5.0 ns", "1.000", "1.000", "0.979", "0.984", "0.987", "1.008"], false),
  ];
  s.addTable(tbl, Object.assign({ x: 5.3, y: 2.15, w: 7.45, colW: [1.75, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95], fontSize: 11 }, TBL_BASE));
  bullets(s, [
    "π 200/10 = retained(200 GeV) / retained(10 GeV): the response-linearity (shaping) metric.",
    { text: "Photons alone would let the ECAL sit at +0.5 ns — but the ECAL half of hadronic showers is also late.", options: {} },
    { text: "At ±0.3 ns the pion ECAL response carries 15% non-linearity; (−0.3, +2.0) ns brings it to 1.6% with photons still exactly lossless.", options: { bold: true } },
    "If BIB forces tighter: +1.0 ns is the compromise (5% non-linearity).",
  ], { x: 5.3, y: 4.7, w: 7.45, h: 1.9, fontSize: 12.5 });
  footer(s, 13);
}

// -------------------------------------------------------- 14 · HCAL tuning
{
  const s = pres.addSlide();
  titleBar(s, "Tuning the HCAL window on its own (pions)");
  configChip(s, "HCAL barrel only; digi window FIXED at default (−0.5, +10) ns; HCAL selector window varied. Denominator: no timing cuts. Photons contribute only ~0.05% leakage here — pions are the relevant probe.", 1.35);
  s.addImage({ path: `${PLOTS}/hcal_pions/retained_vs_selwidth_hcal.png`, x: M, y: 2.15, w: 4.4, h: 4.4 });
  const tbl = [
    headerRow(["HCAL window", "π 10", "π 50", "π 200", "π 200/10"]),
    row(["±0.3 ns (cur.)", "0.515 (0.37)", "0.613 (0.16)", "0.762 (0.10)", "1.479"], true),
    row(["±1.0 ns", "0.785 (0.31)", "0.839 (0.11)", "0.886 (0.05)", "1.129"], false),
    row(["±2.0 ns", "0.846 (0.19)", "0.886 (0.07)", "0.919 (0.04)", "1.087"], false),
    row(["±5.0 ns", "0.934 (0.10)", "0.947 (0.03)", "0.960 (0.02)", "1.028"], false),
    row(["open", "0.983 (0.04)", "0.978 (0.02)", "0.982 (0.01)", "1.000"], false),
  ];
  s.addTable(tbl, Object.assign({ x: 5.3, y: 2.15, w: 7.45, colW: [1.65, 1.55, 1.55, 1.55, 1.15], fontSize: 11 }, TBL_BASE));
  bullets(s, [
    { text: "The current ±0.3 ns keeps only 52–76% of the HCAL energy, with 48% non-linearity and up to 37% event spread — this cut, not the digi window, drives the hadronic shaping.", options: { bold: true, color: ACCENT } },
    "Each doubling of the late edge roughly halves loss and spread; the response flattens around +5 ns (2.8% non-linearity).",
    "Loosening the digi window 10 → 25 ns adds 1–2% more, relevant mainly once the selector is at +5 ns.",
  ], { x: 5.3, y: 4.7, w: 7.45, h: 1.9, fontSize: 12.5 });
  footer(s, 14);
}

// --------------------------------- 15 · why the 10 GeV spread is non-monotonic
{
  const s = pres.addSlide();
  titleBar(s, "The 10 GeV HCAL spread bump, explained",
    "The 68% band on the previous slide is a percentile of a distribution that is far from Gaussian at low energy — here are the distributions themselves");
  configChip(s, "HCAL barrel, π⁻; digi window FIXED at default (−0.5, +10) ns; selector windows in the legends. f = per-event E(window) / E(no timing cuts) — the same quantity whose median and 68% band appear on the previous slide.", 1.35);
  s.addImage({ path: `${PLOTS}/hcal_pions/fraction_dist_HcalBarrel.png`, x: M + 0.35, y: 2.08, w: 11.45, h: 3.44 });
  bullets(s, [
    { text: "10 GeV: a distinct population retains f = 0 exactly — 17 events at the current ±0.3 ns, still 9 at ±2 ns: no HCAL cell has an early-enough FIRST deposit (the HCAL piece is small — median ~2 GeV, ~15 cells — late and neutron-dominated). The 16th percentile sits at 0 up to w = 0.3 ns.", options: { bold: true } },
    "Opening the window pulls prompt-rich events to f ≈ 1 while late-dominated events stay near 0: the distribution stretches across [0, 1], so the 68% spread PEAKS (~0.40 at w ≈ 0.4–0.75 ns) where the median crosses mid-range — a p(1−p)-type maximum — then collapses as everything converges to f ≈ 1. The rise-then-fall of the blue curve is expected, not a glitch.",
    "50 / 200 GeV: showers average over many more cells and sub-showers — a single peak that narrows and moves right, no f = 0 spike, monotonic spread. At 10 GeV and ±0.3 ns the HCAL flips between seeing the shower and seeing nothing, event by event — not a smeared measurement.",
  ], { x: M, y: 5.62, w: W - 2 * M, h: 1.6, fontSize: 10.5 });
  footer(s, 15);
}

// ------------------------------------------------------- 16 · conclusions
{
  const s = pres.addSlide();
  s.background = { color: NAVY };
  s.addText("Suggested timing cuts (no-BIB baseline)", {
    x: M, y: 0.45, w: W - 2 * M, h: 0.7, fontSize: 30, bold: true, color: WHITE, fontFace: "Calibri", margin: 0 });
  const panels = [
    ["ECAL selector", "(−0.3, +2.0) ns", "photons stay exactly lossless; hadronic-ECAL non-linearity 15% → 1.6%. Compromise under BIB pressure: +1.0 ns. EM-only use case: +0.5 ns suffices."],
    ["HCAL selector", "(−0.3, +5.0) ns", "response flattens (non-linearity 2.8%, 93–96% retained). Floor: +2.0 ns (8.7%). The current ±0.3 ns should not be used for hadronic calorimetry unless BIB demands it."],
    ["digi windows", "ECAL 2–3 ns  ·  HCAL ~25 ns", "secondary knobs: ECAL can tighten for free (BIB in the integration tail); HCAL loosening recovers 1–2%, relevant once the selector is at +5 ns."],
  ];
  panels.forEach((p, i) => {
    const cx = M + i * 4.15;
    s.addShape("roundRect", { x: cx, y: 1.45, w: 3.9, h: 3.0, fill: { color: "2A3775" }, rectRadius: 0.08 });
    s.addText([{ text: p[0] + "\n", options: { fontSize: 14, bold: true, color: ICE } },
               { text: p[1] + "\n", options: { fontSize: 24, bold: true, color: WHITE } },
               { text: p[2], options: { fontSize: 11.5, color: "B9C8E8" } }], {
      x: cx + 0.22, y: 1.62, w: 3.46, h: 2.7, fontFace: "Calibri", margin: 0, valign: "top" });
  });
  const closing = [
    "Early edges stay at −0.3 ns everywhere — they remove nothing from signal at any energy, for either particle type.",
    "Any selector retune shifts the hadronic energy scale by 10–25%: Pandora's hadronic calibration and software-compensation weights must be re-derived alongside.",
    "Caveats: no BIB (the other side of the trade-off), θ = 90° barrel only, 100 events/point.",
    "Next steps: repeat with BIB overlay to plot signal-retained vs BIB-admitted per window on the same axes; θ scan for the endcaps; upstream reports for the two code findings.",
  ];
  s.addText(closing.map((t, i) => ({ text: t, options: { bullet: { code: "2022", indent: 12 }, paraSpaceAfter: 7, breakLine: true, fontSize: 12.5, color: ICE } })), {
    x: M, y: 4.85, w: W - 2 * M, h: 2.2, fontFace: "Calibri", margin: 0, valign: "top" });
  s.addText("code, plots and full write-up: timing_studies/ on branch claude/maia-timing-cuts-study-wbtv1n (PLAN.md)", {
    x: M, y: H - 0.42, w: W - 2 * M, h: 0.3, fontSize: 10, color: "8FA8D8", fontFace: "Calibri", margin: 0 });
}

pres.writeFile({ fileName: "MAIA_timing_cuts_study.pptx" }).then(() => console.log("written"));
