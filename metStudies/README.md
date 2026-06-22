# metStudies

Study of the MET distribution in events with full BIB included.

The event loop lives in `plotMET.cpp` rather than python, since the BIB files are
huge and the pyLCIO loops in the rest of this repo are slow. Two things make it fast:

1. **Compiled event loop** — no python interpreter overhead per PFO.
2. **`setReadCollectionNames`** — the LCIO reader is told to only unpack the PFO
   collection, so it skips decompressing the calo/tracker hit collections that
   dominate the size of BIB files. This is the big win; if you want it in your
   python scripts too, it's one line on a pyLCIO reader:
   `reader.setReadCollectionNames(["PandoraPFOs"])`

## Building

This is meant to run in the MAIA software container (v2.11), as described on the
[MAIA howto](https://mcd-wiki.web.cern.ch/software/howto/maia/):

```
# Docker (local)
docker pull ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.11

# Apptainer (OSG/cluster)
apptainer run -B /ospool/uc-shared/project/futurecolliders/data/ \
  /cvmfs/unpacked.cern.ch/ghcr.io/muoncollidersoft/mucoll-sim-ubuntu24:v2.11-amd64
```

The container sets up the environment automatically; then just run:

```
make
```

This needs `root-config` on your PATH and `$LCIO` set, both of which the container
environment provides. If `$LCIO` isn't set, try
`export LCIO=$(dirname $(dirname $(which anajob)))`.

## Running

```
./plotMET [-n maxEvents] [-o out.root] [-c collectionName] file1.slcio [file2.slcio ...]
```

For example, over the MAIA v0 samples on OSG:

```
./plotMET -o met_BIB.root /ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MAIA_v0/v7/recoBIB/muonGun_pT_0_50/*.slcio
./plotMET -o met_noBIB.root /ospool/uc-shared/project/futurecolliders/data/fmeloni/DataMuC_MAIA_v0/v7/reco/muonGun_pT_0_50/*.slcio
```

(Adjust the path to wherever the data is mounted in your container, and check the
subdirectory layout with `ls` — the exact sample naming varies between productions.)

MET is computed as the magnitude of the negative vector sum of px,py over all
PFOs in `PandoraPFOs` (override with `-c` if your samples use a different
collection name — check with `anajob <file>`). The full 3D missing momentum
p_miss (negative vector sum including pz) and its z-component are also filled,
along with the total visible energy sum E — since sqrt(s) is known at a muon
collider, missing energy is just E_miss = sqrt(s) - sum E. All histograms (MET,
p_miss, low-range zooms of each, MET phi/x/y, p_miss_z, sum ET, sum E, nPFOs)
are written to the output ROOT file, and quick-look PNGs go to `plots/`.

To overlay results from different samples (normalized to unit area, using the
repo's `plotHelper` styling):

```
python compareMET.py "no BIB":met_noBIB.root "BIB":met_BIB.root
```
