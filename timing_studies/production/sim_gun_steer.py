"""ddsim steering for gun samples: identical to the mucoll-benchmarks
baseline steering, but with the input file cleared so --enableGun works.
Needs BENCHMARKS_DIR in the environment (set by run_chain.sh).
"""

import os

_base = os.path.join(os.environ["BENCHMARKS_DIR"], "simulation", "steer_baseline.py")
exec(compile(open(_base).read(), _base, "exec"))

SIM.inputFiles = []  # noqa: F821  (SIM comes from the baseline steering)
