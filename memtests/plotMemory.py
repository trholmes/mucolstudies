import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

pid = sys.argv[1]
filename = f"memlog_{pid}.txt"
waittime = 5

file_vals = np.loadtxt(filename, dtype=str)
mem_vals = []
for v in file_vals:
    if v.endswith("g"):
        mem_vals.append(float(np.char.strip(v,"g")))
    elif v.endswith("m"):
        mem_vals.append(float(np.char.strip(v,"m"))/1000.)
    else:
        mem_vals.append(float(v)/1000000.)
t_vals = np.linspace(0, waittime*len(mem_vals), len(mem_vals))

plt.plot(t_vals, mem_vals)
plt.xlabel("Time [s]")
plt.ylabel("Memory Usage [GB]")
plt.savefig(f"memplot_{pid}.png")




