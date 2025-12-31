import numpy as np
import matplotlib.pyplot as plt

# Load CSV (skips header automatically)
data = np.genfromtxt("timeseries.csv", delimiter=",", names=True)

t = data["t_s"]        # time [s]
Tinf = data["Tinf_K"]  # redshifted internal temperature [K]

# Guard against non-positive entries
mask = (t > 0.0) & (Tinf > 0.0)
t = t[mask]
Tinf = Tinf[mask]

plt.figure(figsize=(6, 4))
plt.loglog(t, Tinf, marker="o", linestyle="-", markersize=4)

plt.xlabel(r"$\log_{10}(t\,[\mathrm{s}])$")
plt.ylabel(r"$\log_{10}(T_\infty\,[\mathrm{K}])$")
plt.title("Cooling Curve: $T_\\infty$ vs Time")

plt.grid(True, which="both", ls=":")
plt.tight_layout()
plt.show()