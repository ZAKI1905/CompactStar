import numpy as np
import matplotlib.pyplot as plt

# Load CSV (skips header automatically)
data = np.genfromtxt("timeseries.csv", delimiter=",", names=True)

t = data["t_s"]             # time [s]
Tinf = data["Tinf_K"]       # redshifted internal temperature [K]
Ts = data["Tsurf_K"]        # surface temperature [K]
Lnu = data["L_nu_inf_erg_s"] # Neutrino luminosity
Lgamma = data["L_gamma_inf_erg_s"] # Neutrino luminosity

# Guard against non-positive entries
mask = (t > 0.0) & (Tinf > 0.0) & (Ts > 0.0)
t = t[mask]
Tinf = Tinf[mask]
Ts = Ts[mask]

plt.figure(figsize=(6, 4))

# plt.loglog(t, Tinf, marker="o", linestyle="-", markersize=4,
#            label=r"$T_\infty$ (interior)")
plt.loglog(t, Ts, marker="s", linestyle="--", markersize=0,
           label=r"$T_s$ (surface)")

# plt.loglog(t, Lnu, marker="s", linestyle="--", markersize=2,
#            label=r"$L_{\nu}$")
# plt.loglog(t, Lgamma, marker="s", linestyle="--", markersize=2,
#            label=r"$L_{\gamma}$")

# plt.loglog(t, Lnu+Lgamma, marker="s", linestyle="--", markersize=2,
#            label=r"$L_{cool}$")

plt.xlabel(r"$t\;[\mathrm{s}]$")
plt.ylabel(r"$T\;[\mathrm{K}]$")
plt.title("Neutron Star Cooling Curve")

plt.legend()
plt.grid(True, which="both", ls=":")
plt.tight_layout()
plt.show()