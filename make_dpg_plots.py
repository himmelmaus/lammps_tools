import numpy as np
from matplotlib import pyplot as plt
import sys

plt.rcParams["figure.figsize"] = (30,10)
plt.rcParams["font.size"] = 28
# diffusion plotting
cutoff=70

pol_bins, pol_ds = np.loadtxt("binned_pol_rabe.txt")
non_pol_bins, non_pol_ds = np.loadtxt("binned_msd_non_pol.txt")

water_bins, water_ds = np.loadtxt("binned_msd_water_test.txt")

bulk_d = np.mean(water_ds)
bulk_std = np.std(water_ds)

print(bulk_d, bulk_std)

non_pol_ds = np.roll(non_pol_ds, -np.argmin(non_pol_ds))
# non_pol_bins -= int(non_pol_bins[0])
# pol_bins -= int(pol_bins[0])

# pol_bins -= np.abs(pol_bins[1] - non_pol_bins[0])

pol_ds = pol_ds[(pol_bins>0) & (pol_bins<=cutoff)]
pol_bins = pol_bins[(pol_bins>0) & (pol_bins<=cutoff)]

non_pol_ds = non_pol_ds[(non_pol_bins>0) & (non_pol_bins<=cutoff)]
non_pol_bins = non_pol_bins[(non_pol_bins>0) & (non_pol_bins<=cutoff)]

# plt.fill_between([0, cutoff], [2.5e-9, 2.5e-9], [2.5e-9, ], label="Bulk Value")

plt.plot([0,cutoff], [2.2e-9,2.2e-9], "--", label="Bulk Value [4]")
plt.plot(non_pol_bins, non_pol_ds,"x-", label="Non-polarisable", color="black")
plt.plot(pol_bins, pol_ds, "x-", label="polarisable", color="red")
plt.xlabel("Height from gold surface (Å)")
plt.ylabel("Self-Diffusion Coefficient (${m^2}s^{-1}$)")

plt.legend()
plt.xlim(0, cutoff)
plt.tight_layout()
plt.savefig("rabe_diffusion.png")
plt.cla()
plt.clf()

# density plotting
cutoff=20

non_pol_bins, non_pol_density = np.loadtxt("non_pol_density_profile.txt").T
pol_bins, pol_density = np.loadtxt("pol_density_profile.txt").T

# pol_bins -= 8.4

filter = ~np.isclose(pol_density, 0)
pol_bins = pol_bins[filter]
pol_density = pol_density[filter]
pol_bins -= 8.5
print(min(pol_bins))

pol_density = pol_density[(pol_bins>0) & (pol_bins<=cutoff)]
pol_bins = pol_bins[(pol_bins>0) & (pol_bins<=cutoff)]

non_pol_density = non_pol_density[(non_pol_bins>0) & (non_pol_bins<=cutoff)]
non_pol_bins = non_pol_bins[(non_pol_bins>0) & (non_pol_bins<=cutoff)]



# plt.plot([0,cutoff], [20.98, 0.98], "--", label="Bulk Value")
plt.plot(non_pol_bins, non_pol_density,"-", label="Non-polarisable", color="black")
plt.plot(pol_bins, pol_density, "-", label="polarisable", color="red")
plt.plot([0,cutoff], [0.98,0.98], "--", label="Bulk Value [4]")
plt.xlabel("Height from gold surface (Å)")
plt.ylabel("Mass Density ($kgm^{-3}$)")
plt.xlim(0, cutoff)
plt.legend()
plt.tight_layout()

plt.savefig("rabe_density.png")
plt.cla()
plt.clf()

