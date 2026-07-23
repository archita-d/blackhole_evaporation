import numpy as np
import matplotlib.pyplot as plt

G = 6.67e-11
c = 3e8
h = 6.626e-34
k = 1.381e-23
sigma = 5.670e-8
MeV_to_joules = 1.602176634e-13

def radius(M):
    return 2 * G * M / c**2

def temperature(M):
    return (h / (2 * 3.14159)) * c**3 / (8 * 3.14159 * G * M * k)

def power(M):
    r_s = radius(M)
    A = 4 * 3.14159 * r_s**2
    T = temperature(M)
    return sigma * A * T**4

def b_function(E_MeV, T):
    E = E_MeV * MeV_to_joules
    x = E / (k * T)
    x = np.minimum(x, 700)
    denominator = np.expm1(x)
    denominator = np.where(denominator == 0, 1e-300, denominator)
    return (2 * E**2) / (h**3 * c**2) / denominator

def step3(M, E_low, E_high, time_sample):
    T = temperature(M)
    r_s = radius(M)
    surface_area_m2 = 4 * 3.14159 * r_s**2   # NO cm conversion - b_function is per m^2
    E_values_MeV = np.linspace(E_low, E_high, 200)
    p_values = b_function(E_values_MeV, T)
    E_values_joules = E_values_MeV * MeV_to_joules
    integrated = np.trapezoid(p_values, E_values_joules)
    return integrated * surface_area_m2 * time_sample * 3.14159

def step4(M, E_low, E_high, time_sample, D, A_detector):
    emitted = step3(M, E_low, E_high, time_sample)
    detector_fraction = A_detector / (4 * 3.14159 * D**2)
    return emitted * detector_fraction

def gamma_array(M, energy_bin, time_sample, D, A_detector):
    # this IS "step_counts" - the count over each of the 32 energy bins,
    # over the time sample, at the earth
    counts = []
    for i in range(len(energy_bin) - 1):
        counts.append(step4(M, energy_bin[i], energy_bin[i + 1], time_sample, D, A_detector))
    return np.array(counts)


# from powerpoint slides basically 

M_test = 1.73e11        # kg - a snapshot mass to test on
D_test = 9.4607e15      # 1 light-year in meters
A_detector = 50e-4     # m^2 (this is just 50 cm^2 converted to m^2)
time_sample = 1e-6    # 1 microsecond

energy_bin = np.logspace(np.log10(0.001), np.log10(500.0), 33)
E_low = energy_bin[:-1]
E_high = energy_bin[1:]
bin_midpoint = (E_low + E_high) / 2

step_counts = gamma_array(M_test, energy_bin, time_sample, D_test, A_detector)

# 1) Make count spectrum: divide each counts by the energy bin width
bin_width = E_high - E_low
count_spectrum = step_counts / bin_width

# 2) Make energy spectrum: multiply each of the count spectrum by the mid-point of the bin
energy_spectrum = count_spectrum * bin_midpoint

# 3) Print out the count array, the count spectrum, and the energy spectrum
print("=== STEP COUNTS (raw, 32 bins) ===")
print(step_counts)

print("\n=== COUNT SPECTRUM (counts / bin width) ===")
print(count_spectrum)

print("\n=== ENERGY SPECTRUM (count spectrum * midpoint) ===")
print(energy_spectrum)

# im essentially checking here  where does the energy spectrum peak and how does that compare to kT?
T_test = temperature(M_test)
kT_MeV = (k * T_test) / MeV_to_joules
peak_idx = np.argmax(energy_spectrum)
peak_energy = bin_midpoint[peak_idx]
print(f"\nTemperature: {T_test:.3e} K")
print(f"kT in MeV: {kT_MeV:.3e}")
print(f"Energy spectrum peaks at bin centered on: {peak_energy:.3e} MeV")
print(f"Ratio of peak energy to kT: {peak_energy / kT_MeV:.3f}  (should be around 3.something)")

# plot it
plt.figure(figsize=(9, 6))
plt.plot(bin_midpoint, energy_spectrum, marker='o')
plt.xlabel("Energy (MeV)")
plt.ylabel("Energy Spectrum")
plt.title(f"Energy Spectrum Check (M={M_test:.3e} kg, T={T_test:.3e} K)")
plt.xscale('log')
plt.yscale('log')
plt.grid(True, which='both', alpha=0.4)
plt.tight_layout()
plt.savefig("step_counts_debug.png", dpi=200)
plt.show()

# redo steps 1 & 2 with CONSISTENT units...
bin_width_joules = bin_width * MeV_to_joules          # convert bin width to Joules
bin_midpoint_joules = bin_midpoint * MeV_to_joules    # convert bin midpoint to Joules

# 1) count spectrum: counts / bin width, now in joules - units: counts per Joule
count_spectrum = step_counts / bin_width_joules

# 2) energy spectrum: count spectrum * midpoint, now in joules
energy_spectrum = count_spectrum * bin_midpoint_joules

print("\n=== REDONE (Joule-consistent) COUNT SPECTRUM ===")
print(count_spectrum)
print("\n=== REDONE (Joule-consistent) ENERGY SPECTRUM ===")
print(energy_spectrum)

# recheck the peak location (energy in MeV still, for readability)
peak_idx = np.argmax(energy_spectrum)
peak_energy = bin_midpoint[peak_idx]
print(f"\nPeak still at: {peak_energy:.3e} MeV, ratio to kT: {peak_energy/kT_MeV:.3f}")

# STEP 4 (redone) - now everything is consistently in Joules...
detector_fraction = A_detector / (4 * 3.14159 * D_test**2)
recovered_spectrum_at_source = energy_spectrum / time_sample / detector_fraction
total_power_from_spectrum = np.sum(recovered_spectrum_at_source * bin_width_joules)

sb_power = power(M_test)

print("\n=== STEP 4 (redone, Joule-consistent) ===")
print(f"Recovered power from integrated spectrum: {total_power_from_spectrum:.3e} W")
print(f"Stefan-Boltzmann power(M):                 {sb_power:.3e} W")
print(f"Ratio (should be close to 1.0):             {total_power_from_spectrum / sb_power:.3f}")

peak_idx = np.argmax(energy_spectrum)
low_region_idx = np.arange(0, peak_idx // 2)   # first half of the way up to the peak, well below it
# checking for slpe here

log_E = np.log(bin_midpoint[low_region_idx])
log_spectrum = np.log(energy_spectrum[low_region_idx])

# fit a straight line to log(energy) vs log(spectrum) - the slope of that line
# is the power-law exponent (should be close to 2, per Rayleigh-Jeans)
slope, intercept = np.polyfit(log_E, log_spectrum, 1)

print(f"Number of points used (well below peak): {len(low_region_idx)}")
print(f"Fitted slope (should be close to 2.0): {slope:.4f}")

# also print a couple of individual point-to-point slopes as a sanity check
print("\nPoint-to-point slopes (consecutive bins, low-energy region):")
for i in range(len(low_region_idx) - 1):
    idx1 = low_region_idx[i]
    idx2 = low_region_idx[i + 1]
    local_slope = (np.log(energy_spectrum[idx2]) - np.log(energy_spectrum[idx1])) / \
                  (np.log(bin_midpoint[idx2]) - np.log(bin_midpoint[idx1]))
    print(f"  Bin {idx1}->{idx2}: slope = {local_slope:.4f}")