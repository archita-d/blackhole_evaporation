import numpy as np
import matplotlib.pyplot as plt

# constants
G = 6.67e-11
c = 3e8
h = 6.626e-34
hbar = h / (2 * 3.14159)
k = 1.381e-23
sigma = 5.670e-8
MeV_to_joules = 1.602176634e-13

def radius(M):
    return 2 * G * M / c**2

def temperature(M):
    return hbar * c**3 / (8 * 3.14159 * G * M * k)

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
    surface_area_m2 = 4 * 3.14159 * r_s**2       # stays in m^2 to  matches b_function's uni
    E_values = np.linspace(E_low, E_high, 200)
    p_values = b_function(E_values, T)
    E_values_joules = E_values * MeV_to_joules    # integrated in Joules now
    integrated = np.trapezoid(p_values, E_values_joules)
    return integrated * surface_area_m2 * time_sample * 3.14159   

def step4(M, E_low, E_high, time_sample, D, A_detector):
    emitted = step3(M, E_low, E_high, time_sample)
    detector_fraction = A_detector / (4 * 3.14159 * D**2)
    return emitted * detector_fraction

def gamma_array(M, energy_bin, time_sample, D, A_detector):
    counts = []
    for i in range(len(energy_bin) - 1):
        counts.append(step4(M, energy_bin[i], energy_bin[i + 1], time_sample, D, A_detector))
    return np.array(counts)

# blackbody validation will first convert counts into a count spectrum
# and energy spectrum, then
# check the peak locaton, the slope, and the total power
# and all of this is to validate our results

M_test = 1.73e11        # test mass
D_test = 9.4607e15      # 1 light-year, in meters
A_detector = 50e-4      # 50 cm^2, in m^2
time_sample = 1e-6      # 1 microsecond

energy_bin = np.logspace(np.log10(0.001), np.log10(500.0), 33)
E_low = energy_bin[:-1]
E_high = energy_bin[1:]
bin_midpoint = (E_low + E_high) / 2

bin_width_joules = (E_high - E_low) * MeV_to_joules
bin_midpoint_joules = bin_midpoint * MeV_to_joules

step_counts = gamma_array(M_test, energy_bin, time_sample, D_test, A_detector)
count_spectrum = step_counts / bin_width_joules
energy_spectrum = count_spectrum * bin_midpoint_joules

print("Step counts (32 bins):")
print(step_counts)
print("\nCount spectrum (per Joule):")
print(count_spectrum)
print("\nEnergy spectrum:")
print(energy_spectrum)

# check 1 is the pek location vs kT 
T_test = temperature(M_test)
kT_MeV = (k * T_test) / MeV_to_joules
peak_idx = np.argmax(energy_spectrum)
peak_energy = bin_midpoint[peak_idx]
print(f"\nTemperature: {T_test:.3e} K")
print(f"kT: {kT_MeV:.3e} MeV")
print(f"Energy spectrum peaks at: {peak_energy:.3e} MeV")
print(f"Ratio of peak energy to kT: {peak_energy / kT_MeV:.3f}  (expected is around 3)")

# check 2 is the slope below the peak (expecte is around 2.0) 
low_region_idx = np.arange(0, peak_idx // 2)
log_E = np.log(bin_midpoint[low_region_idx])
log_spectrum = np.log(energy_spectrum[low_region_idx])
slope, _ = np.polyfit(log_E, log_spectrum, 1)
print(f"Slope below the peak: {slope:.4f}  (expected: ~2.0)")

# check 3 is the our total power vs. Stefan-Boltzmann
detector_fraction = A_detector / (4 * 3.14159 * D_test**2)
raw_counts_at_source = step_counts / detector_fraction / time_sample
energy_per_bin = raw_counts_at_source * bin_midpoint_joules
total_power_from_spectrum = np.sum(energy_per_bin * bin_width_joules)
sb_power = power(M_test)
print(f"Recovered power: {total_power_from_spectrum:.3e} W")
print(f"Stefan-Boltzmann power: {sb_power:.3e} W")
print(f"Ratio: {total_power_from_spectrum / sb_power:.3f}  (expected: ~1.0)")

# plot the validated specrum
plt.figure(figsize=(9, 6))
plt.plot(bin_midpoint, energy_spectrum, marker='o')
plt.xlabel("Energy (MeV)")
plt.ylabel("Energy Spectrum")
plt.title(f"Blackbody Validation (M = {M_test:.2e} kg)")
plt.xscale('log')
plt.yscale('log')
plt.grid(True, which='both', alpha=0.4)
plt.tight_layout()
plt.savefig("step_counts_debug.png", dpi=200)
plt.show()
