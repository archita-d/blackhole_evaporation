import numpy as np
import matplotlib.pyplot as plt

# constants
G = 6.67e-11
c = 3e8
hbar = 1.055e-34   # h bar, aka plancks constant/2pi
k = 1.381e-23
MeV_to_joules = 1.602176634e-13

def radius(M):
    return 2 * G * M / c**2

def temperature(M):
    return hbar * c**3 / (8 * 3.14159 * G * M * k)

def b_function(E_MeV, T): #aka black body photon function
    E = E_MeV * MeV_to_joules
    x = E / (k * T)
    x = np.minimum(x, 700)
    denominator = np.expm1(x)
    denominator = np.where(denominator == 0, 1e-300, denominator)
    return (2 * E**2) / (h**3 * c**2) / denominator

# this is the detector pipeline,  steps 1-4
# takes the blackbody function and turns it into actual counts a
# detector would see

# step one from slide eight where we multiplying by the sq-cm of the blackbody surface
def step1(M, E):
    T = temperature(M)
    r_s = radius(M)
    surface_area_m2 = 4 * 3.14159 * r_s**2
    p_function_value = b_function(E, T)
    return p_function_value * surface_area_m2

# step 2 is to multiply by the time sample of the gamma-ray detectr (1 microsec)
def step2(M, E, time_sample):
    step1_result = step1(M, E)
    return step1_result * time_sample

# step 3 is to integrate the blackbody function for each of the energy ranges
def step3(M, E_low, E_high, time_sample):
    T = temperature(M)
    r_s = radius(M)
    surface_area_m2 = 4 * 3.14159 * r_s**2

    E_values = np.linspace(E_low, E_high, 200) # 200 is the amount of slices made
    p_function_values = b_function(E_values, T)
    E_values_joules = E_values * MeV_to_joules  # gotta integrate in joules not mev
    integrated_range = np.trapezoid(p_function_values, E_values_joules)

    gamma_rays = integrated_range * surface_area_m2 * time_sample * 3.14159

    return gamma_rays

# step 4 is to multiply by how much area the detector subtends at distance D
def step4(M, E_low, E_high, time_sample, D, A_detector):
    gamma_rays = step3(M, E_low, E_high, time_sample)
    detector = A_detector / (4 * 3.14159 * D**2)
    return gamma_rays * detector

# the energy bin array
energy_bin = np.logspace(np.log10(0.001), np.log10(1.0), 33)  # 0.001 mev = 1 kev, 33 edges

def gamma_array(M, energy_bin, time_sample, D, A_detector):
    counts = []
    for i in range(len(energy_bin) - 1):
        E_low = energy_bin[i]
        E_high = energy_bin[i + 1]
        count = step4(M, E_low, E_high, time_sample, D, A_detector)
        counts.append(count)
    return np.array(counts)


# tesing it out

M_test = 1e11
E_test = 0.5          # mev, placeholder
E_low, E_high = 0.25, 0.75
time_sample = 1e-6
AUconversion = 1.495978707e11   # 1 au in METERS now
D_test = 100 * AUconversion
A_detector = 50e-4   # 50 cm^2 into  m^2

result_step2 = step2(M_test, E_test, time_sample)
print("Step 2 result:", result_step2)

result_step3 = step3(M_test, E_low, E_high, time_sample)
print("Step 3 result:", result_step3)

result_step4 = step4(M_test, E_low, E_high, time_sample, D_test, A_detector)
print("Step 4 (detected) result:", result_step4)

gamma_array_result = gamma_array(M_test, energy_bin, time_sample, D_test, A_detector)
print("Gamma-ray array (32 elements):")
print(gamma_array_result)


# just a bar chart to see what the spectrum looks like
def plot_spectrum(energy_bin, gamma_array_result):
    bincenters = (energy_bin[:-1] + energy_bin[1:]) / 2
    binwidths = np.diff(energy_bin)

    plt.figure(figsize=(10, 6))
    plt.bar(bincenters, gamma_array_result, width=binwidths, color='blue', alpha=0.6, edgecolor='black')
    plt.xlabel("Photon Energy (MeV)")
    plt.ylabel("Detected Gamma-Rays")
    plt.title("Detected Gamma-Ray Counts Across Energy Bins")
    plt.xscale('log')  # log scales show more detail here
    plt.yscale('log')
    plt.grid(True, alpha=0.5)
    plt.savefig("gamma_ray_spectrum.png")
    plt.show()

plot_spectrum(energy_bin, gamma_array_result)


# plotting the photon function curve with the integrated bin shaded in
# so u can actually see what the integration in step3 is doing
def plot_integration(M, E_low, E_high, plot_range_low=0.01, plot_range_high=2.0):
    T = temperature(M)

    E_full = np.linspace(plot_range_low, plot_range_high, 500)
    photon_full = b_function(E_full, T)

    E_bin = np.linspace(E_low, E_high, 200)
    photon_bin = b_function(E_bin, T)

    plt.figure(figsize=(8, 6))
    plt.plot(E_full, photon_full, color='blue', label='Blackbody photon function')
    plt.fill_between(E_bin, photon_bin, color='blue', alpha=0.3,
                      label=f'Integrated: {E_low}-{E_high} MeV')
    plt.xlabel("Photon Energy (MeV)")
    plt.ylabel("Photons / sec / MeV / steradian / sq-m")
    plt.title(f"Blackbody Photon Function at T = {T:.3e} K")
    plt.legend()
    plt.grid(True, alpha=0.5)
    plt.savefig("photon_function_integration.png")
    plt.show()

plot_integration(M_test, E_low, E_high)