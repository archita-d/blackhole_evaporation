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

C_const = power(1.0)

def mass_at_time_remaining(t_remaining):
    # this basically lets us jump stright to the mass at any specific time we want before
    if t_remaining <= 0:
        return 0.0
    return (3 * C_const * t_remaining / c**2) ** (1/3)

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
    surface_area_m2 = 4 * 3.14159 * r_s**2
    E_values = np.linspace(E_low, E_high, 200)
    p_values = b_function(E_values, T)
    E_values_joules = E_values * MeV_to_joules
    integrated = np.trapezoid(p_values, E_values_joules)
    return integrated * surface_area_m2 * time_sample * 3.14159

def step4(M, E_low, E_high, time_sample, D, A_detector):
    emitted = step3(M, E_low, E_high, time_sample)
    detector_fraction = A_detector / (4 * 3.14159 * D**2)
    return emitted * detector_fraction


# so we are going to zoom into a specific window right before evaporation
# you can also adjust window_seconds to zoom in/out if you want

signal_bands = {
    "1e9-2e9 MeV":   (1e9, 2e9),
    "5e9-6e9 MeV":   (5e9, 6e9),
    "1e10-2e10 MeV": (1e10, 2e10),
    "5e10-6e10 MeV": (5e10, 6e10),
    "1e11-2e11 MeV": (1e11, 2e11),
}

D_test = 9.4607e15      # 1 light-year, in meters
A_detector = 50e-4      # 50 cm^2, in m^2
time_sample = 1e-6    # 1 microsecond

def plot_zoomed_window(window_seconds, n_points=500, margin_fraction=0.1):
    margin = window_seconds * margin_fraction
    t_axis = np.linspace(-window_seconds, margin, n_points)

    signals = {name: [] for name in signal_bands}
    temps = []
    masses = []

    for t in t_axis:
        t_remaining = -t
        M = mass_at_time_remaining(t_remaining)

        if M <= 0:
            temps.append(np.nan)
            masses.append(np.nan)
            for name in signal_bands:
                signals[name].append(0.0)
            continue

        T = temperature(M)
        temps.append(T)
        masses.append(M)

        for name, (E_low, E_high) in signal_bands.items():
            count = step4(M, E_low, E_high, time_sample, D_test, A_detector)
            rate = count / time_sample
            signals[name].append(rate)

    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t_axis, temps, color='darkred')
    axes[0].set_ylabel("Temperature (K)")
    axes[0].set_title(f"Final {window_seconds:.1e} s before evaporation "
                       f"(mass shrinks to ~{masses[0]:.2e} kg by start of window)")
    axes[0].grid(True, alpha=0.4)

    for name, values in signals.items():
        axes[1].plot(t_axis, values, marker='o', markersize=3, label=name)
    axes[1].set_xlabel("Time before evaporation (s)  [0 = evaporation, on the right]")
    axes[1].set_ylabel("Detected Rate (counts/sec)")
    axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig(f"zoomed_window_{window_seconds:.0e}s.png", dpi=200)
    plt.show()

plot_zoomed_window(window_seconds=1e-5)
