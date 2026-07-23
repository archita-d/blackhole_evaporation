import numpy as np
import matplotlib.pyplot as plt
G = 6.67e-11
c = 3e8
hbar = 1.055e-34
k = 1.381e-23
sigma = 5.670e-8
def radius(M):
    return 2 * G * M / c**2
def temperature(M):
    return hbar * c**3 / (8 * 3.14159 * G * M * k)
def power(M):
    r_s = radius(M)
    A = 4 * 3.14159 * r_s**2
    T = temperature(M)
    return sigma * A * T**4
C = power(1.0) # this is for storing the constants
def mass_vs_time(M0, t):
    inside = M0**3 - (3 * C / c**2) * t
    M = np.where(inside > 0, np.cbrt(np.maximum(inside, 0)), 0)
    return M
age_universe_s = 13.8e9 * 365.25 * 24 * 3600
t_years = np.linspace(0, 13.8e9, 1000)
t_seconds = t_years * 365.25 * 24 * 3600

# attempting to solve for a mass that evaporates today                e
M0_today = (3 * C * age_universe_s / c**2) ** (1/3)
print("Mass that evaporates today:", M0_today, "kg")

masses_to_test = {
    "Small (1e8 kg) - evaporated long ago": 1e8,
   f"Medium ({M0_today:.2e} kg) - evaporating now": M0_today,
    "Large (1e20 kg) - barely lost mass": 1e20,
}
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=True, sharey=False)
for ax, (label, M0) in zip(axes, masses_to_test.items()):
    M_t = mass_vs_time(M0, t_seconds)
    ax.plot(t_years, M_t / M0, color='blue')
    ax.set_title(label)
    ax.set_xlabel("Time since Big Bang (years)")
    ax.set_ylabel("Fraction of the mass remaining")
    ax.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig("three_mass_regimes.png")
plt.show()