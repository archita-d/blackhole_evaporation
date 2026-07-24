import numpy as np

# constants 
G = 6.67e-11          # gravitational constant, m^3/(kg s^2)
c = 3e8                # speed of light, m/s
hbar = 1.055e-34   # h bar, aka plancks constant
k = 1.381e-23          # boltzanns constant, J/K
sigma = 5.670e-8       # stefan boltzmann constant
MeV_to_joules = 1.602176634e-13   # need this bc our bins are in MeV but the constants are all in joules

def radius(M):
    # schwarzchild radius
    return 2 * G * M / c**2

def temperature(M):
    # hawking temp 
    return hbar * c**3 / (8 * 3.14159 * G * M * k)

def power(M):
    # stefan boltzmann law which calclates total power radiated
    r_s = radius(M)
    A = 4 * 3.14159 * r_s**2
    T = temperature(M)
    return sigma * A * T**4

# testing these formulas
MIN_MASS = 1e3     # kg
MAX_MASS = 1e25    # kg

M = float(input(f"Enter a black hole mass in kg (between {MIN_MASS:.0e} and {MAX_MASS:.0e}): "))

if M < MIN_MASS or M > MAX_MASS:
    print(f"Mass must be between {MIN_MASS:.0e} and {MAX_MASS:.0e} kg")
else:
    r = radius(M)
    T = temperature(M)
    P = power(M)
    print(f"\nMass = {M:.3e} kg")
    print(f"  Radius:      {r:.3e} m")
    print(f"  Temperature: {T:.3e} K")
    print(f"  Power:       {P:.3e} W")
