import numpy as np                 
import matplotlib.pyplot as plt     

G = 6.67e-11                           # gravitational constant
c = 3e8                                # speed of light
h = 6.626e-34                          # planck's constant
k = 1.381e-23                          # boltzmann's constant

def radius(M):
    return 2 * G * M / c**2            # schwarzschild radius

def temperature(M):
    return (h / (2 * 3.14159)) * c**3 / (8 * 3.14159 * G * M * k)
    # hawking temp 
def power(M):
    r_s = radius(M)                   
    A = 4 * 3.14159 * r_s**2           
    T = temperature(M)                 
    return 5.670e-8 * A * T**4        

MeV_to_joules = 1.602176634e-13        
def b_function(E_MeV, T):
    E = E_MeV * MeV_to_joules          # convert to joules
    x = E / (k*T)                      # ratio of photon energy to thermal energy
    x = np.minimum(x, 700)             # cap it so it never ovrflows
    denominator = np.expm1(x)     
    denominator = np.where(denominator == 0, 1e-300, denominator)
   
    return (2 * E**2) / (h**3 * c**2) / denominator
    # blackbody photon-count spectrum at this energy/temp

def peak_energy_MeV(T, n_points=2000):
    # the peak of a blackbody spectrum sits at roughly 3x kT so we cn use that to set a generous search window
    kT_MeV = (k * T) / MeV_to_joules   
    search_low = kT_MeV * 0.01        
    search_high = kT_MeV * 20          
    
    E_values = np.logspace(np.log10(search_low), np.log10(search_high), n_points)
    # 2000 sample points across that range
    p_values = b_function(E_values, T)   
    peak_idx = np.argmax(p_values)       
    return E_values[peak_idx]            

def run_dt_and_peak_tracking(M0, max_steps=100000, mass_cutoff=1e-6, save_every=50, xraythres=1e6):
    M = M0                              # starting mass
    time_elapsed = 0                    # clock starts at zero

    history_time = []                   # will hold saved time snapshots
    history_dt = []                     # will hold saved dt values
    history_mass = []                   # will hold saved mass values
    history_temp = []                   # will hold saved temp values
    history_peak_energy = []            # will hold saved peak-energy values

    for i in range(max_steps):         
        P = power(M)                
        dt = 0.001 * M * c**2 / P      
        T = temperature(M)            

        if M < mass_cutoff:             # basically fully evapoated
            print("Mass reached 0 at step", i, "time:", time_elapsed, "s")
            break

        dM = (P * dt) / c**2            # mass lost this step
        M = M - dM                      # update mass
        time_elapsed += dt              # update clock

        if i % save_every == 0:       
            history_time.append(time_elapsed)
            history_dt.append(dt)
            history_mass.append(M)
            history_temp.append(T)
            if T >= xraythres:           # only bothr finding the peak if its hot enough to mater
                history_peak_energy.append(peak_energy_MeV(T))
            else:
                history_peak_energy.append(np.nan)   # not meaningful before it's hot enough
               

    return (np.array(history_time), np.array(history_dt), np.array(history_mass),
            np.array(history_temp), np.array(history_peak_energy))
    # hand back all 5 histories


M0_run = 1.73e11                       
ht, hdt, hmass, htemp, hpeak = run_dt_and_peak_tracking(M0_run)
# run the simulation and unpack all 5 returned arrays

# Plot 1 is dt vs time
plt.figure(figsize=(9, 6))
plt.plot(ht, hdt, marker='o', linestyle='none', color='darkgreen')
# points onlye
plt.xlabel("Time (s)")
plt.ylabel("dt (s)")
plt.title(f"Adaptive Time Step (dt) vs. Time (M0 = {M0_run:.2e} kg)")
plt.xscale('log')                      
plt.yscale('log')
plt.grid(True, which='both', alpha=0.4)
plt.tight_layout()
plt.savefig("dt_vs_time.png", dpi=200)
plt.show()

# Plot 2: peak emission energy vs time
plt.figure(figsize=(9, 6))
plt.plot(ht, hpeak, marker='o', linestyle='none', color='darkorange')
plt.xlabel("Time (s)")
plt.ylabel("Peak Emission Energy (MeV)")
plt.title(f"Peak Emission Energy vs. Time (M0 = {M0_run:.2e} kg)")
plt.xscale('log')
plt.yscale('log')
plt.grid(True, which='both', alpha=0.4)
plt.tight_layout()
plt.savefig("peak_energy_vs_time.png", dpi=200)
plt.show()

print("dt values:", hdt)                 # list the raw dt numbers to check them
print("Peak energy values (MeV):", hpeak) # list the raw peak-enrgy numbers to
