# Primordial Black Hole Evaporation
This repo models Hawking radiation and black hole evaporation over time 
## Files
### Code
#### `formulacalculations`, `timestepmass_nonadaptive`, `massradiustempplots`, `microsecondmass`, `microsecondtemp`,`microsecondpower`
All of these files showcase core Hawking radiation formulas (radius, temperature, power) as well as fixed-timestep simulations of mass/temperature/power, at 1-second and microsecond resolution, mostly done for testing purposes. 

#### `adaptivetemp_test`
This file successfully simulates temperature and mass loss for a 1e6 kg black hole using an adaptive time step (`dt = 0.001 * M * c**2 / P`) instead of a fixed one, so dt is automatically shrinking as the black hole nears the end of its life instead of using the previous constant microsecond step.

The original microsecond simulation used a constant `dt = 1e-6 s` for every step, which wastes steps early on when it was not needed (the black hole is barely changing). I was able to recompute `dt` each step based on the current mass and power which allows the simulation take large steps early (when change is slow) and small steps near the end (when change is fast). 

#### `lifetime_comparison`
this file derives a closed-form formula for black hole mass as a function of time directly from the Hawking radiation law, and then uses it to compare evaporation across a wide range of masses without simulating each one.

Since power `P ∝ 1/M²` (from combining area `∝ M²` and temperature⁴ `∝ 1/M⁴`) -  came from `formulacalculations`, and `dM/dt = -P/c²`, solving this rate equation gives:
M(t) = cube_root(M₀³ - 3k·t)        where k = C/c², C = power(1.0) 
note: C gives us all the constants from the power function (see `formulacalculations`)

and total lifetime:
t_life = M₀³ · c² / (3C)

then we can rearrange for the mass that evaporates exactly today (age of universe ≈ 13.8 Gyr):
M₀ = cube_root(3·C·t_life / c²) = 1.73e11 kg

then we plot `mass_vs_time` for three regimes (1e8 kg, 1.73e11 kg, 1e20 kg) over the universe's entire age, to show why only the medium mass is evaporating now and as shown the smalll ones vanished almost instantly, large ones have barely lost mass.

### Plots
#### `mass_vs_time.png`, `massvradius.png`, `temperature_plot.png`, `mass_vs_time_micro_full.png`, `temperature_vs_time_micro_full.png`
Plots of mass, radius, and temperature vs. mass or time for a 1e6 kg black hole, at both 1-second and microsecond resolution.
