# primordial_blackholeevaporation

## What is this?

This repository models the evaporation of primordial black holes (PBHs)
via Hawking radiation. It derives the mass of a PBH that would be
completing its evaporation in the present epoch, builds an adaptive
numerical simulation of that evaporation, and predicts the photon
spectrum a real gamma-ray detector would observe, validating the
result against known blackbody physics (Stefan-Boltzmann law and the
Wien-like peak location).

There is also a paper (in `report/`) describing the physics and
methods in full, and presenting the final results.

Please note that the paper is still a work in progress!

## Getting the code

First, clone this repository to your own computer:

```bash
git clone https://github.com/archita-d/primordial_blackholeevaporation.git
cd primordial_blackholeevaporation
```

All commands in the rest of this README assume you are inside this
`primordial_blackholeevaporation` folder.

## Prerequisites

You will need Python 3 and the following packages:

```bash
pip install numpy matplotlib --break-system-packages
```

(Omit `--break-system-packages` if you are using a virtual environment.)

To build the paper you will also need a LaTeX distribution, for example:

```bash
sudo apt install texlive-latex-base texlive-latex-extra
```

## How to run the code and reproduce all results

All commands below should be run from the root of this repository.
Each script is self-contained and can be run on its own.

### 1. Core physics formulas
File: `code/formulacalculations.py`

This computes radius, temperature, and radiated power for a range of
black hole masses, as a sanity check on the underlying physics.

```bash
python3 code/formulacalculations.py
```

Expected output: printed radius, temperature, and power for several
test masses.

### 2. Deriving the crossover mass
File: `code/lifetimecomparison.py`

This derives the analytical mass-loss and lifetime formulas from
first principles, solves for the mass of a PBH that would be
completing its evaporation today, and produces a comparison plot
across three mass regimes.

```bash
python3 code/lifetimecomparison.py
```

Expected output: prints the crossover mass (approximately 1.73e11 kg)
and saves `three_mass_regimes.png` in `plots/`.

### 3. Adaptive evaporation simulation
File: `code/adaptivesimulation.py`

This runs a full adaptive time-step simulation of the crossover-mass
black hole from formation to complete evaporation, and cross-validates
the simulated lifetime against the analytical formula from step 2.

```bash
python3 code/adaptivesimulation.py
```

Expected output: prints the simulated lifetime (should closely match
the analytical value from step 2, within about 1%).

### 4. Detector pipeline
File: `code/slide8_photon.py`

This builds the four-step process for converting the blackbody photon
spectrum into a predicted 32-bin detector count array, for a given
mass, distance, and detector area. It also produces a bar-chart
visualization of the resulting spectrum.

```bash
python3 code/slide8_photon.py
```

Expected output: prints step-by-step intermediate results and the
final 32-element count array, and saves `gamma_ray_spectrum.png`
and `photon_function_integration.png` in `plots/`.

### 5. Peak energy and time-step tracking
File: `code/dtpeakenergy.py`

This tracks the adaptive time step size and the black hole's peak
emission energy as functions of time leading up to full evaporation.

```bash
python3 code/dtpeakenergy.py
```

Expected output: saves `dt_vs_time.png` and `peak_energy_vs_time.png`
in `plots/`, and prints the raw dt and peak-energy values.

### 6. Blackbody validation
File: `code/validating.py`

This is the key validation of the entire pipeline. It converts the
32-bin detector counts into a count spectrum and an energy spectrum,
checks that the energy spectrum peaks at the expected location
(approximately 3 times kT), checks that the spectrum's slope below
the peak matches the expected Rayleigh-Jeans value of 2.0, and
checks that the total power recovered from the spectrum matches the
independent Stefan-Boltzmann formula.

```bash
python3 code/validating.py
```

Expected output:
- Energy spectrum peak at approximately 2.99 times kT
- Slope below the peak of approximately 2.0
- Recovered-power to Stefan-Boltzmann ratio of approximately 0.989
- Saves `step_counts_debug.png` in `plots/`, showing the full
  blackbody-shaped spectrum

## Key results to expect

- Crossover mass (evaporating in the present epoch): approximately 1.73e11 kg
- Cross-validated lifetime (analytical vs. simulated): agreement within approximately 1%
- Blackbody validation: peak at 2.99 x kT, slope of 2.0, power ratio of 0.989
- Peak emission energy climbs from approximately 100 MeV to over 1e19 MeV over the final moments of evaporation

## Known limitations

This model treats the black hole as an ideal blackbody (Stefan-Boltzmann
law) and does not include greybody factors or the emission of particle
species other than photons (for example, electron-positron pairs, which
become relevant below approximately 2e13 kg).

Please note that this project is still a work in progress and additional sections will be added shortly.

## Acknowledgements
I would like to thank my mentor, Ed Fenimore of Los Alamos National Laboratory, for his guidance and feedback throughout this project. I would also like to thank the Institute for Computing inResearch for organizing this program.

Some short snippets of code (<5 lines) were written with the assistance of generative AI.
