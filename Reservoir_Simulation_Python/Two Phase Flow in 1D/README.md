# Buckley-Leverett Analytical Solution

This repository contains a Python implementation of the Buckley-Leverett semi-analytical solution for modeling two-phase (water-oil) displacement in reservoirs. It demonstrates my expertise in multiphase flow dynamics, relative permeability modeling using Corey-Brooks correlations, and numerical optimization for shock front saturation. Developed as part of my coursework and research at UT Austin, this tool validates numerical simulators for applications in enhanced oil recovery and reservoir management.

## Project Overview
The code computes fractional flow, relative permeabilities, and saturation profiles, solving for the shock front saturation using nonlinear optimization. It supports visualization of fractional flow curves and saturation profiles at specified times, handling variable viscosities, residual saturations, and exponents.

Key features:
- Corey-Brooks model for water and oil relative permeabilities.
- Fractional flow and derivative calculations.
- Nonlinear solve for shock front saturation using SciPy.
- Plotting functions for fractional flow and saturation profiles (static and interactive).
- Tested for accuracy in front saturation under homogeneous and varied conditions.

This project complements my portfolio in reservoir simulation, where I've developed finite difference-based multiphase models in Python and MATLAB for CO2 sequestration and hydraulic fracturing studies.

## Installation
1. Create the Conda environment: `conda env create -f environment.yml`
2. Activate the environment: `conda activate buckley_leverett`

## Usage
- The core class is in `buckley_leverett.py`.
- Example: Load inputs from a YAML file and plot a saturation profile:
  ```python
  import matplotlib.pyplot as plt
  from buckley_leverett import BuckleyLeverett
  bl = BuckleyLeverett('inputs.yml')
  bl.plot_saturation_profile(0.25)
  plt.show()