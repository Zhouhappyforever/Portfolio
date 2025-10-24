# 2D Reservoir Simulation with Nechelik Data
This repository contains a Python-based simulator for two-dimensional reservoir pressure modeling using finite difference methods. It demonstrates my skills in numerical methods for petroleum engineering, including implicit solvers for transient flow problems in heterogeneous reservoirs with multiple wells. This project was developed as part of my coursework and research at UT Austin, building on reservoir simulation techniques for applications in hydraulic fracturing and enhanced oil recovery.
## Project Overview
The simulator models pressure distribution in a 2D reservoir discretized into grid blocks, accounting for heterogeneous permeability and porosity from real data (Nechelik reservoir), variable grid spacing, boundary conditions (prescribed pressure or flux), and multiple wells (rate or BHP controlled). It uses sparse matrices for efficiency, handles inactive cells, and supports contour plotting of pressure profiles.
Key features:
- Reads permeability and porosity from data files for realistic heterogeneous simulations.
- Supports implicit solvers with well models including skin factors and Peaceman's equation.
- Uses SciPy for sparse linear algebra and conjugate gradient solvers.
- Includes visualization with Matplotlib for pressure contours.
- Tested with unit tests for transmissibility, accumulation, solver accuracy, and well integrations.
This aligns with my expertise in reservoir simulation, as highlighted in my portfolio, where I've implemented similar models in C++ and Julia for multiphase flow and geomechanics.
## Installation
1. Create the Conda environment: `conda env create -f environment.yml`
2. Activate the environment: `conda activate twodreservoir`
## Usage
- The core classes are in `two_d_reservoir_simulator.py` and `project1.ipynb`.
- Example: Load inputs from a YAML file and run the simulation:
  ```python
  import matplotlib.pyplot as plt
  from two_d_reservoir_simulator import TwoDimReservoir
  from project1 import Project1
  simulator = Project1('inputs.yml')
  simulator.solve()
  simulator.plot()
  plt.show()

Customize inputs in inputs.yml (e.g., wells, boundaries, solver type).
Data files: Nechelik_perm.dat and Nechelik_poro.dat for permeability and porosity.
Run unit tests: 

python unit_tests.py

## Files

two_d_reservoir_simulator.py: Base 2D simulation class.
two_d_reservoir_simulator.ipynb: Jupyter notebook for 2D class development.
project_simulator.py: Project-specific class with file reading and plotting.
project_simulator.ipynb: Jupyter notebook for the main project execution and visualization.
environment.yml: Conda environment specification.
unit_tests.py: Unit tests for validation.
inputs.yml: Input configuration file.
Nechelik_perm.dat, Nechelik_poro.dat: Reservoir property data files.

## Example Results
After running the simulation with Nechelik data, pressure contours show flow from injectors to producers, influenced by heterogeneous permeability. Inactive cells are masked, and pressures stabilize over time steps.
## Contact
For questions or collaborations, reach me at andrewmiller@utexas.edu.
Andrew Z. Miller
Master's Student in Petroleum Engineering, UT Austin