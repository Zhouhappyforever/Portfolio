# 1D Reservoir Simulation

This repository contains a Python-based simulator for one-dimensional reservoir pressure modeling using finite difference methods. It demonstrates my skills in numerical methods for petroleum engineering, including implicit, explicit, and mixed solvers for transient flow problems. This project was developed as part of my coursework and research at UT Austin, building on reservoir simulation techniques for applications in hydraulic fracturing and enhanced oil recovery.

## Project Overview
The simulator models pressure distribution in a 1D reservoir discretized into grid blocks, accounting for variable permeability, porosity, grid spacing, and boundary conditions (prescribed pressure or flux). It uses sparse matrices for efficiency and supports plotting of pressure profiles over time.

Key features:
- Handles homogeneous and heterogeneous reservoir properties.
- Supports explicit, implicit, and theta-mixed methods for time-stepping.
- Uses SciPy for sparse linear algebra and conjugate gradient solvers.
- Tested with unit tests for transmissibility, accumulation, and solver accuracy.

This aligns with my expertise in reservoir simulation, as highlighted in my portfolio, where I've implemented similar models in C++ and Julia for multiphase flow and geomechanics.

## Installation
1. Create the Conda environment: `conda env create -f environment.yml`
2. Activate the environment: `conda activate onedreservoir`

## Usage
- The core class is in `one_d_reservoir_simulator.py`.
- Example: Load inputs from a YAML file and run the simulation:
  `import matplotlib.pyplot as plt`
  `from one_d_reservoir_simulator import OneDimReservoir`
  `simulator = OneDimReservoir('inputs.yml')`
  `simulator.solve()`
  `simulator.plot()`
  `plt.show()`
- Customize inputs in a YAML file (e.g., boundary conditions, solver type).
- Run unit tests: `python unit_tests.py`

## Files
- `one_d_reservoir_simulator.py`: Main simulation class.
- `one_d_reservoir_simulator.ipynb`: Jupyter notebook for interactive execution and visualization.
- `environment.yml`: Conda environment specification.
- `unit_tests.py`: Unit tests for validation.
- `LICENSE`: Apache License 2.0.

## Example Results
After running the simulation with default inputs, pressure profiles evolve as fluid flows from high-pressure boundaries, stabilizing over time steps. Heterogeneous cases show varied transmissibility impacts.

## Contact
For questions or collaborations, reach me at andrewmiller@utexas.edu.

Andrew Z. Miller  
Master's Student in Petroleum Engineering, UT Austin  
LinkedIn: [Andrew Z. Miller](https://www.linkedin.com/in/andrew-z-miller)