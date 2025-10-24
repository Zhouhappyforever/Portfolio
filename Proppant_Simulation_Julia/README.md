# Proppant Distribution Simulation in Hydraulic Fracturing

## Overview

This repository contains a Julia-based simulation for modeling two-phase proppant transport in a horizontal pipe, a critical component in hydraulic fracturing operations. Developed as part of my research in reservoir simulation and geomechanics at The University of Texas at Austin, this interactive demonstration leverages finite element methods to predict proppant distribution, fluid-particle interactions, and pressure dynamics. The simulation is built using the Gridap.jl package for high-fidelity numerical solutions, showcasing my expertise in computational modeling for petroleum engineering applications.

As a Master's student in Petroleum Engineering with a focus on hydraulic fracturing, reservoir simulation, and AI-driven optimization (as detailed in my [resume](Resume_Andrew.docx)), this project demonstrates my ability to implement advanced numerical schemes for multiphase flow problems. It aligns with my research under Dr. Foster on neural networks for proppant distribution modeling and my prior experience developing C++ hydraulic fracturing simulations, where I improved computational efficiency by 16%. This Julia implementation serves as a flexible, open-source tool for exploring proppant behavior, which is essential for optimizing fracturing designs and enhancing well productivity.

## Key Features

- **Two-Phase Flow Modeling**: Simulates particle (proppant) and fluid phases using mass and momentum conservation equations, incorporating drag forces, collision dispersive pressure, and virtual mass effects.
- **Finite Element Discretization**: Utilizes Gridap.jl for 1D domain meshing, multi-field finite element spaces, and transient solving with backward Euler time stepping.
- **Stability Enhancements**: Includes skew-symmetric forms for mass balance and repulsive force corrections to handle high particle fractions and prevent instabilities.
- **Interactive Visualization**: In-line plotting with Plots.jl for particle fraction, velocities, and lambda at each time step, plus VTK export for advanced post-processing.
- **Customizable Parameters**: Easily adjustable domain length, material properties, time steps, and boundary conditions for sensitivity analysis in fracturing scenarios.
- **Debugging and Validation**: Built-in print statements and type checks for robust development, ensuring reliable results in production environments.

This simulation is particularly relevant for hydraulic fracturing positions, as it addresses proppant settling, transport efficiency, and pressure gradients—key factors in fracture propagation and reservoir stimulation.

## Requirements

- Julia v1.6 or higher
- Required packages (install via Julia's Pkg manager):
  - Gridap.jl
  - Gridap.Geometry
  - Gridap.FESpaces
  - Gridap.MultiField
  - Gridap.ODEs
  - Gridap.CellData
  - Gridap.Visualization
  - Plots.jl

To install dependencies, run the following in Julia:
```julia
using Pkg
Pkg.add(["Gridap", "Plots"])
```

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/proppant-simulation.git
   cd proppant-simulation
   ```

2. Run the simulation script:
   ```bash
   julia Simulation_2Phase.jl
   ```

3. Customize parameters in Section 2 of the script (e.g., domain length `L`, particle diameter `d_part`, inlet conditions) to simulate different fracturing scenarios.

4. Outputs:
   - Console logs for debugging and progress.
   - `final_plot.png`: Side-by-side plots of particle fraction, velocities, and lambda at the final time step.
   - `final_solution.vtk`: Export for visualization in tools like Paraview.

Example run: The default parameters simulate a 0.1 m horizontal pipe with 100 elements over 0.01 seconds, demonstrating proppant influx at the inlet and transport dynamics.

## Results and Interpretation

The simulation captures the evolution of proppant distribution, showing how particles propagate under fluid drag while accounting for gravity and collisions. Sample outputs from a default run:

- **Particle Fraction (φ_p)**: Increases near the inlet due to specified boundary conditions, illustrating proppant placement in fractures.
- **Velocities (v_p, v_f)**: Fluid velocity drives particle movement, with relative slip influenced by drag coefficients.
- **Lambda (λ)**: Represents pressure-like terms in the momentum equations, highlighting stress distributions.

These results can inform fracturing job designs by predicting proppant settling risks and optimizing slurry compositions. For advanced use, integrate with machine learning models (e.g., neural networks) to calibrate parameters against field data, as explored in my UT Austin research.

## Limitations and Future Work

- Current model assumes 1D flow in a horizontal pipe; extensions to 2D/3D fractures are planned.
- Simplified constants (e.g., drag and bpress) for testing; variable formulations are implemented but can be toggled for realism.
- Potential integration with AI: Use as a physics-based baseline for generative AI uncertainty quantification, building on my experience with OpenAI APIs and ML certifications (Stanford Machine Learning).

## Acknowledgments

This work draws from my research at UT Austin under Drs. Foster, Daigle, and Sharma, as well as practical experience at Beijing Ennosoft Corp. in reservoir simulation algorithms. Inspired by literature on multiphase flow in porous media and hydraulic fracturing mechanics.

## Contact

For collaborations or inquiries related to hydraulic fracturing simulations, reach out at andrewmiller@utexas.edu or connect on LinkedIn. Open to discussions on integrating this with high-performance computing or AI-enhanced modeling.

## License

Feel free to use and modify for academic or professional purposes.