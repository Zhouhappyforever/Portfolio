# Andrew Z. Miller's Petroleum Engineering Portfolio

Welcome to my GitHub repository! This portfolio showcases my skills and projects in petroleum engineering, with a focus on numerical simulations for reservoir modeling, hydraulic fracturing, and sustainable energy systems. As a Master's student in Petroleum Engineering at UT Austin, I have hands-on experience developing simulation codes in Julia, Python, C++, and MATLAB, optimizing algorithms for high-performance computing, and applying them to real-world challenges like proppant distribution in hydraulic fractures and carbon management.

My work aligns with industry needs in hydraulic fracturing and reservoir simulation, including code development for multiphase flow models and techno-economic analysis. For more details on my background, please refer to my resume.

## Key Skills
- **Programming & Tools**: Julia (for high-performance simulations), Python (with PyTorch for ML integration), C++ (for reservoir algorithms), MATLAB.
- **Technical Expertise**: Reservoir simulation (black oil and multiphase models), hydraulic fracturing modeling (proppant transport, finite element methods), geomechanics, petrophysics, techno-economic analysis for CCUS/CDR.
- **Software & Methods**: Finite element modeling, neural networks, dynamic programming optimization (e.g., CPLEX), mixed finite element methods.
- **Industry Experience**: Internships in reservoir engineering, AI R&D for energy applications, field operations in hydraulic fracturing and solar electrification.

## Projects

### 1. Reservoir Numerical Simulation Projects
This folder (`Reservoir_Simulation_Practice`) contains interactive demonstrations of reservoir modeling techniques. Key highlights:
- Jupyter notebooks (.ipynb) for black oil and two-phase flow simulations, including data visualization and parameter sensitivity analysis.
- Datasets for testing implicit/explicit schemes and optimization algorithms.
- Demonstrates efficiency improvements (e.g., 16% via algorithmic optimizations) in simulating pressure distributions and fluid dynamics.

These projects build on my research assistant role at UT Austin, where I implemented C++ algorithms for reservoir simulations using mixed finite element methods.

[View Folder](Reservoir_Simulation_Practice/)

### 2. Two-Phase Proppant Flow Simulation in Hydraulic Fracturing
This project (`Reservoir_Simulation_2Phase`) features a Julia-based numerical simulator for modeling proppant distribution in hydraulic fracturing operations. Using Gridap for finite element discretization, it solves coupled equations for particle and fluid phases, incorporating drag, collision pressure, and virtual mass effects.

- **Key File**: [Simulation_2Phase.jl](Simulation_2Phase.jl) – Core code for the transient solver, weak form residuals, and visualization.
- **Simulation Results**: Includes plots of particle fraction (phi_p), velocities (v_p, v_f), and lambda at final time step. Example output:
  - Particle fraction stabilizes with inlet conditions propagating through the domain.
  - Velocities show phase interactions under gravity and drag forces.
- **VTK Export**: For 3D visualization in tools like Paraview (file: final_solution.vtk).
- **Relevance**: This work extends prior research by incorporating virtual mass effects and repulsive force corrections for stability, directly applicable to hydraulic fracture design optimization.

Run the simulation: `julia Simulation_2Phase.jl` (requires Gridap, Plots packages).

[View Folder](Reservoir_Simulation_2Phase/)


### 3. Field Experience in Hydraulic Fracturing and Team Projects
To complement my simulation work, I've included photos from hands-on field experiences:
- Hydraulic fracturing operations during internships at Hilong Petroleum and Chang’s Group, where I managed engineering operations and coordinated international teams.
- Group work in ASME L.E.A.D. cohort, involving cross-university subsurface projects and site visits.

These experiences honed my skills in applying simulations to practical scenarios, such as optimizing fracture designs and ensuring HSE protocols in remote settings.



## Additional Repositories
- **Wellbore Simulation (Forked)**: Julia-based high-performance computing for wellbore modeling, adapted from PGE383-HPC, where I'm the TA for that course. [Link](https://github.com/Zhouhappyforever/wellbore)

## Contact
- Email: andrewmiller@utexas.edu
- Phone: (512) 203-4555

I'm eager to contribute to teams developing or applying advanced simulators for hydraulic fracturing and reservoir optimization. Feel free to reach out for collaborations or discussions!
