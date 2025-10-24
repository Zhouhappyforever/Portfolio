# --------------------------------------------------------------------------
# SECTION 1: PREAMBLE AND PACKAGE IMPORTS
# --------------------------------------------------------------------------
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.MultiField
using Gridap.ODEs
using Gridap.CellData
using Plots

# Debugging: Check package import
println("Section 1: Packages imported successfully.")

# --------------------------------------------------------------------------
# SECTION 2: SIMULATION PARAMETERS
# --------------------------------------------------------------------------
# Domain parameters
L = 0.1            # Length of the domain (m)
n = 100            # Number of elements
order = 2          # Polynomial order for finite elements

# Material and fluid properties
rho_p = 2650.0     # Particle density (kg/m³)
rho_f = 1000.0     # Fluid density (kg/m³)
mu_f = 0.001       # Fluid viscosity (Pa·s)
g = 9.81           # Gravity (m/s²)
theta = 0.0        # Inclination angle (radians, 0 for horizontal)
d_part = 100e-6    # Particle diameter (m)
a = 3.0            # Coefficient for collision pressure
C_vm = 0.5         # Virtual mass coefficient
phi_max = 0.63     # Maximum packing fraction

# Simplification constants
drag_const = 1.0   # Constant drag for testing
bpress_const = 0.0 # Zero bpress for testing

# Time stepping parameters
T = 0.01           # Final time (s)
num_steps = 10     # Reduced for testing
dt = T / num_steps # Time step size (s)

# Initial and boundary conditions
phi_p0 = 0.0       # Initial particle fraction
v_p0 = 0.0         # Initial particle velocity (m/s)
v_f0 = 0.1         # Initial fluid velocity (m/s)
lam0 = 0.0         # Initial lambda

phi_inlet = 0.5    # Inlet particle fraction
vp_inlet = 0.1     # Inlet particle velocity (m/s)
vf_inlet = 0.1     # Inlet fluid velocity (m/s)
lam_outlet = 0.0   # Outlet lambda

# 1D direction vector for scalar derivatives
dir = VectorValue(1.0)

# Debugging: Print parameters
println("Section 2: Parameters defined successfully.")
println("Domain length L = $L m, elements n = $n, order = $order")

# --------------------------------------------------------------------------
# SECTION 3: AUXILIARY FUNCTIONS
# --------------------------------------------------------------------------
# Repulsive force correction for scalars
function f_rep(p::Number)
  clamped_p = max(p, 0.0)
  clamped_p = min(clamped_p, 0.99)
  rep_term = (clamped_p / max(1 - clamped_p, 1e-10))^2.05
  1 + (rep_term - 1) * (1 / (1 + exp(-100 * (clamped_p - 0.1))))
end

# Repulsive force correction for fields
function f_rep(phi::CellField)
  Operation(f_rep)(phi)
end

# Drag coefficient (variable version; use drag_const for testing)
function dp(phi)
  18 * mu_f / (d_part^2) * phi * f_rep(phi)
end

# Collision dispersive pressure (variable version; use bpress_const for testing)
function bp(phi, dvp_ds)
  rho_p * (a^2) * (phi^3) * (dvp_ds^2)
end

# Sin theta (constant for horizontal pipe)
sin_theta = sin(theta)

# Debugging: Test auxiliary functions
test_phi = 0.2
test_dvp_ds = 1.0
println("Section 3: Auxiliary functions defined.")
println("Test dp(0.2) = $(dp(0.2))")
println("Test dp(0.05) = $(dp(0.05))")
println("Test dp(-0.01) = $(dp(-0.01))")
println("Test dp(0.99) = $(dp(0.99))")
println("Test bp(0.2, 1.0) = $(bp(0.2, 1.0))")

# --------------------------------------------------------------------------
# SECTION 4: MESH AND FUNCTION SPACES
# --------------------------------------------------------------------------
# Create 1D model
model = CartesianDiscreteModel((0.0, L), (n,))

# Labels for boundaries
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "left", [1])
add_tag_from_tags!(labels, "right", [2])

# Quadrature degree
degree = 2 * order

# Triangulation and measure
Omega = Triangulation(model)
dOmega = Measure(Omega, degree)

# Reference FE
reffe = ReferenceFE(lagrangian, Float64, order)

# Test spaces
V_phi = TestFESpace(model, reffe, dirichlet_tags=["left"])
V_vp = TestFESpace(model, reffe, dirichlet_tags=["left"])
V_vf = TestFESpace(model, reffe, dirichlet_tags=["left"])
V_lam = TestFESpace(model, reffe, dirichlet_tags=["right"])

# Trial spaces with Dirichlet conditions
Phi = TrialFESpace(V_phi, x -> phi_inlet)
Vp = TrialFESpace(V_vp, x -> vp_inlet)
Vf = TrialFESpace(V_vf, x -> vf_inlet)
Lam = TrialFESpace(V_lam, x -> lam_outlet)

# Multi-field test and trial spaces
Y = MultiFieldFESpace([V_phi, V_vp, V_vf, V_lam])
X = MultiFieldFESpace([Phi, Vp, Vf, Lam])

# Transient multi-field space
X_t = TransientMultiFieldFESpace([Phi, Vp, Vf, Lam])

# Debugging: Check spaces
println("Section 4: Mesh and function spaces created successfully.")
println("Number of dofs: phi=$(num_free_dofs(V_phi)), vp=$(num_free_dofs(V_vp)), vf=$(num_free_dofs(V_vf)), lam=$(num_free_dofs(V_lam))")

# --------------------------------------------------------------------------
# SECTION 5: INITIAL CONDITIONS
# --------------------------------------------------------------------------
# Initial fields
phi_init(x) = phi_p0
vp_init(x) = v_p0
vf_init(x) = v_f0
lam_init(x) = lam0

# Interpolate initial conditions
xh0 = interpolate_everywhere([phi_init, vp_init, vf_init, lam_init], X_t(0.0))

# Debugging: Check initial conditions
println("Section 5: Initial conditions set successfully.")

# --------------------------------------------------------------------------
# SECTION 6: WEAK FORM RESIDUAL
# --------------------------------------------------------------------------
function res(t, xt, yt)
    phi_t, vp_t, vf_t, lam_t = xt
    psi, wp, wf, chi = yt
    dphi_t = ∂t(phi_t)
    dvp_t = ∂t(vp_t)
    dvf_t = ∂t(vf_t)
    phi_f = 1.0 - phi_t
    drag = drag_const   # Constant for testing
    dvp_ds = ∇(vp_t) ⋅ dir
    bpress = bpress_const  # Zero for testing
    rel_vel = vf_t - vp_t

    # Scalar derivatives
    dphi_ds = ∇(phi_t) ⋅ dir
    dpsi_ds = ∇(psi) ⋅ dir
    dphi_f_ds = ∇(phi_f) ⋅ dir
    dchi_ds = ∇(chi) ⋅ dir
    dvf_ds = ∇(vf_t) ⋅ dir
    dwp_ds = ∇(wp) ⋅ dir
    dwf_ds = ∇(wf) ⋅ dir

    # Particle mass balance (skew-symmetric form for stability)
    particle_mass = ∫(psi * dphi_t - 0.5 * (psi * dphi_ds * vp_t - phi_t * dpsi_ds * vp_t)) * dOmega

    # Fluid mass balance (skew-symmetric)
    fluid_mass = ∫(chi * (-dphi_t) - 0.5 * (chi * dphi_f_ds * vf_t - phi_f * dchi_ds * vf_t)) * dOmega

    # Particle momentum balance
    particle_mom = ∫(wp * rho_p * (dvp_t + vp_t * dvp_ds)) * dOmega -
                   ∫(wp * rho_p * g * sin_theta) * dOmega -
                   ∫(wp * drag * rel_vel) * dOmega -
                   ∫(wp * C_vm * rho_f * phi_t * (dvf_t + vf_t * dvf_ds - dvp_t - vp_t * dvp_ds)) * dOmega +
                   ∫(dwp_ds * (phi_t * lam_t)) * dOmega +
                   ∫(dwp_ds * bpress) * dOmega

    # Fluid momentum balance
    fluid_mom = ∫(wf * rho_f * (dvf_t + vf_t * dvf_ds)) * dOmega -
                ∫(wf * rho_f * g * sin_theta) * dOmega +
                ∫(wf * drag * rel_vel) * dOmega +
                ∫(wf * C_vm * rho_f * phi_t * (dvf_t + vf_t * dvf_ds - dvp_t - vp_t * dvp_ds)) * dOmega +
                ∫(dwf_ds * (phi_f * lam_t)) * dOmega

    # Debug: Print types of selected terms (only for first time step to reduce output)
    if t == 0.0
        println("DEBUG: Type of particle_mass: $(typeof(particle_mass))")
        println("DEBUG: Type of fluid_mass: $(typeof(fluid_mass))")
        println("DEBUG: Type of particle_mom: $(typeof(particle_mom))")
        println("DEBUG: Type of fluid_mom: $(typeof(fluid_mom))")
    end

    return particle_mass + fluid_mass + particle_mom + fluid_mom
end

# Debugging: Weak form defined
println("Section 6: Weak form residual defined successfully.")

# --------------------------------------------------------------------------
# SECTION 7: TRANSIENT OPERATOR AND SOLVER
# --------------------------------------------------------------------------
# Transient FE operator
op = TransientFEOperator(res, X_t, Y)

# Nonlinear solver
ls = LUSolver()
nls = NLSolver(ls; method=:newton, iterations=10, show_trace=true)

# ODE solver (Backward Euler)
theta_val = 1.0
ode_solver = ThetaMethod(nls, dt, theta_val)

# Time range
t0 = 0.0
tF = T

# Debugging: Solver setup
println("Section 7: Transient operator and solver configured successfully.")

# --------------------------------------------------------------------------
# SECTION 8: SOLVE THE TRANSIENT PROBLEM
# --------------------------------------------------------------------------
sol_t = Gridap.solve(ode_solver, op, xh0, t0, tF)

# Collect solutions for plotting
all_solutions = []
let iter = 0
    for sol_item in sol_t
        if iter < 3  # Limit to first 3 for brevity
            println("DEBUG: Type of sol_item: $(typeof(sol_item))")
        end
        xh, t = sol_item  # Direct unpack: first MultiFieldFEFunction, second Float64
        if iter < 3  # Limit to first 3 for brevity
            println("DEBUG: Type of xh: $(typeof(xh))")
            println("DEBUG: Type of t: $(typeof(t))")
        end
        iter += 1
        push!(all_solutions, (xh, t))
        if iter <= 3 || iter % 10 == 0  # Print for first 3 and every 10th
            println("DEBUG: Solved time step $iter, t = $t")
        end
    end
end

# Debugging: Solution completed
println("Section 8: Transient problem solved successfully. Total steps: $(length(all_solutions))")

# --------------------------------------------------------------------------
# SECTION 9: IN-LINE PLOTTING OF FINAL RESULT
# --------------------------------------------------------------------------
println("DEBUG: Checking if final solution was captured...")
if !isempty(all_solutions)
    println("Generating final plots...")
    
    # Extract the final solution
    final_xh, final_t = all_solutions[end]
    println("DEBUG: Type of final_xh: $(typeof(final_xh))")
    println("DEBUG: Type of final_t: $(typeof(final_t))")
    
    # Unpack fields
    final_phi, final_vp, final_vf, final_lam = final_xh
    
    # Get nodal coordinates and evaluate fields
    coords = Gridap.Geometry.get_node_coordinates(model)
    println("DEBUG: Type of coords: $(typeof(coords))")
    println("DEBUG: Length of coords: $(length(coords))")
    x = [p[1] for p in coords]
    
    # Evaluate fields at coordinates
    phi_vals = final_phi.(coords)
    vp_vals = final_vp.(coords)
    vf_vals = final_vf.(coords)
    lam_vals = final_lam.(coords)
    
    # Debug: Print sample values
    println("DEBUG: Sample phi_vals: $(phi_vals[1:5])")
    println("DEBUG: Sample vp_vals: $(vp_vals[1:5])")
    println("DEBUG: Sample vf_vals: $(vf_vals[1:5])")
    println("DEBUG: Sample lam_vals: $(lam_vals[1:5])")
    
    p1 = plot(x, phi_vals, label="phi_p", title="Particle Fraction at t = $(round(final_t, digits=5)) s", xlabel="s (m)", ylabel="phi_p")
    p2 = plot(x, vp_vals, label="v_p", title="Velocities at t = $(round(final_t, digits=5)) s", xlabel="s (m)", ylabel="v (m/s)")
    plot!(p2, x, vf_vals, label="v_f")
    p3 = plot(x, lam_vals, label="lam", title="Lambda at t = $(round(final_t, digits=5)) s", xlabel="s (m)", ylabel="lam")

    # Display side-by-side
    display(plot(p1, p2, p3, layout=(1,3), size=(1200, 400)))
    
    # Save plot for inspection
    savefig("final_plot.png")
    println("Plot saved as final_plot.png")
else
    println("Final solution not found for plotting because no time steps were collected.")
end

# Debugging: Confirm plots
println("Section 9: In-line plotting completed successfully.")