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
using Gridap.Visualization: createpvd, createvtk
using LineSearches: BackTracking
using LinearAlgebra

println("Section 1: Packages imported successfully.")

# --------------------------------------------------------------------------
# SECTION 2: SIMULATION PARAMETERS
# --------------------------------------------------------------------------
L = 20           # Length of the domain (m)
n = 10           # Number of elements
order = 1        # Linear for simplicity

# Material and fluid properties
rho_p = 2650.0   # Particle density (kg/m³)
rho_f = 1000.0   # Fluid density (kg/m³)
mu_f  = 0.002    # Fluid viscosity (Pa·s)
g = 9.81         # Gravity (m/s²)
theta = 0.0      # Horizontal pipe
d_part = 0.001   # Particle diameter (m)
a = 0.0          # No collision pressure
C_vm = 0.0       # Virtual mass disabled
phi_max = 0.63

visc_p     = 0.0
visc_f     = mu_f
visc_f_art = 0.0
visc_mass  = 0.0

D_pipe = 0.1
f_fric = 0.02
h = L / n
v_target = 1.5          # target / inlet velocity (used as a scale, not BC)
v_inlet  = v_target

# artificial compressibility not used (no lambda unknown)
beta = 0.0

# Time stepping
T  = 20.0                 # Final time
dt = 0.01                # Time step

# BC parameters
phi_p0   = 1e-6          # Initial proppant fraction in the pipe
v_p0     = 0.0
v_f0     = 0.0

phi_inlet   = 0.1        # Inlet proppant fraction

# Interpret pressure as not solved but used via constant gradient.
# lam_inlet and lam_outlet define the pressure at boundaries for reference.
lam_inlet   = 5.0e4      # Inlet pressure (Pa)
lam_outlet  = 0.0        # Outlet pressure (Pa)

p_grad = (lam_inlet - lam_outlet) / L   # ≈ -∂p/∂s

println("Section 2: Parameters defined successfully.")
println("Domain length L = $L m, elements n = $n, order = $order")
println("Time stepping: dt = $dt, T = $T")
println("Driving pressure gradient p_grad = $p_grad Pa/m")
println("Inlet: phi = $phi_inlet, pin = $lam_inlet; pout = $lam_outlet")

println("Velocity scale (not fixed BC): v_target = $v_target m/s")

# Direction vector in 1D for computing directional derivatives (∂/∂s)
const dir = VectorValue(1.0)

# --------------------------------------------------------------------------
# SECTION 3: AUXILIARY FUNCTIONS
# --------------------------------------------------------------------------
function f_rep(p::Number)
  1.0
end

function f_rep(phi::CellField)
  Operation(f_rep)(phi)
end

# keep φ ≥ 0 so drag never vanishes, avoiding singular Jacobians
clamp_nonneg(p) = max(p, 0.0)

function dp(phi)
  # Drag coefficient β(φ) ~ 18 μ_f / d^2 * f_rep(φ_clamped)
  clamped_phi = isa(phi, Number) ? clamp_nonneg(phi) : Operation(clamp_nonneg)(phi)
  18 * mu_f / (d_part^2) * f_rep(clamped_phi)
end

bp(phi, dvp_ds) = 0.0

sin_theta = sin(theta)

function fric_term(vf)
  abs_vf = sqrt(vf * vf + 1e-4)
  f_fric * rho_f * abs_vf * vf / (2 * D_pipe)
end

# Positive-part cutoff for velocities (used in inflow BC if needed)
cutoff_pos(x) = max(x, 0.0)

println("Section 3: Auxiliary functions defined.")
println("Test dp(0.2)  = $(dp(0.2))")
println("Test dp(0.05) = $(dp(0.05))")
println("Test dp(-0.01)= $(dp(-0.01))")
println("Test dp(0.99) = $(dp(0.99))")
println("Test bp(0.2, 1.0) = $(bp(0.2, 1.0))")

# --------------------------------------------------------------------------
# SECTION 4: MESH AND FUNCTION SPACES
# --------------------------------------------------------------------------
model = CartesianDiscreteModel((0.0, L), (n,))

labels = get_face_labeling(model)
add_tag_from_tags!(labels, "left",  [1])
add_tag_from_tags!(labels, "right", [2])

degree = 2 * order
Omega  = Triangulation(model)
dOmega = Measure(Omega, degree)

Gamma_left  = BoundaryTriangulation(model, tags = "left")
Gamma_right = BoundaryTriangulation(model, tags = "right")

dGamma_left  = Measure(Gamma_left, degree)
dGamma_right = Measure(Gamma_right, degree)

reffe_1 = ReferenceFE(lagrangian, Float64, 1)

# BCs:
#  phi : no essential BC; inlet concentration imposed weakly via flux on Γ_left
#  vp  : free (driven by pressure gradient + drag)
#  vf  : free (driven by pressure gradient + drag + friction)
V_phi = TestFESpace(model, reffe_1)
V_vp  = TestFESpace(model, reffe_1)
V_vf  = TestFESpace(model, reffe_1)

# Trial space for phi (no essential BC; inflow concentration handled in weak form)
Phi_trial = TrialFESpace(V_phi)

# vp and vf have no essential boundary conditions
Vp_trial  = TrialFESpace(V_vp)
Vf_trial  = TrialFESpace(V_vf)

Phi = TransientTrialFESpace(Phi_trial)
Vp  = TransientTrialFESpace(Vp_trial)
Vf  = TransientTrialFESpace(Vf_trial)

Y   = MultiFieldFESpace([V_phi, V_vp, V_vf])
X_t = TransientMultiFieldFESpace([Phi, Vp, Vf])

println("Section 4: Mesh and function spaces created successfully.")
println("Number of dofs: phi=$(num_free_dofs(V_phi)), vp=$(num_free_dofs(V_vp)), vf=$(num_free_dofs(V_vf))")

# --------------------------------------------------------------------------
# SECTION 5: INITIAL CONDITIONS
# --------------------------------------------------------------------------
phi_init(x) = phi_p0
vp_init(x)  = v_p0
vf_init(x)  = v_f0

t0  = 0.0
xh0 = interpolate_everywhere([phi_init, vp_init, vf_init], X_t(t0))

println("Section 5: Initial conditions set successfully.")

# --------------------------------------------------------------------------
# SECTION 6: WEAK FORM RESIDUAL
# --------------------------------------------------------------------------
function res(t, xt, yt)
  # Unknowns and tests
  phi_t, vp_t, vf_t = xt
  psi,   wp,   wf   = yt

  # Time derivatives
  dphi_t = ∂t(phi_t)
  dvp_t  = ∂t(vp_t)
  dvf_t  = ∂t(vf_t)

  # Volume fraction of fluid
  phi_f = 1.0 - phi_t

  # Drag and relative velocity
  drag    = dp(phi_t)
  rel_vel = vf_t - vp_t

  # Spatial derivatives in 1D
  dphi_ds = ∇(phi_t) ⋅ dir
  dpsi_ds = ∇(psi)   ⋅ dir

  # Proppant mass balance:
  #   ∂t φ + ∂s(φ v_f) = 0  in Ω,
  # with an imposed inflow concentration φ_inlet at Γ_left through the
  # flux term. The weak form (1D) is
  #   ∫ ψ ∂t φ dΩ - ∫ (∂ψ/∂s) φ v_f dΩ
  #   + ∫_{Γ_right} ψ (φ v_f) dΓ - ∫_{Γ_left} ψ (φ_inlet v_f) dΓ = 0.
  #
  # Note: φ has no essential (Dirichlet) boundary condition now; the inlet
  # concentration is enforced purely via the boundary flux on Γ_left.

  particle_mass =
      ∫( psi * dphi_t - dpsi_ds * (phi_t * vf_t) ) * dOmega +
      ∫( psi * (phi_t * vf_t) ) * dGamma_right -
      ∫( psi * (phi_inlet * vf_t) ) * dGamma_left

  # Particle momentum
  particle_mom =
      ∫( wp * (rho_p * dvp_t - p_grad - drag * rel_vel) ) * dOmega

  # Fluid momentum
  fluid_mom =
      ∫( wf * (rho_f * dvf_t - p_grad + drag * rel_vel + Operation(fric_term)(vf_t)) ) * dOmega

  if t == 0.0
    println("DEBUG: Residual types at t=0:")
    println("  particle_mass :: ", typeof(particle_mass))
    println("  particle_mom  :: ", typeof(particle_mom))
    println("  fluid_mom     :: ", typeof(fluid_mom))
  end

  return particle_mass + particle_mom + fluid_mom
end

println("Section 6: Weak form residual defined successfully.")

# --------------------------------------------------------------------------
# SECTION 7: TRANSIENT OPERATOR AND SOLVER
# --------------------------------------------------------------------------
op = TransientFEOperator(res, X_t, Y)

ls  = LUSolver()
nls = NLSolver(ls; method=:newton,
                   linesearch=BackTracking(),
                   iterations=50,
                   ftol=1e-7,
                   xtol=1e-8,
                   show_trace=true)

theta_val  = 1.0
ode_solver = ThetaMethod(nls, dt, theta_val)

tF = T

println("Section 7: Transient operator and solver configured successfully.")

# --------------------------------------------------------------------------
# SECTION 8: SOLVE THE TRANSIENT PROBLEM AND EXPORT VTK
# --------------------------------------------------------------------------
vtk_dir = joinpath(homedir(), "Downloads/vtk")
isdir(vtk_dir) || mkdir(vtk_dir)

iter     = Ref(0)
final_xh = xh0
final_t  = t0

sol_t = Gridap.solve(ode_solver, op, xh0, t0, tF)

createpvd(joinpath(vtk_dir, "proppant_timeseries")) do pvd
  for (xh, t) in sol_t
    iter[] += 1
    final_xh = xh
    final_t  = t

    if iter[] <= 3 || iter[] % 10 == 0
      println("DEBUG: Solved time step $(iter[]), t = $t")
    end

    if iter[] % 100 == 0 || iter[] == 1
      phi_k, vp_k, vf_k = xh

      # Optional: analytic pressure field for VTK
      p_fun  = x -> lam_inlet - p_grad * x[1]
      p_cell = CellField(p_fun, Omega)

      vtk = createvtk(Omega, joinpath(vtk_dir, "proppant_step_$(lpad(iter[], 5, '0'))");
        cellfields = [
          "phi" => phi_k,
          "vp"  => vp_k,
          "vf"  => vf_k,
          "p"   => p_cell
        ]
      )
      pvd[t] = vtk
    end
  end
  println("Wrote PVD time series: proppant_timeseries.pvd")
  println("Section 8: Transient problem solved successfully. Total steps: $(iter[])")
end

# --------------------------------------------------------------------------
# SECTION 9: IN-LINE PLOTTING OF FINAL RESULT
# --------------------------------------------------------------------------
println("DEBUG: Checking if final solution was captured...")
println("Generating final plots...")

final_phi, final_vp, final_vf = final_xh

coords = Gridap.Geometry.get_node_coordinates(model)
x = [p[1] for p in coords]

phi_vals = final_phi.(coords)
vp_vals  = final_vp.(coords)
vf_vals  = final_vf.(coords)

# Analytic pressure for plotting (optional)
p_vals = [lam_inlet - p_grad * xi for xi in x]

println("DEBUG: Sample phi_vals: $(phi_vals[1:5])")
println("DEBUG: Sample vp_vals:  $(vp_vals[1:5])")
println("DEBUG: Sample vf_vals:  $(vf_vals[1:5])")

println("DEBUG: phi range = ", (minimum(phi_vals), maximum(phi_vals)))
println("DEBUG: vp  range = ", (minimum(vp_vals),  maximum(vp_vals)))
println("DEBUG: vf  range = ", (minimum(vf_vals),  maximum(vf_vals)))

p1 = plot(x, phi_vals, label="phi_p",
          title="Particle Fraction at t = $(round(final_t, digits=5)) s",
          xlabel="s (m)", ylabel="phi_p")

p2 = plot(x, vp_vals, label="v_p",
          title="Velocities at t = $(round(final_t, digits=5)) s",
          xlabel="s (m)", ylabel="v (m/s)")
plot!(p2, x, vf_vals, label="v_f")

display(plot(p1, p2, layout=(1,2), size=(900, 400)))
savefig("final_plot.png")
println("Plot saved as final_plot.png")

println("Section 9: In-line plotting completed successfully.")