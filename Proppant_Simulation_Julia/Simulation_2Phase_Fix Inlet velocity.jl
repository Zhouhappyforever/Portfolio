# --------------------------------------------------------------------------
# FINAL ROBUST PROPPANT SIMULATION (Fixed BCs & Residual)
# --------------------------------------------------------------------------
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.MultiField
using Gridap.ODEs
using Gridap.CellData
using Gridap.Visualization: createpvd, createvtk
using LinearAlgebra

println("Starting Simulation...")

# --------------------------------------------------------------------------
# 1. PARAMETERS
# --------------------------------------------------------------------------
# Increased element count to 20 for better resolution
L, n, order = 20.0, 20, 1  
rho_p, rho_f, mu_f = 2650.0, 1000.0, 0.002
d_part, phi_max = 0.001, 0.63
D_pipe, f_fric = 0.1, 0.02

# Regularization:
# visc_p and visc_f damp velocity oscillations,
# visc_lam strongly smooths the pressure field so that
# it does not dominate the dynamics.
visc_p, visc_f, visc_mass, visc_lam = 1.0e-2, 1.0e-2, 0.1, 1.0

# Time Control
T, dt = 20.0, 0.01

# Boundary Values
phi_inlet, v_inlet = 0.1, 1.5
q_in = phi_inlet * v_inlet   # inlet proppant flux φ_in * v_in
lam_inlet, lam_outlet = 5.0e4, 0.02
beta = 1.0
# Relaxation rate driving phi toward phi_inlet in the inlet cell
k_inlet = v_inlet / (L/n)  # characteristic rate ~ velocity / cell size
dir_1d = VectorValue(1.0) # Explicit Direction

# --------------------------------------------------------------------------
# 2. MESH & MANUAL TAGGING (The Robust Method)
# --------------------------------------------------------------------------
model = CartesianDiscreteModel((0.0, L), (n,))

# Get the Labeling Object
labels = get_face_labeling(model)

# 1. Access the Internal Tag Array for Vertices (Dimension 0)
#    This bypasses 'add_tag!' issues by writing directly to the memory.
vertex_tags = get_face_tag(labels, 0)

# 2. Force Node 1 (Inlet) and Node N+1 (Outlet) to have unique IDs
const INLET_ID  = 100
const OUTLET_ID = 101

# In 1D Cartesian, Node 1 is always Left, Node 'end' is always Right.
vertex_tags[1]   = INLET_ID
vertex_tags[end] = OUTLET_ID

# 3. Register these IDs as String Tags
add_tag!(labels, "inlet_fix",  [INLET_ID])
add_tag!(labels, "outlet_fix", [OUTLET_ID])

println("Boundary Tags Applied: Node 1 -> 'inlet_fix', Node $(length(vertex_tags)) -> 'outlet_fix'")

# --------------------------------------------------------------------------
# 3. FE SPACES
# --------------------------------------------------------------------------
reffe = ReferenceFE(lagrangian, Float64, 1)

# Inlet: Fixed Velocity & Phi
V_phi = TestFESpace(model, reffe, dirichlet_tags="inlet_fix")
V_vp  = TestFESpace(model, reffe, dirichlet_tags="inlet_fix")
V_vf  = TestFESpace(model, reffe, dirichlet_tags="inlet_fix")

# Outlet: Fixed Pressure (reference level). Inlet pressure is free.
V_lam = TestFESpace(model, reffe, dirichlet_tags="outlet_fix")

# Boundary Functions
phi_bc(x) = phi_inlet
v_bc(x)   = v_inlet
lam_bc(x) = lam_outlet

# Trial Spaces
Phi = TransientTrialFESpace(TrialFESpace(V_phi, phi_bc))
Vp  = TransientTrialFESpace(TrialFESpace(V_vp, v_bc))
Vf  = TransientTrialFESpace(TrialFESpace(V_vf, v_bc))
Lam = TransientTrialFESpace(TrialFESpace(V_lam, lam_bc))

Y   = MultiFieldFESpace([V_phi, V_vp, V_vf, V_lam])
X_t = TransientMultiFieldFESpace([Phi, Vp, Vf, Lam])

println("Spaces Configured. DOFs should be locked.")

# --------------------------------------------------------------------------
# 4. INITIAL CONDITIONS
# --------------------------------------------------------------------------
degree = 2 * order
Omega  = Triangulation(model)
dOmega = Measure(Omega, degree)

# Boundary triangulation for flux terms (inlet only)
Γ_in   = BoundaryTriangulation(model, tags="inlet_fix")
dΓ_in  = Measure(Γ_in, degree)

# In 1D, we can use a scalar inlet normal (sign can be adjusted if needed)
const n_in = 1.0

# We define h based on mesh size for the plug
h_mesh = L/n
phi_bg = 1.0e-6

phi_init(x) = x[1] <= h_mesh ? phi_inlet : phi_bg
vp_init(x)  = v_inlet
vf_init(x)  = v_inlet
lam_init(x) = lam_outlet

xh0 = interpolate_everywhere([phi_init, vp_init, vf_init, lam_init], X_t(0.0))

# Indicator field for the inlet region (first cell): chi_in ≈ 1 near inlet, 0 elsewhere
chi_in = interpolate_everywhere(x -> (x[1] <= h_mesh ? 1.0 : 0.0), V_phi)

# --------------------------------------------------------------------------
# 5. RESIDUAL (Robust Version)
# --------------------------------------------------------------------------
clamp_nonneg(p) = max(p, 1.0e-6)
f_rep(phi) = 1.0 
f_rep(phi::CellField) = Operation(f_rep)(phi)
function dp(phi)
  # Drag coefficient per unit volume: β(φ) ∝ φ, consistent with
  # φ ρ_p ∂t v_p = β(φ) (v_f - v_p). The effective coupling β(φ)/φ
  # is then O(1) and independent of φ.
  c_phi = isa(phi, Number) ? clamp_nonneg(phi) : Operation(clamp_nonneg)(phi)
  18 * mu_f / (d_part^2) * c_phi * f_rep(c_phi)
end
function fric_term(vf)
  abs_vf = sqrt(vf * vf + 1e-4)
  f_fric * rho_f * abs_vf * vf / (2 * D_pipe)
end

function res(t, xt, yt, dOmega)
  phi_t, vp_t, vf_t, lam_t = xt
  psi,  wp,   wf,   q      = yt

  phi_phys = Operation(clamp_nonneg)(phi_t)

  # Derivatives
  dphi_t, dvp_t, dvf_t, dlam_t = ∂t(phi_t), ∂t(vp_t), ∂t(vf_t), ∂t(lam_t)
  dvp_ds = ∇(vp_t)⋅dir_1d
  dvf_ds = ∇(vf_t)⋅dir_1d
  dlam_ds = ∇(lam_t)⋅dir_1d
  dphi_ds = ∇(phi_t)⋅dir_1d

  # Physics
  drag    = dp(phi_phys)
  rel_vel = vf_t - vp_t
  mix_flux = phi_t*vp_t + (1.0-phi_t)*vf_t

  # Residual term for mass conservation in conservative form with an inlet source:
  #   ∂t φ + ∂s(φ vp) - visc_mass ∂ss φ = S_in
  # S_in approximates the effect of specifying φ v_p = q_in at x = 0 over the first cell.
  S_in = k_inlet * chi_in * (phi_inlet - phi_phys)
  res_mass = ∫( psi*dphi_t
                - (∇(psi)⋅dir_1d)*(phi_t*vp_t)
                + visc_mass*(∇(psi)⋅dir_1d)*dphi_ds
                - psi*S_in )*dOmega

  res_vp_term = phi_phys*rho_p*dvp_t - drag*rel_vel
  res_vf_term = rho_f*dvf_t + drag*rel_vel + dlam_ds + Operation(fric_term)(vf_t)
  
  res_mom_p = ∫( wp*res_vp_term + visc_p*(∇(wp)⋅dir_1d)*dvp_ds )*dOmega
  res_mom_f = ∫( wf*res_vf_term + visc_f*(∇(wf)⋅dir_1d)*dvf_ds )*dOmega

  # Continuity / pressure equation in quasi-steady form:
  #   ∂s (mix_flux) = 0, regularized by a strong pressure Laplacian.
  res_cont  = ∫( -(∇(q)⋅dir_1d)*mix_flux
                 + visc_lam*(∇(q)⋅dir_1d)*dlam_ds )*dOmega

  return res_mass + res_mom_p + res_mom_f + res_cont
end

# --------------------------------------------------------------------------
# 6. SOLVER
# --------------------------------------------------------------------------
op = TransientFEOperator((t,x,y) -> res(t,x,y,dOmega), X_t, Y)

nls = NLSolver(LUSolver(); method=:newton, iterations=50, show_trace=true)
ode_solver = ThetaMethod(nls, dt, 1.0) 

vtk_dir = joinpath(homedir(), "Downloads/vtk")
isdir(vtk_dir) || mkdir(vtk_dir)

println("Solving...")
createpvd(joinpath(vtk_dir, "proppant_timeseries")) do pvd
  iter = 0
  for (xh, t) in Gridap.solve(ode_solver, op, xh0, 0.0, T)
    iter += 1
    if iter % 10 == 0 || iter == 1
      pvd[t] = createvtk(Omega, joinpath(vtk_dir, "step_$(lpad(iter,4,'0'))");
                         cellfields = ["phi"=>xh[1], "vp"=>xh[2], "vf"=>xh[3], "p"=>xh[4]])
      println("Step $iter: t=$t")
    end
  end
end
println("Simulation Complete.")