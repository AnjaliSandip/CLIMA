#### Heat model

####
#### 1) Import/Export Needed Functions
####

# Add necessary CliMA functions and sub-routines
using StaticArrays
using CLIMA.VariableTemplates
import CLIMA.DGmethods: BalanceLaw,
                        vars_aux, vars_state, vars_gradient, vars_diffusive,
                        flux_nondiffusive!, flux_diffusive!, source!,
                        gradvariables!, diffusive!, update_aux!, nodal_update_aux!,
                        init_aux!, init_state!,
                        boundary_state!, wavespeed, LocalGeometry


####
#### 2) Define Structs
####

# Introduce needed variables into HeatModel struct
Base.@kwdef struct HeatModel{Fρc, Fα, FiT, Fst} <: BalanceLaw

  # Define heat capacity. This is an input to the model now.
  ρc::Fρc       = (state, aux, t) -> 1   # [ Sand: ρc = 2.49e6 J m-3 K-1 ; Clay: ρc = 2.61e6 J m-3 K-1 ]
  # Replace this with a function that calculates heat capacity (based on liquid+ice)
  # OR Replace this with tabulated values of heat capacity (based on liquid+ice)

  # Define kappa (thermal conductivity). This is an input to the model now.
  α::Fα         = (state, aux, t) -> 10     # [ Sand: λ = 2.42 W m-1 K-1 ; Clay: λ = 1.17 W m-1 K-1 ]

  # Define initial and boundary condition parameters
  initialT::FiT = (aux, t) -> 273.15 + 2.0 # Initial Temperature. This is an input to the model now.
  surfaceT::Fst = (state, aux, t) -> (273.15 + 2.0) # Surface boundary condition. This is an input to the model now.
end


# --------------------------------- 3) Define CliMA vars ---------------------------------------

# Stored in the aux state are:
#   `coord` coordinate points (needed for BCs)
#   `u` advection velocity
#   `D` Diffusion tensor
vars_aux(::HeatModel, FT) = @vars(z::FT, T::FT) # stored dg.auxstate
vars_state(::HeatModel, FT) = @vars(ρcT::FT, q_tot::FT) # stored in Q
vars_gradient(::HeatModel, FT) = @vars(T::FT) # not stored
vars_diffusive(::HeatModel, FT) = @vars(∇T::SVector{3,FT}) # stored in dg.diffstate

# --------------------------------- 4) CliMA functions needed for simulation -------------------

# ---------------- 4a) Update states

# Update all auxiliary variables
function update_aux!(
    dg::DGModel,
    m::HeatModel,
    Q::MPIStateArray,
    t::Real,
    elems::UnitRange,
)
  nodal_update_aux!(soil_nodal_update_aux!, dg, m, Q, t, elems)
  return true
end
# Update all auxiliary nodes
function soil_nodal_update_aux!(
  m::HeatModel,
  state::Vars,
  aux::Vars,
  t::Real)
  aux.T = state.ρcT / m.ρc(state, aux, t)
end

# ---------------- 4b) Calculate state and derivative of T

# Calculate T based on internal energy state variable
function gradvariables!(
    m::HeatModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
  transform.T = state.ρcT / m.ρc(state, aux, t)
end
# Gradient of T calculation
function diffusive!(
    m::HeatModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
  diffusive.∇T = ∇transform.T
end
# Calculate thermal flux (non-diffusive (?))
function flux_nondiffusive!(
    m::HeatModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
  )
end
# Calculate thermal flux (diffusive (?))
function flux_diffusive!(
    m::HeatModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
  )
   flux.ρcT -= m.α(state, aux, t) * diffusive.∇T
end

# ---------------- 4c) Extra Sources

# Introduce sources of energy (e.g. Metabolic heat from microbes)
function source!(m::HeatModel, state::Vars, _...)
end

# ---------------- 4d) Initialization

# Initialize z-Profile
function init_aux!(m::HeatModel, aux::Vars, geom::LocalGeometry)
  aux.z = geom.coord[3]
  aux.T = m.initialT(aux, 0)
end

# Initialize State variables from T to internal energy
function init_state!(m::HeatModel, state::Vars, aux::Vars, coords, t::Real)
  state.ρcT = m.ρc(state, aux, t) * aux.T
end

# ---------------- 4e) Boundary Conditions

# Boundary condition function
function boundary_state!(nf, m::HeatModel, state⁺::Vars, aux⁺::Vars,
                         nM, state⁻::Vars, aux⁻::Vars, bctype, t, _...)
  if bctype == 1
    # surface
    state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    nothing
  end
end
# Boundary condition function - repeated?
function boundary_state!(nf, m::HeatModel, state⁺::Vars, diff⁺::Vars,
                         aux⁺::Vars, nM, state⁻::Vars, diff⁻::Vars, aux⁻::Vars,
                         bctype, t, _...)
  if bctype == 1
    # surface
    state⁺.ρcT = m.ρc(state⁻, aux⁻, t) * m.surfaceT(state⁻, aux⁻, t)
  elseif bctype == 2
    # bottom
    diff⁺.∇T = -diff⁻.∇T
  end
end
