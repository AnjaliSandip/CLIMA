module Courant


export Advective_CFL, Diffusive_CFL

using LinearAlgebra, StaticArrays
using ..VariableTemplates
using ..MoistThermodynamics
using ..PlanetParameters
import ..MoistThermodynamics: internal_energy
using ..SubgridScaleParameters
using GPUifyLoops
using ..MPIStateArrays: MPIStateArray
using CLIMA.Atmos: AtmosModel, soundspeed
import CLIMA.DGmethods: BalanceLaw, vars_aux, vars_state,
                        vars_diffusive
using ..Mesh.Topologies
using ..Mesh.Grids
using ..Mesh.Grids: VerticalDirection, HorizontalDirection, EveryDirection
using Logging
using Printf
using LinearAlgebra





function Advective_CFL(m::AtmosModel, state::Vars, aux::Vars,
                               diffusive::Vars, Δx, Δt, direction=VerticalDirection())
  if direction isa VerticalDirection
    w = state.ρu[3]/state.ρ
          return Δt * (abs(w) + soundspeed(m.moisture, m.orientation, state,
                                            aux)) / Δx
  elseif direction isa HorizontalDirection
    u = max(state.ρu[1]/state.ρ,state.ρu[2]/state.ρ)
          return Δt * (abs(u) + soundspeed(m.moisture, m.orientation, state,
                                            aux)) / Δx
  end
end

function Diffusive_CFL(m::AtmosModel, state::Vars, aux::Vars,
                               diffusive::Vars, Δx, Δt, direction=VerticalDirection())
  if direction isa VerticalDirection
    s_diff = diffusive.sdiff_v
          return Δt * (s_diff + soundspeed(m.moisture, m.orientation, state,
                                            aux)) / Δx^2
  elseif direction isa HorizontalDirection
    s_diff = diffusive.sdiff_h
          return Δt * (s_diff + soundspeed(m.moisture, m.orientation, state,
                                            aux)) / Δx^2
  end
end
end
