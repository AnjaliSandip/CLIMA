#### Heat equation tutorial

#=
Heat equation

```
∂T
-- + ∇⋅(-α ∇T) = 0
∂t
```

Or,

```
∂T
-- + ∇⋅(F(T,t)) = 0
∂t
```

where
 - `α` is the thermal conductivity (W/(m K))
 - `F(T,t)` is prescribed in `flux_diffusive!`

=#

####
#### 1) Import/Export Needed Functions
####
println("1) Import/Export Needed Functions")

# Load necessary CliMA subroutines
using MPI
using Test
using CLIMA
using Logging
using Printf
using NCDatasets
using LinearAlgebra
using OrderedCollections
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Writers
using CLIMA.VTK
using CLIMA.Mesh.Elements: interpolationmatrix
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.GenericCallbacks: EveryXWallTimeSeconds, EveryXSimulationSteps
using CLIMA.GenericCallbacks
using CLIMA.ODESolvers

using Interpolations
using DelimitedFiles

using Plots

ENV["CLIMA_GPU"] = "false"

# Initialize CliMA
CLIMA.init()

FT = Float64

# Change output directory and save plots there
output_dir = joinpath(dirname(dirname(pathof(CLIMA))), "output", "land")
mkpath(output_dir)

# Add soil model
include("heat_equation_model.jl")

####
#### Include helper and plotting functions (to be refactored/moved into CLIMA src)
####

function get_z(grid)
    # TODO: this currently uses some internals: provide a better way to do this
    return reshape(grid.vgeo[(1:(N+1)^2:(N+1)^3),CLIMA.Mesh.Grids.vgeoid.x3id,:],:)*100
end
include(joinpath("..","helper_funcs.jl"))
include(joinpath("..","plotting_funcs.jl"))

####
#### 2) Set up domain
####
println("2) Set up domain...")

# NOTE: this is using 5 vertical elements, each with a 5th degree polynomial,
# giving an approximate resolution of 5cm
velems = collect(0:10) # Elements at: [0.0 -0.2 -0.4 -0.6 -0.8 -1.0] (m)
velems = velems / 10
N = 5 # Order of polynomial function between each element

# Set domain using Stached Brick Topology
topl = StackedBrickTopology(MPI.COMM_WORLD, (0.0:1,0.0:1,velems);
    periodicity = (true,true,false),
    boundary=((0,0),(0,0),(1,2)))
grid = DiscontinuousSpectralElementGrid(topl, FloatType = Float64, DeviceArray = Array, polynomialorder = N)

# Create Model struct
m = HeatModel(
    # Define heat capacity of soil
    ρc = (state, aux, t) ->  1,

    # Define thermal conductivity of soil
    α  = (state, aux, t) ->  0.01,

    # Define initial temperature of soil
    initialT = (aux, t) -> (273.15 + 22.0),

    # Define surface boundary condition
    surfaceT = (state, aux, t) -> 300.0 # replace with T_data
)

# Set up DG scheme
dg = DGModel( #
  m, # "PDE part"
  grid,
  CentralNumericalFluxNonDiffusive(), # penalty terms for discretizations
  CentralNumericalFluxDiffusive(),
  CentralNumericalFluxGradient())

# Minimum spatial and temporal steps
Δ = min_node_distance(grid)

given_CFL = 0.08
CFL_bound = given_CFL*Δ^2 / m.α(1,1,1)
dt = CFL_bound

####
#### 3) Define variables for simulation
####
println("3) Define variables for simulation...")

# Define time variables
const minute = 60
const hour = 60*minute
const day = 24*hour
const n_outputs = 5
const timeend = 20

# Output frequency:
const every_x_simulation_time = ceil(Int, timeend/n_outputs)


####
#### 4) Prep ICs, and time-stepper and output configurations
####
println("4) Prep ICs, and time-stepper and output configurations...")

# state variable
Q = init_ode_state(dg, Float64(0))

# initialize ODE solver
lsrk = LSRK54CarpenterKennedy(dg, Q; dt = dt, t0 = 0)

# Plot initial state
p = get_plot(grid, Q, dg.auxstate, 0)
export_plots(p, joinpath(output_dir, "initial_state.png"))

mkpath(output_dir)

plots = []
dims = OrderedDict("z" => collect(get_z(grid)))
# run for 8 days (hours?) to get to steady state

output_data = DataFile(joinpath(output_dir, "output_data"))

step = [0]
stcb = GenericCallbacks.EveryXSimulationTime(every_x_simulation_time, lsrk) do (init = false)
  state_vars = get_vars_from_stack(grid, Q, m, vars_state; exclude=["θi"])
  aux_vars = get_vars_from_stack(grid, dg.auxstate, m, vars_aux; exclude=["z","ψ"])
  all_vars = OrderedDict(state_vars..., aux_vars...)
  write_data(NetCDFWriter(), output_data(step[1]), dims, all_vars, gettime(lsrk))
  step[1]+=1
  nothing
end


####
#### 5) Solve the equations
####
println("5) Solve the equations...")

solve!(Q, lsrk; timeend=timeend, callbacks=(stcb,))

#####
##### 6) Post-processing
#####
println("6) Post-processing...")

all_data = collect_data(output_data, step[1])

# To get "T" at timestep 0:
# all_data[0]["T"][:]

plot_solution(all_data, ("ρcT",), joinpath(output_dir, "solution_vs_time.png"))


