using MPI, Test
include("../testhelpers.jl")

@testset "DGmethods" begin
  tests = [(1, "compressible_Navier_Stokes/rising_bubble-model.jl")
           (1, "compressible_Navier_Stokes/density_current-model.jl")
          ]

  runmpi(tests, @__FILE__)
end