Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[]) # JuliaLang/julia/pull/28625

using CLIMA, Documenter, Literate

# TODO: Add generated examples back
# include("generate.jl")

GENERATED_BL_EXAMPLES = [
    joinpath("examples", "DGmethods_old", "generated", f)
    for
    f in ("ex_001_periodic_advection.md", "ex_002_solid_body_rotation.md")
]
GENERATED_BL_EXAMPLES = filter!(x -> isfile(x), GENERATED_BL_EXAMPLES)

# Literate examples created in their own literate_example directory for now
const examples_directory = joinpath(@__DIR__, "..", "literate_examples")
const output_directory = joinpath(@__DIR__, "src", "generated")
# read the examples in the directory
examples = readdir(examples_directory)

for example in examples
    example_filepath = joinpath(examples_directory, example)
    Literate.markdown(example_filepath, output_directory, documenter = true)
end

pages = Any[
    "Home" => "index.md",
    "Common" => Any[
        "MoistThermodynamics" => "Common/MoistThermodynamics.md",
        "SurfaceFluxes" => "Common/SurfaceFluxes.md",
    ],
    "Utilites" => Any[],
    "Atmos" => Any[
        "EDMF Equations" => "Atmos/EDMFEquations.md",
        "Microphysics" => "Atmos/Microphysics.md",
        "AtmosModel" => "Atmos/Model/AtmosModel.md",
        "Turbulence" => "Atmos/Model/turbulence.md",
        "Tracers" => "Atmos/Model/tracers.md",
    ],
    "Diagnostics" => "Diagnostics.md",
    "Numerics" => Any[
        "ODESolvers" => "Numerics/ODESolvers.md",
        "LinearSolvers" => "Numerics/LinearSolvers.md",
        "Mesh" => "Numerics/Mesh.md",
        "Arrays" => "Numerics/Arrays.md",
        "DGmethods_old" => "Numerics/DGmethods_old.md",
    ],
    "InputOutput.md",
    "Flow chart" => "FlowChart.md",
    "Examples" => [
        "Conjugate Gradient" => "generated/example_cg.md",
        "Notes on Literate" => "generated/literate_markdown.md",
    ],
    "Developer docs" => Any[
        "Contribution Guides" => Any[
            "ContributionGuides/CONTRIBUTING.md",
            "ContributionGuides/IterativeSolvers.md",
        ],
        "CodingConventions.md",
        "AcceptableUnicode.md",
        "VariableList.md",
        "DiagnosticVariables.md",
    ],
]

if !isempty(GENERATED_BL_EXAMPLES)
    push!(
        pages,
        "Balance Law Examples" =>
            ["BalanceLawOverview.md", GENERATED_BL_EXAMPLES...],
    )
end


makedocs(
    sitename = "CLIMA",
    doctest = false,
    strict = false,
    linkcheck = false,
    checkdocs = :exports,
    # checkdocs = :all,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(),
            ),
        )),
        # prettyurls = !("local" in ARGS),
        # canonical = "https://climate-machine.github.io/CLIMA/stable/",
    ),
    clean = true,
    modules = [Documenter, CLIMA],
    pages = pages,
)

# make sure there are no *.vtu files left around from the build
p = joinpath(@__DIR__, "build", "examples", "DGmethods_old", "generated")
if ispath(p)
    cd(p) do
        foreach(file -> endswith(file, ".vtu") && rm(file), readdir())
    end
end

deploydocs(
    repo = "github.com/climate-machine/CLIMA.git",
    target = "build",
    push_preview = true,
)
