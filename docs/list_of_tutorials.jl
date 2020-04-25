####
#### Defines list of tutorials given `tutorials_dir` and `generated_dir`
####

# generate tutorials
import Literate
using AutoPages: gather_pages, replace_reverse

tutorials_jl = [
    joinpath(root, f)
    for (root, dirs, files) in Base.Filesystem.walkdir(tutorials_dir)
    for f in files
]

filter!(x -> !occursin("topo.jl", x), tutorials_jl)                       # currently broken, TODO: Fix me!
filter!(x -> !occursin("dry_rayleigh_benard.jl", x), tutorials_jl)        # currently broken, TODO: Fix me!
filter!(x -> !occursin("ex_001_periodic_advection.jl", x), tutorials_jl)  # currently broken, TODO: Fix me!
filter!(x -> !occursin("ex_002_solid_body_rotation.jl", x), tutorials_jl) # currently broken, TODO: Fix me!
filter!(x -> !occursin("ex_003_acoustic_wave.jl", x), tutorials_jl)       # currently broken, TODO: Fix me!
filter!(x -> !occursin("ex_004_nonnegative.jl", x), tutorials_jl)         # currently broken, TODO: Fix me!
filter!(x -> !occursin("KinematicModel.jl", x), tutorials_jl)             # currently broken, TODO: Fix me!
filter!(x -> !occursin("ex_1_saturation_adjustment.jl", x), tutorials_jl) # currently broken, TODO: Fix me!
filter!(x -> !occursin("ex_2_Kessler.jl", x), tutorials_jl)               # currently broken, TODO: Fix me!

println("Building literate tutorials:")
for tutorial in tutorials_jl
    println("    $(tutorial)")
end

for tutorial in tutorials_jl
    gen_dir = joinpath(generated_dir, relpath(dirname(tutorial), tutorials_dir))
    input = abspath(tutorial)
    script = Literate.script(input, gen_dir)
    code = strip(read(script, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    Literate.markdown(input, gen_dir, postprocess = mdpost)
    # Literate.notebook(input, gen_dir, execute = true)
end

tutorials, _ = gather_pages(;
  filenames=relpath.(tutorials_jl, dirname(@__DIR__)),
  extension_filter=x->endswith(x, ".jl"),
  transform_extension=x->replace_reverse(x, ".jl" => ".md"; count=1),
  remove_first_level=true,
  transform_path=x->replace(x, "tutorials"=>"generated", count=1),
  )

# Allow flag to skip generated
# tutorials since this is by
# far the slowest part of the
# docs build.
if !generate_tutorials
    tutorials = Any[]
end
