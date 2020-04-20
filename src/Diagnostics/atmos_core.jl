using ..Atmos
using ..Atmos: thermo_state, turbulence_tensors
using ..Mesh.Topologies
using ..Mesh.Grids
using ..MoistThermodynamics
using LinearAlgebra

include("atmos_common.jl")

"""
    atmos_core_init(bl, currtime)

Initialize the 'AtmosCore' diagnostics group.
"""
function atmos_core_init(dgngrp::DiagnosticsGroup, currtime)
    atmos_collect_onetime(Settings.mpicomm, Settings.dg, Settings.Q)

    return nothing
end

include("thermo.jl")

"""
    atmos_core_collect(bl, currtime)

Perform a global grid traversal to compute various diagnostics.
"""
function atmos_core_collect(dgngrp::DiagnosticsGroup, currtime)
    mpicomm = Settings.mpicomm
    dg = Settings.dg
    Q = Settings.Q
    mpirank = MPI.Comm_rank(mpicomm)
    current_time = string(currtime)

    # extract grid information
    bl = dg.balancelaw
    grid = dg.grid
    topology = grid.topology
    N = polynomialorder(grid)
    Nq = N + 1
    Nqk = dimensionality(grid) == 2 ? 1 : Nq
    npoints = Nq * Nq * Nqk
    nrealelem = length(topology.realelems)
    nvertelem = topology.stacksize
    nhorzelem = div(nrealelem, nvertelem)

    # get the state, auxiliary and geo variables onto the host if needed
    if Array ∈ typeof(Q).parameters
        localQ = Q.realdata
        localaux = dg.auxstate.realdata
        localvgeo = grid.vgeo
        localdiff = dg.diffstate.realdata
    else
        localQ = Array(Q.realdata)
        localaux = Array(dg.auxstate.realdata)
        localvgeo = Array(grid.vgeo)
        localdiff = Array(dg.diffstate.realdata)
    end
    FT = eltype(localQ)

    zvals = AtmosCollected.zvals

    core_repdvsr = zeros(FT, Nqk * nvertelem)
    thermo_array =
        [zeros(FT, num_thermo(bl, FT)) for _ in 1:npoints, _ in 1:nrealelem]
    ql_w_gt_0 = [zeros(FT, (Nq * Nq * nhorzelem)) for _ in 1:(Nqk * nvertelem)]
    @visitQ nhorzelem nvertelem Nqk Nq begin
        evk = Nqk * (ev - 1) + k

        state = extract_state(dg, localQ, ijk, e)
        aux = extract_aux(dg, localaux, ijk, e)
        MH = localvgeo[ijk, grid.MHid, e]

        thermo = thermo_vars(bl, thermo_array[ijk, e])
        compute_thermo!(bl, state, aux, thermo)

        if thermo.moisture.q_liq > 0 && state.w > 0
            idx = (Nq * Nq * (eh - 1)) + (Nq * (j - 1)) + i
            ql_w_gt_0[evk][idx] = one(FT)
            core_repdvsr[evk] += MH
        end
    end
    MPI.Allreduce!(core_repdvsr, +, mpicomm)
    core_frac = zeros(FT, Nqk * nvertelem)
    for evk in 1:(Nqk * nvertelem)
        tot_ql_w_gt_0 = MPI.Reduce(sum(ql_w_gt_0[evk]), +, 0, mpicomm)
        tot_horz = MPI.Reduce(length(ql_w_gt_0[evk]), +, 0, mpicomm)
        if mpirank == 0
            core_frac[evk] = tot_ql_w_gt_0 / tot_horz
        end
    end

    if mpirank == 0
        dprefix = @sprintf(
            "%s_%s_%s_num%04d",
            dgngrp.out_prefix,
            dgngrp.name,
            Settings.starttime,
            dgngrp.num
        )
        dfilename = joinpath(Settings.output_dir, dprefix)
    end

    MPI.Barrier(mpicomm)
    return nothing
end # function collect

function atmos_core_fini(dgngrp::DiagnosticsGroup, currtime) end
