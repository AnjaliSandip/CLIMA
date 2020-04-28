"""
    RemainderModel(main::BalanceLaw, subcomponents::Tuple)

Compute the "remainder" contribution of the `main` model, after subtracting
`subcomponents`.

Currently only the `flux_nondiffusive!` and `source!` are handled by the
remainder model
"""
struct RemainderModel{M, S} <: BalanceLaw
    main::M
    subs::S
end

# Inherit most of the functionality from the main model
vars_state(rem::RemainderModel, FT) = vars_state(rem.main, FT)

vars_gradient(rem::RemainderModel, FT) = vars_gradient(rem.main, FT)

vars_diffusive(rem::RemainderModel, FT) = vars_diffusive(rem.main, FT)

vars_aux(rem::RemainderModel, FT) = vars_aux(rem.main, FT)

vars_integrals(rem::RemainderModel, FT) = vars_integrals(rem.main, FT)

vars_reverse_integrals(rem::RemainderModel, FT) = vars_integrals(rem.main, FT)

vars_gradient_laplacian(rem::RemainderModel, FT) =
    vars_gradient_laplacian(rem.main, FT)

vars_hyperdiffusive(rem::RemainderModel, FT) = vars_hyperdiffusive(rem.main, FT)

update_aux!(dg::DGModel, rem::RemainderModel, args...) =
    update_aux!(dg, rem.main, args...)

update_aux_diffusive!(dg::DGModel, rem::RemainderModel, args...) =
    update_aux_diffusive!(dg, rem.main, args...)

integral_load_aux!(rem::RemainderModel, args...) =
    integral_load_aux!(rem.main, args...)

integral_set_aux!(rem::RemainderModel, args...) =
    integral_set_aux!(rem.main, args...)

reverse_integral_load_aux!(rem::RemainderModel, args...) =
    reverse_integral_load_aux!(rem.main, args...)

reverse_integral_set_aux!(rem::RemainderModel, args...) =
    reverse_integral_set_aux!(rem.main, args...)

hyperdiffusive!(rem::RemainderModel, args...) =
    hyperdiffusive!(rem.main, args...)

flux_diffusive!(rem::RemainderModel, args...) =
    flux_diffusive!(rem.main, args...)

gradvariables!(rem::RemainderModel, args...) = gradvariables!(rem.main, args...)

diffusive!(rem::RemainderModel, args...) = diffusive!(rem.main, args...)

boundary_state!(nf, rem::RemainderModel, x...) =
    boundary_state!(nf, rem.main, x...)

normal_boundary_flux_diffusive!(nf, rem::RemainderModel, args...) =
    normal_boundary_flux_diffusive!(nf, rem.main, args...)

init_aux!(rem::RemainderModel, args...) = init_aux!(rem.main, args...)

init_state!(rem::RemainderModel, args...) = init_state!(rem.main, args...)

function wavespeed(rem::RemainderModel, nM, state::Vars, aux::Vars, t::Real)
    ref = aux.ref_state
    return wavespeed(rem.main, nM, state, aux, t) -
           sum(sub -> wavespeed(sub, nM, state, aux, t), rem.subs)
end

function flux_nondiffusive!(
    rem::RemainderModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    m = getfield(flux, :array)
    flux_nondiffusive!(rem.main, flux, state, aux, t)

    flux_s = similar(flux)
    m_s = getfield(flux_s, :array)

    for sub in rem.subs
        fill!(m_s, 0)
        flux_nondiffusive!(sub, flux_s, state, aux, t)
        m .-= m_s
    end
    nothing
end

function source!(
    rem::RemainderModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    m = getfield(source, :array)
    source!(rem.main, source, state, diffusive, aux, t, direction)

    source_s = similar(source)
    m_s = getfield(source_s, :array)

    for sub in rem.subs
        fill!(m_s, 0)
        source!(sub, source_s, state, diffusive, aux, t, direction)
        m .-= m_s
    end
    nothing
end
