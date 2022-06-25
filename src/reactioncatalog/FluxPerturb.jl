module FluxPerturb

import Interpolations


import PALEOboxes as PB

using ..DocStrings

const LINEARINTERPOLATION_TEMPLATE = Interpolations.LinearInterpolation(
    [0.0, 1.0], 
    [NaN, NaN],
    extrapolation_bc = Interpolations.Throw()
)


"""
    ReactionFluxPerturb
 
Provide a scalar flux `F`, linearly interpolated from a table of values vs time `tforce`.

The table of values is set by parameters `perturb_times` and `perturb_totals`, and optionally `perturb_deltas`.

The input time Variable is `tforce`, with default linking to the `global.tforce` Variable.

Use the configuration file to rename output variable `F` (and if necessary, the input Variable `tforce`).

NB: no extrapolation ! (so eg set guard values for `perturb_times` at -1e30, 1e30)

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionFluxPerturb{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),

        PB.ParDoubleVec("perturb_times", [-1e30, 1e30], units="yr",
            description="interpolated perturbation times"),
        PB.ParDoubleVec("perturb_totals",[1.0, 1.0], units="mol yr-1",
            description="interpolated perturbation totals"),
        PB.ParDoubleVec("perturb_deltas",[0.0, 0.0], units="per mil",
            description="interpolated perturbation deltas"),
    )

    interp_Ftotal::typeof(LINEARINTERPOLATION_TEMPLATE) = LINEARINTERPOLATION_TEMPLATE
    interp_Fdelta::typeof(LINEARINTERPOLATION_TEMPLATE) = LINEARINTERPOLATION_TEMPLATE
end


function PB.register_methods!(rj::ReactionFluxPerturb)
    
    @info "ReactionFluxPerturb.register_methods! $(PB.fullname(rj))"

    IsotopeType = rj.pars.field_data.v
    PB.setfrozen!(rj.pars.field_data)

    vars = [
        PB.VarDepScalar(        "global.tforce",        "yr",       "time at which to apply perturbation"),
        PB.VarContribScalar(    "F",                    "mol yr-1", "interpolated flux perturbation",
            attributes=(:field_data=>IsotopeType,)),
        PB.VarPropScalar(       "%reaction%FApplied",   "mol yr-1", "flux perturbation applied, for diagnostic output",
            attributes=(:field_data=>IsotopeType,)),
    ]

    PB.add_method_setup!(
        rj,
        setup_flux_perturb,
        Tuple{}(), 
        p=IsotopeType,
    )

    PB.add_method_do!(
        rj,
        do_flux_perturb,
        (PB.VarList_namedtuple(vars), ),
        p=IsotopeType,
    )

    return nothing
end

function setup_flux_perturb(
    m::PB.ReactionMethod, 
    _,
    cellrange::PB.AbstractCellRange,
    attribute_name,
)
    attribute_name == :setup || return

    @info "$(PB.fullname(m)):"

    rj = m.reaction
    IsotopeType = m.p

    rj.interp_Ftotal = Interpolations.LinearInterpolation(
        rj.pars.perturb_times.v, 
        rj.pars.perturb_totals.v,
        extrapolation_bc = Interpolations.Throw()
    )

    if IsotopeType <: PB.AbstractIsotopeScalar
        rj.interp_Fdelta = Interpolations.LinearInterpolation(
            rj.pars.perturb_times.v, 
            rj.pars.perturb_deltas.v,
            extrapolation_bc = Interpolations.Throw()
        )
    end

    return nothing
end

function do_flux_perturb(
    m::PB.ReactionMethod, 
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    IsotopeType = m.p

    tforce              = vars.tforce[]    
    F                   = @PB.isotope_totaldelta(IsotopeType, rj.interp_Ftotal(tforce), rj.interp_Fdelta(tforce))

    vars.F[]            += F
    vars.FApplied[]     = F
  
    return nothing
end


"""
    ReactionRestore
 
Adds `RestoringFlux` in response to discrepancy between Variable `WatchLevel` and Parameter `RequiredLevel`

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionRestore{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),

        PB.ParDouble("RequiredLevel", 0.0, units="mol",
            description="target value of WatchLevel"),
        PB.ParDouble("trestore", 1.0, units="yr",
            description="restoring timescale"),
        PB.ParDouble("RequiredDelta", 0.0, units="per mil", 
            description="target value of WatchLevel delta"),

        PB.ParBool("source_only", false, 
            description="false to allow input and output, true to allow input only"),
    )

end

function PB.register_methods!(rj::ReactionRestore)
  
    PB.setfrozen!(rj.pars.field_data)
    IsotopeType = rj.pars.field_data.v   

    vars = [
        PB.VarDepScalar(    "WatchLevel",                   "mol",      "level to observe and restore",
            attributes=(:field_data=>IsotopeType,)),
        PB.VarContribScalar("RestoringFlux",                "mol yr-1", "restoring flux",
            attributes=(:field_data=>IsotopeType,)),
        PB.VarPropScalar(   "%reaction%RestoringApplied",   "mol yr-1", "restoring flux for diagnostic output",
            attributes=(:field_data=>IsotopeType,)),    
    ]

  
    PB.add_method_do!(
        rj,
        do_restore,
        (PB.VarList_namedtuple(vars), ), 
        p=IsotopeType,
    )

    return nothing
end

function do_restore(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction    
    IsotopeType = m.p

    requiredlevelisotope = @PB.isotope_totaldelta(IsotopeType, rj.pars.RequiredLevel.v, rj.pars.RequiredDelta.v)

    restoreflux = -(vars.WatchLevel[] - requiredlevelisotope)/rj.pars.trestore.v
   
    if !rj.pars.source_only.v || PB.get_total(restoreflux) > 0
        vars.RestoringApplied[] = restoreflux
        vars.RestoringFlux[]    += restoreflux
    else
        vars.RestoringApplied[] = 0.0*restoreflux
        vars.RestoringFlux[]    += 0.0*restoreflux # multiply by zero so AD sparsity detection sees the restoring term
    end

    return nothing
end

end
