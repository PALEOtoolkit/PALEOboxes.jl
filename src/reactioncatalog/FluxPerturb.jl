module FluxPerturb

import Interpolations


import PALEOboxes as PB

"""
    ReactionFluxPerturb
 
 Add scalar flux `F` interpolated from parameters (use configuration file to rename output variable `F`)
 NB: no extrapolation ! (so eg set guard values at -1e30, 1e30)
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
end

function PB.register_methods!(rj::ReactionFluxPerturb)
    
    @info "ReactionFluxPerturb.register_methods! $(PB.fullname(rj))"

    PB.setfrozen!(rj.pars.field_data)

    vars = [
        PB.VarDepScalar(        "global.tforce",        "yr",       "time at which to apply perturbation"),
        PB.VarContribScalar(    "F",                    "mol yr-1", "interpolated flux perturbation",
            attributes=(:field_data=>rj.pars.field_data.v,)),
        PB.VarPropScalar(       "%reaction%FApplied",   "mol yr-1", "flux perturbation applied, for diagnostic output",
            attributes=(:field_data=>rj.pars.field_data.v,)),
    ]

    PB.add_method_do!(
        rj,
        do_flux_perturb,
        (PB.VarList_namedtuple(vars), ), 
        preparefn=prepare_do_flux_perturb,
    )

    return nothing
end

# We need to allow parameters to be changed after 'register_methods!'.
# But Interpolations encodes lengths of vectors in its type, which means
# it can't be added to ReactionMethod p and then changed later (without type instability).
# So in order to keep type stability, calculate it here and add to end of variables Tuple.
function prepare_do_flux_perturb(m::PB.ReactionMethod, (vars, ))
    rj = m.reaction

    IsotopeType = rj.pars.field_data.v

    interp_Ftotal = Interpolations.LinearInterpolation(
        rj.pars.perturb_times.v, 
        rj.pars.perturb_totals.v,
        extrapolation_bc = Interpolations.Throw()
    )

    if IsotopeType <: PB.AbstractIsotopeScalar
        interp_Fdelta = Interpolations.LinearInterpolation(
            rj.pars.perturb_times.v, 
            rj.pars.perturb_deltas.v,
            extrapolation_bc = Interpolations.Throw()
        )
    else
        interp_Fdelta = nothing
    end

    return (vars, IsotopeType, interp_Ftotal, interp_Fdelta)
end

function do_flux_perturb(
    m::PB.ReactionMethod, 
    (vars, IsotopeType, interp_Ftotal, interp_Fdelta),
    cellrange::PB.AbstractCellRange,
    deltat
)

    tforce              = vars.tforce[]    
    F                   = @PB.isotope_totaldelta(IsotopeType, interp_Ftotal(tforce), interp_Fdelta(tforce))

    vars.F[]            += F
    vars.FApplied[]     = F
  
    return nothing
end


"""
    ReactionRestore
 
 Adds `RestoringFlux` in response to discrepancy between Variable `WatchLevel` and Parameter `RequiredLevel`
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



"Install create_reactionXXX factories when module imported"
function __init__()    
    PB.add_reaction_factory(ReactionFluxPerturb)
    PB.add_reaction_factory(ReactionRestore)

    return nothing
end

end
