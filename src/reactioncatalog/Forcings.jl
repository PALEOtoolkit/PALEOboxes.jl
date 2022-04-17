module Forcings

import PALEOboxes as PB

"""
    ReactionForceInterp
 
 `F` forcing interpolated from parameters (use configuration file to rename output variable `F`)
 NB: no extrapolation ! (so eg set guard values at -1e30, 1e30)
 """
Base.@kwdef mutable struct ReactionForceInterp{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("force_times", [-1e30, 1e30], units="yr", 
            description="interpolated forcing times"),
        PB.ParDoubleVec("force_values",[1.0, 1.0], units="",
            description="interpolated forcing values"),
    )
 
end

function PB.register_methods!(rj::ReactionForceInterp)
    @info "ReactionForceInterp.register_methods! $(PB.fullname(rj))"

    vars = [
        PB.VarDepScalar("global.tforce",    "yr",   "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("F",               "",     "interpolated forcing"),
    ]

    # placeholder - interp_F is calculated in prepare_forceinterp
    interp_F = PB.LinInterp([NaN])
    PB.add_method_do!(
        rj,
        do_forceinterp,
        (PB.VarList_namedtuple(vars), ),
        p=interp_F,
        preparefn=prepare_forceinterp,
    )

    return nothing
end

function prepare_forceinterp(m::PB.ReactionMethod, vardata)
    rj = m.reaction

    interp_F = PB.LinInterp(rj.pars.force_times.v)

    m.p = interp_F

    return vardata
end

function do_forceinterp(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    interp_F = m.p

    vars.F[] = PB.interp(interp_F, vars.tforce[], rj.pars.force_values.v)
   
    return nothing
end




"Install create_reactionXXX factories when module imported"
function __init__()
    PB.add_reaction_factory(ReactionForceInterp)
    return nothing
end


end # module
