module Forcings

import PALEOboxes as PB

using ..DocStrings

"""
    ReactionForceInterp
 
Provide a scalar Property `F`, linearly interpolated from a table of values vs time `tforce`.

The table of values is set by parameters `force_times` and `force_value`.

The input time Variable is `tforce`, with default linking to the `global.tforce` Variable.

Use the configuration file to rename the output variable `F` (and if necessary, the input Variable `tforce`).

NB: no extrapolation ! (so eg set guard values for `force_times` at -1e30, 1e30)

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
 """
Base.@kwdef mutable struct ReactionForceInterp{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDoubleVec("force_times", [-1e30, 1e30], units="yr", 
            description="interpolated forcing times"),
        PB.ParDoubleVec("force_values", [1.0, 1.0], units="",
            description="interpolated forcing values"),
    )
 
    # placeholder - interp_F is calculated in setup_forceinterp
    interp_F::PB.LinInterp = PB.LinInterp([NaN])
end

function PB.register_methods!(rj::ReactionForceInterp)
    @info "ReactionForceInterp.register_methods! $(PB.fullname(rj))"

    PB.add_method_setup!(rj, setup_forceinterp, Tuple{}())

    vars = [
        PB.VarDepScalar("global.tforce",    "yr",   "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("F",               "",     "interpolated forcing"),
    ]

    PB.add_method_do!(rj, do_forceinterp, (PB.VarList_namedtuple(vars), ))

    return nothing
end

function setup_forceinterp(m::PB.ReactionMethod, pars, _, cellrange::PB.AbstractCellRange, attribute_name)
    attribute_name == :setup || return

    @info "$(PB.fullname(m)):"

    rj = m.reaction

    rj.interp_F = PB.LinInterp(pars.force_times.v)

    return nothing
end

function do_forceinterp(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    vars.F[] = PB.interp(rj.interp_F, vars.tforce[], pars.force_values.v)
   
    return nothing
end



end # module
