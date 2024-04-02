module Forcings

import PALEOboxes as PB

using ..DocStrings

"""
    ReactionForceInterp
 
Provide a scalar Property `F`, linearly interpolated from a table of values vs time `tforce`.

The table of values is set by parameters `force_times` and `force_values`.

The input time Variable is `tforce`, with default linking to the `global.tforce` Variable.

Use the configuration file to rename the output variable `F` (and if necessary, the input Variable `tforce`).

Set `extrapolate = "extrapolate"` to use `extrapolate_before`, `extrapolate_after` to set constant values
when 'tforce' is out-of-range of `force_times` (or alternatively, set `extrapolate = "throw"` and use 
guard values for `force_times` at -1e30, 1e30).

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
        PB.ParString("extrapolate", "throw", allowed_values=["throw", "constant"],
            description="behaviour if tforce is out of range"),
        PB.ParDouble("extrapolate_before", NaN,
            description="value to use if 'extrapolate=constant' and tforce < first(perturb_times)"),
        PB.ParDouble("extrapolate_after", NaN,
            description="value to use if 'extrapolate=constant' and tforce > last(perturb_times)"),
    )
 
end

function PB.register_methods!(rj::ReactionForceInterp)
    @info "ReactionForceInterp.register_methods! $(PB.fullname(rj))"

    vars = [
        PB.VarDepScalar("global.tforce",    "yr",   "historical time at which to apply forcings, present = 0 yr"),
        PB.VarPropScalar("F",               "",     "interpolated forcing"),
    ]

    PB.add_method_do!(rj, do_forceinterp, (PB.VarList_namedtuple(vars), ))

    return nothing
end

function do_forceinterp(m::PB.ReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    tforce              = vars.tforce[]    

    idx_after = searchsortedfirst(pars.force_times.v, tforce)
    if tforce == first(pars.force_times.v) # extra check as searchsortedfirst uses >=
        idx_after += 1
    end

    if idx_after == 1        
        (pars.extrapolate[] == "constant") || error("tforce $tforce out-of-range, < first(force_times) = $(first(pars.force_times.v))")
        vars.F[] = pars.extrapolate_before[]
    elseif idx_after > length(pars.force_times.v)
        (pars.extrapolate[] == "constant") || error("tforce $tforce out-of-range, > last(force_times) = $(last(pars.force_times.v))")
        vars.F[] = pars.extrapolate_after[]
    else
        # linearly interpolate
        t_l, t_h = pars.force_times.v[idx_after-1], pars.force_times.v[idx_after]

        x_l = (t_h-tforce)/(t_h - t_l)
        x_h = (tforce-t_l)/(t_h - t_l)

        vars.F[] = x_l*pars.force_values.v[idx_after-1] + x_h*pars.force_values.v[idx_after]
    end
   
    return nothing
end



end # module
