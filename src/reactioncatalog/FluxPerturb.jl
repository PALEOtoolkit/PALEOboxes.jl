module FluxPerturb

import PALEOboxes as PB

using ..DocStrings


"""
    ReactionFluxPerturb
 
Provide a scalar flux `F`, linearly interpolated from a table of values vs time `tforce`.

The table of values is set by parameters `perturb_times` and `perturb_totals`, and optionally `perturb_deltas`.

The input time Variable is `tforce`, with default linking to the `global.tforce` Variable.

Use the configuration file to rename output variable `F` (and if necessary, the input Variable `tforce`).

Set `extrapolate = "extrapolate"` to use `extrapolate_before`, `extrapolate_after` to set constant values
when 'tforce' is out-of-range of `perturb_times` (or alternatively, set `extrapolate = "throw"` and use 
guard values for `perturb_times` at -1e30, 1e30).

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
        PB.ParString("extrapolate", "throw", allowed_values=["throw", "constant"],
            description="behaviour if tforce is out of range"),
        PB.ParDouble("extrapolate_before", NaN,
            description="value to use if 'extrapolate=constant' and tforce < first(perturb_times)"),
        PB.ParDouble("extrapolate_before_delta", NaN,
            description="value to use if 'extrapolate=constant' and tforce < first(perturb_times)"),
        PB.ParDouble("extrapolate_after", NaN,
            description="value to use if 'extrapolate=constant' and tforce > last(perturb_times)"),
        PB.ParDouble("extrapolate_after_delta", NaN,
            description="value to use if 'extrapolate=constant' and tforce > last(perturb_times)"),
    )

end


function PB.register_methods!(rj::ReactionFluxPerturb)
    
    @info "ReactionFluxPerturb.register_methods! $(PB.fullname(rj))"

    IsotopeType = rj.pars.field_data[]
    PB.setfrozen!(rj.pars.field_data)

    vars = [
        PB.VarDepScalar(        "global.tforce",        "yr",       "time at which to apply perturbation"),
        PB.VarContribScalar(    "F",                    "mol yr-1", "interpolated flux perturbation",
            attributes=(:field_data=>IsotopeType,)),
        PB.VarPropScalar(       "%reaction%FApplied",   "mol yr-1", "flux perturbation applied, for diagnostic output",
            attributes=(:field_data=>IsotopeType,)),
    ]

    PB.add_method_do!(
        rj,
        do_flux_perturb,
        (PB.VarList_namedtuple(vars), ),
        p=IsotopeType,
    )

    return nothing
end


function do_flux_perturb(
    m::PB.ReactionMethod,
    pars, 
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    IsotopeType = m.p

    tforce              = vars.tforce[]    

    f = zero(first(pars.perturb_totals.v)*vars.tforce[]) # ensure zero is of correct type for perturb_totals, tforce
    f_delta = zero(first(pars.perturb_deltas.v)*vars.tforce[])

    idx_after = searchsortedfirst(pars.perturb_times.v, tforce)
    if tforce == first(pars.perturb_times.v) # extra check as searchsortedfirst uses >=
        idx_after += 1
    end

    if idx_after == 1        
        (pars.extrapolate[] == "constant") || error("tforce $tforce out-of-range, < first(perturb_times) = $(first(pars.perturb_times.v))")
        f += pars.extrapolate_before[]
        if IsotopeType != PB.ScalarData
            f_delta += pars.extrapolate_before_delta[]
        end
    elseif idx_after > length(pars.perturb_times.v)
        (pars.extrapolate[] == "constant") || error("tforce $tforce out-of-range, > last(perturb_times) = $(last(pars.perturb_times.v))")
        f += pars.extrapolate_after[]
        if IsotopeType != PB.ScalarData
            f_delta += pars.extrapolate_after_delta[]
        end
    else
        # linearly interpolate
        t_l, t_h = pars.perturb_times.v[idx_after-1], pars.perturb_times.v[idx_after]

        x_l = (t_h-tforce)/(t_h - t_l)
        x_h = (tforce-t_l)/(t_h - t_l)

        f += x_l*pars.perturb_totals.v[idx_after-1] + x_h*pars.perturb_totals.v[idx_after]
        if IsotopeType != PB.ScalarData
            f_delta += x_l*pars.perturb_deltas.v[idx_after-1] + x_h*pars.perturb_deltas.v[idx_after]
        end
    end

    F = @PB.isotope_totaldelta(IsotopeType, f, f_delta)

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
    IsotopeType = rj.pars.field_data[]

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

    requiredlevelisotope = @PB.isotope_totaldelta(IsotopeType, rj.pars.RequiredLevel[], rj.pars.RequiredDelta[])

    restoreflux = -(vars.WatchLevel[] - requiredlevelisotope)/rj.pars.trestore[]
   
    if !rj.pars.source_only[] || PB.get_total(restoreflux) > 0
        vars.RestoringApplied[] = restoreflux
        vars.RestoringFlux[]    += restoreflux
    else
        vars.RestoringApplied[] = 0.0*restoreflux
        vars.RestoringFlux[]    += 0.0*restoreflux # multiply by zero so AD sparsity detection sees the restoring term
    end

    return nothing
end

end
