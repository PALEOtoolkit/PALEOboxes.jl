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


"""
    ReactionFluxPerturbColumn
 
 Add distributed flux `F` to a single column (use configuration file to rename output variable `F`)
 Flux `perturb_total` is distributed over region `zmin` < `z` < `zmax` of column index `icol`, with functional form
 `abs(z - distribute_z)^distribute_n`.
 """
Base.@kwdef mutable struct ReactionFluxPerturbColumn{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParIntVec("icol", [1],
            description="indices of column(s) to apply to"),
        PB.ParDouble("zmin", 0.0, units="m",
            description="lower z of region to apply flux to"),
        PB.ParDouble("zmax", 1e5, units="m",
            description="upper z of region to apply flux to"),
        PB.ParDouble("distribute_z", 1e5, units="m",
            description="flux scales as abs(z - distribute_z)^distribute_n"),
        PB.ParDouble("distribute_n", 1.0,
            description="flux scales as abs(z - distribute_z)^distribute_n"),

        PB.ParDouble("perturb_total", 1.0, units="mol yr-1",
            description="total flux to add"),
    )
   
end


function PB.register_methods!(rj::ReactionFluxPerturbColumn)

    vars = [
        PB.VarDepStateIndep("zmid",                 "m",        "mean z coord of box"),
        PB.VarDepStateIndep("volume",               "m^3",      "volume of box"),    
        PB.VarContrib(      "F",                    "mol yr-1", "distributed flux perturbation"),
        PB.VarContrib(      "(F_total)",            "mol yr-1", "total distributed flux perturbation"),
        PB.VarProp(         "%reaction%FApplied",   "mol yr-1", "flux perturbation applied, for diagnostic output"),    
    ]
    
    PB.add_method_do!(
        rj,
        do_flux_perturb_column,
        (PB.VarList_namedtuple(vars), ),
    )

    return nothing
end

function do_flux_perturb_column(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    # TODO this should be split into two methods, with scaleF and scaleFtot calculated in setup

    scaletot = 0.0
    for (icol, colindices) in cellrange.columns
        if icol in rj.pars.icol.v
            # first pass: define and scale functional form, using FApplied as a temporary variable
            scaletot = 0.0
            @inbounds for i in colindices
                if (vars.zmid[i] >= rj.pars.zmin.v) && (vars.zmid[i] <= rj.pars.zmax.v)
                    scaleF = vars.volume[i]*abs(vars.zmid[i] - rj.pars.distribute_z.v)^rj.pars.distribute_n.v
                    scaletot += scaleF
                    vars.FApplied[i] = scaleF
                else
                    vars.FApplied[i] = 0.0
                end
            end
        end
    end
    for (icol, colindices) in cellrange.columns
        if icol in rj.pars.icol.v
            # second pass: apply flux
            @inbounds for i in colindices
                if (vars.zmid[i] >= rj.pars.zmin.v) && (vars.zmid[i] <= rj.pars.zmax.v)
                    vars.FApplied[i] = vars.FApplied[i]/scaletot * rj.pars.perturb_total.v
                    vars.F[i] += vars.FApplied[i]
                end
            end
            vars.F_total[icol] += rj.pars.perturb_total.v
        end
    end
  
    return nothing
end


"""
    ReactionRedoxBalance
 
 Add `O2eqRequired` O2 equivalents as either O2 (if +ve) or CH4 (if -ve)
 """
Base.@kwdef mutable struct ReactionRedoxBalance{P} <: PB.AbstractReaction
    base::PB.ReactionBase
   
    pars::P = PB.ParametersTuple(
        PB.ParBool("include_ad", true,
            description="true to include automatic differentation in output fluxes"),
    )
end


function PB.register_methods!(rj::ReactionRedoxBalance)

    vars = [
        PB.VarDep(      "O2eqRequired", "mol O2eq yr-1", "O2 equivalents to add"),
        PB.VarContrib(  "F_O2",         "mol yr-1", "+ve redox flux"),
        PB.VarContrib(  "(F_O2_2)",         "mol yr-1", "+ve redox flux"),
        PB.VarContrib(  "F_CH4",        "mol yr-1", "-ve redox flux"),
        PB.VarContrib(  "(F_CH4_2)",    "mol yr-1", "-ve redox flux"),
    ]
 
    if rj.pars.include_ad.v
        ad_func = x -> x
    else
        ad_func = PB.value_ad
    end

    PB.add_method_do!(
        rj,
        do_flux_perturb_column,
        (PB.VarList_namedtuple(vars), ),
        p=ad_func,
    )

    return nothing
end

function do_redox_balance(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    ad_func = m.p

    for i in cellrange.indices
        O2eq = ad_func(vars.O2eqRequired[i])
        flux_O2 = zero(O2eq)
        flux_CH4 = zero(O2eq)

        if vars.O2eqRequired[i] > 0
            flux_O2 = O2eq
        else
            flux_CH4 = -0.5*O2eq
        end

        vars.F_O2[i] += flux_O2
        PB.add_if_available(vars.F_O2_2, i, flux_O2)
        vars.F_CH4[i] += flux_CH4
        PB.add_if_available(vars.F_CH4_2, i, flux_CH4)

    end
  
    return nothing
end


"Install create_reactionXXX factories when module imported"
function __init__()    
    PB.add_reaction_factory(ReactionFluxPerturb)
    PB.add_reaction_factory(ReactionRestore)
    PB.add_reaction_factory(ReactionFluxPerturbColumn)
    PB.add_reaction_factory(ReactionRedoxBalance)

    return nothing
end

end
