# -*- coding: utf-8 -*-
module Reservoirs

import ...PALEOboxes as PB
using ...PALEOboxes: @public

using ..DocStrings

@public ReactionReservoirScalar, ReactionReservoir, ReactionReservoirTotal, ReactionReservoirConst, 
    ReactionReservoirForced, ReactionReservoirWellMixed, ReactionConst, ReactionScalarConst

"""
    ReactionReservoirScalar

A single scalar biogeochemical reservoir with optional paired isotope reservoir, for use in a 0D Domain 
(eg sedimentary or ocean reservoirs for COPSE [Bergman2004](@cite)).

Creates State and associated Variables, depending on parameter settings:
- `const=false`: usual case
  - `state_norm=false` create state variable `R` (units mol, with attribute `vfunction=VF_StateExplicit`)
    and `R_sms` (units mol yr-1, with attribute `vfunction=VF_Deriv`).
  - `state_norm=true` create state variable `R_solve` (`R` normalized by the values of attribute `R:norm_value`, with attribute `vfunction=VF_StateExplicit`)
    and `R_solve_sms` (units yr-1, with attribute `vfunction=VF_Deriv`).
- `const=true`: a constant value, create `R` (a Property), and `R_sms` (a Target)

In addition:
- a Property `R_norm` (normalized value) is always created.
- if parameter `field_data <: AbstractIsotopeScalar` (eg `IsotopeLinear`), a Property `R_delta` is created.

The local name prefix `R` should then be renamed using `variable_links:` in the configuration file.

# Initialisation
Initial and norm value is set in the `variable_attributes:` section in 
  the configuration file, using `R:initial_value`, `R:initial_delta`, and `R:norm_value`.

# Example configuration in .yaml file
                reservoir_P:  # 0D ocean Phosphorus
                    class: ReactionReservoirScalar
                    parameters:
                        # field_data: ScalarData        # change to IsotopeLinear to represent an isotope
                        # const: false                  # true to fix to constant value
                    variable_links:
                        R*: P                           # rename to represent Phosphorus
                    variable_attributes:
                        R:norm_value:           3.1e15  # mol 
                        R:initial_value:        3.1e15  # mol

# See also
[`ReactionReservoir`](@ref) (one value per cell for a spatially resolved Domain eg ocean),
[`ReactionReservoirWellMixed`](@ref) (one value for a whole spatially resolved Domain).

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReservoirScalar{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
        PB.ParBool("const", false,
            description="true to provide a constant value: R is not a state variable, fluxes in R_sms Variable are ignored"),
        PB.ParBool("state_norm", false,
            description="true to provide solver with normalized values"),
    )

    norm_value::Float64  = NaN
end


function PB.register_methods!(rj::ReactionReservoirScalar)

    do_vars = PB.VariableReaction[PB.VarPropScalar("R_norm", "", "scalar reservoir normalized")]
    if rj.pars.const[]
        R            = PB.VarPropScalar(      "R", "mol", "scalar constant reservoir", attributes=(:field_data =>rj.pars.field_data[],))
        push!(do_vars, R)

        PB.add_method_setup!(
            rj,
            setup_reactionreservoirscalar,
            (PB.VarList_fields([R]), PB.VarList_nothing() ),
        )

        R_sms        = PB.VarTargetScalar(     "R_sms", "mol yr-1", "scalar reservoir source-sinks", attributes=(:field_data =>rj.pars.field_data[],))      
        # sms variable not used by us, but must appear in a method to be linked and created
        PB.add_method_do_nothing!(rj, [R_sms])
    else
        if rj.pars.state_norm[]
            R           = PB.VarPropScalar("R", "mol", "scalar reservoir", attributes=(:field_data =>rj.pars.field_data[],))            
            R_solve     = PB.VarStateExplicitScalar("R_solve", "", "normalized scalar reservoir", attributes=(:field_data =>rj.pars.field_data[],))
            append!(do_vars, [R, R_solve])
            PB.add_method_setup!(
                rj,
                setup_reactionreservoirscalar,
                (PB.VarList_fields([R]), PB.VarList_fields([R_solve]) ),
            )

            R_sms       = PB.VarTargetScalar(     "R_sms", "mol yr-1", "scalar reservoir source-sinks", attributes=(:field_data =>rj.pars.field_data[],))
            R_solve_sms = PB.VarDerivScalar(     "R_solve_sms", "yr-1", "normalized scalar reservoir source-sinks", attributes=(:field_data =>rj.pars.field_data[],))
        
            PB.add_method_do!(
                rj,
                do_reactionreservoirscalar_sms,
                (PB.VarList_namedtuple([R_sms, R_solve_sms]), ),
            )
        else 
            R           = PB.VarStateExplicitScalar("R", "mol", "scalar reservoir", attributes=(:field_data =>rj.pars.field_data[],))
            push!(do_vars, R)
            PB.add_method_setup!(
                rj,
                setup_reactionreservoirscalar,
                (PB.VarList_fields([R]), PB.VarList_nothing() ),
            )

            R_sms       = PB.VarDerivScalar(     "R_sms", "mol yr-1", "scalar reservoir source-sinks", attributes=(:field_data =>rj.pars.field_data[],))
            # sms variable not used by us, but must appear in a method to be linked and created
            PB.add_method_do_nothing!(rj, [R_sms])
        end
        PB.setfrozen!(rj.pars.state_norm)
    end
    PB.setfrozen!(rj.pars.const)

    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        push!(do_vars, PB.VarPropScalar("R_delta", "per mil", "scalar reservoir isotope delta"))
    end
    PB.setfrozen!(rj.pars.field_data)

    PB.add_method_do!(
        rj,
        do_reactionreservoirscalar,
        (PB.VarList_namedtuple(do_vars), ),
    )

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function setup_reactionreservoirscalar(m::PB.AbstractReactionMethod, pars, (R, R_solve, ), cellrange::PB.AbstractCellRange, attribute_name)
    rj = m.reaction

    # VariableReactions corresponding to (R, R_solve)
    R_vars, R_solve_vars = PB.get_variables_tuple(m)
    R_var = only(R_vars)
    R_domvar = R_var.linkvar

    rj.norm_value = PB.get_attribute(R_domvar, :norm_value)

    if pars.const[]
        if attribute_name == :setup
            PB.init_field!(
                only(R), :initial_value, R_domvar, (_, _)->1.0, [], cellrange, (PB.fullname(R_domvar), "", "")
            )
        end
    else
        if attribute_name in (:norm_value, :initial_value)
            if pars.state_norm[]
                R_solve_var = only(R_solve_vars)
                R_solve_domvar = R_solve_var.linkvar
                PB.init_field!(
                    only(R_solve), attribute_name, R_domvar, (_, _)->1/rj.norm_value, [], cellrange, (PB.fullname(R_solve_domvar), " / $(rj.norm_value)", " [from $(PB.fullname(R_domvar))]")
                )
            else
                PB.init_field!(
                    only(R), attribute_name, R_domvar, (_, _)->1.0, [], cellrange, (PB.fullname(R_domvar), "", "")
                )
            end
        end
    end

    return nothing
end

function do_reactionreservoirscalar(m::PB.AbstractReactionMethod, pars, (vars, ), cr::PB.AbstractCellRange, deltat)
    rj = m.reaction

    if pars.state_norm[] && !pars.const[]
        vars.R[] = vars.R_solve[]*rj.norm_value
        vars.R_norm[] = PB.get_total(vars.R_solve[])
    else
        vars.R_norm[]  = PB.get_total(vars.R[])/rj.norm_value
    end

    if hasfield(typeof(vars), :R_delta)
        vars.R_delta[] = PB.get_delta(vars.R[])
    end
    return nothing
end

function do_reactionreservoirscalar_sms(m::PB.AbstractReactionMethod, pars, (vars, ), cr::PB.AbstractCellRange, deltat)
    rj = m.reaction

    vars.R_solve_sms[]  += vars.R_sms[]/rj.norm_value
  
    return nothing
end


"""
    ReactionReservoir, ReactionReservoirTotal

A single (vector) reservoir (state variable) representing a biogeochemical tracer, ie one value per cell in a spatially-resolved
Domain (eg ocean).

State Variables can represent either per-cell concentration, or per-cell moles, set by parameter `state_conc`:
- `state_conc=false` (default): create `R` (mol) and `R_sms` (mol yr-1) as state variable and source-sink,
  calculate `R_conc` (mol m-3)
- `state_conc=true`: create `R_conc` (mol m-3) and `R_conc_sms` (mol m-3 yr-1) as state variable and source-sink,
  calculate `R` (mol), NB: `R_sms` (mol yr-1) is still available and is added to `R_conc_sms`.

In addition:
- if parameter `field_data <: AbstractIsotopeScalar` (eg `IsotopeLinear`), a Property `R_delta` is created.
- `ReactionReservoirTotal` or `ReactionReservoirConcTotal` also calculates the Domain total `R_total` (units mol), eg to check budgets.

Local name prefix `R` should then be renamed using `variable_links:` in the configuration file.

# Initialisation
Initial value is set using `variable_attributes:` in the configuration file, using `R:initial_value` and `R:initial_delta` 
(for 'state_conc=false') or `R_conc:initial_value` and `R_conc:initial_delta` (for a 'state_conc=true').
(NB: even if `initial_value` is set on `R`, it *sets concentration in `mol m-3`*.)

Transport is defined by attributes `:advect`, `:vertical_movement` (m d-1) set on the concentration variable `R_conc`. Optical
extinction is defined by the `:specific_light_extinction` (m^2 mol-1) attribute set on the concentration variable `R_conc`.

# Example configuration in .yaml file
                reservoir_P:  # ocean Phosphorus
                    class: ReactionReservoirTotal       # include _total (mol)
                    parameters:
                        # field_data: ScalarData        # change to IsotopeLinear to represent an isotope                       
                    variable_links:
                        R*: P                           # rename to represent Phosphorus
                    variable_attributes:
                        R:norm_value:           3e-3    # mol m-3, normalisation value (used by some solvers)
                        R:initial_value:        2e-3    # mol m-3, initial concentration

# See also
[`ReactionReservoirWellMixed`](@ref) (one value for the whole Domain), [`ReactionReservoirScalar`](@ref) (one value for
a reservoir in a 0D Domain eg for COPSE [Bergman2004](@cite)), [`ReactionReservoirConst`](@ref) (constant time-independent value
ie no state variable). 

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReservoir{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),

        PB.ParBool("total", false,
            description="true to calculate R_total"),

        PB.ParDouble("limit_delta_conc", 0.0, units="mol m-3",
            description="**EXPERIMENTAL** attempt to limit delta for low/-ve concentrations (0.0 to disable)"),

        PB.ParBool("state_conc", false,
            description="true to define R_conc, R_sms_conc as the state variable pair and calculate R, false to define R, R_sms and calculate R_conc"),
    )
end

abstract type ReactionReservoirTotal <: PB.AbstractReaction end
function PB.create_reaction(::Type{ReactionReservoirTotal}, base::PB.ReactionBase)
    rj = ReactionReservoir(base=base)
    PB.setvalueanddefault!(rj.pars.total, true)
    PB.setfrozen!(rj.pars.total)
    return rj
end


function PB.register_methods!(rj::ReactionReservoir)

    PB.setfrozen!(rj.pars.field_data)
    PB.setfrozen!(rj.pars.total)

    R_attributes=(
        :field_data=>rj.pars.field_data[],
        :calc_total=>rj.pars.total[],
    )
    R_conc_attributes = (
        :field_data=>rj.pars.field_data[], 
        :advect=>true,
        :vertical_movement=>0.0,
        :specific_light_extinction=>0.0,
        :vphase=>PB.VP_Undefined,
        :diffusivity_speciesname=>"",
        :diffusivity=>missing,
        :charge=>missing,
        :gamma=>missing,
    )

    volume  = PB.VarDep("volume",   "m3",       "cell volume (or cell phase volume eg for a sediment with solid and liquid phases)")

    if rj.pars.state_conc[]
        R = PB.VarProp("R",        "mol",      "vector reservoir"; attributes=R_attributes)            
        R_conc = PB.VarStateExplicit("R_conc",   "mol m-3",  "concentration"; attributes=R_conc_attributes)
        PB.add_method_setup_initialvalue_vars_default!(rj, [R_conc])
    else
        R = PB.VarStateExplicit("R",        "mol",      "vector reservoir"; attributes=R_attributes)
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [R],
            convertvars = [volume],
            convertfn = ((volume,), i) -> volume[i],
            convertinfo = " * volume",
        )
        R_conc = PB.VarProp("R_conc",   "mol m-3",  "concentration"; attributes=R_conc_attributes)
    end
    PB.setfrozen!(rj.pars.state_conc)

    do_vars = PB.VariableReaction[R, R_conc, volume,]

    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        push!(do_vars, PB.VarProp("R_delta", "per mil", "isotopic composition"))
        if rj.pars.limit_delta_conc[] > 0.0
            @warn "$(PB.fullname(rj)) using experimental limit_delta_conc $(rj.pars.limit_delta_conc[]) mol m-3"
        end
    end

    PB.add_method_do!(rj, do_reactionreservoir, (PB.VarList_namedtuple(do_vars),))

    if rj.pars.state_conc[]
        R_conc_sms = PB.VarDeriv("R_conc_sms",    "mol m-3 yr-1", "vector reservoir source-sinks",
            attributes=(:field_data=>rj.pars.field_data[], ),
        )
        R_sms = PB.VarTarget("R_sms",    "mol yr-1", "vector reservoir source-sinks",
            attributes=(:field_data=>rj.pars.field_data[], ),
        )
        PB.add_method_do!(rj, do_reactionreservoirconc_sms, (PB.VarList_namedtuple([volume, R_conc_sms, R_sms]),))
    else
        R_sms = PB.VarDeriv("R_sms",    "mol yr-1", "vector reservoir source-sinks",
            attributes=(:field_data=>rj.pars.field_data[], ),
        )
        # sms variable not used by us, but must appear in a method to be linked and created
        PB.add_method_do_nothing!(rj, [R_sms])
    end

    if rj.pars.total[]
        PB.add_method_do_totals_default!(rj, [R])
    end
    PB.setfrozen!(rj.pars.total)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_reactionreservoir(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    if pars.state_conc[]
        @inbounds for i in cellrange.indices
            vars.R[i]  = vars.R_conc[i]*vars.volume[i]
        end
    else
        @inbounds for i in cellrange.indices
            vars.R_conc[i]  = vars.R[i]/vars.volume[i]
        end
    end

    if hasfield(typeof(vars), :R_delta)
        limit_value = pars.limit_delta_conc[]
        if limit_value > 0.0
            # norm_value = PB.get_attribute(rj.var_R, :norm_value)::Float64
            # limit_value = 1e-6  # mol m-3

            @inbounds for i in cellrange.indices
                vars.R_delta[i]  = PB.get_delta_limit(vars.R[i], limit_value*vars.volume[i], 100.0)
            end
        else
            @inbounds for i in cellrange.indices
                vars.R_delta[i]  = PB.get_delta(vars.R[i])
            end
        end
    end

    return nothing
end

function do_reactionreservoirconc_sms(m::PB.AbstractReactionMethod, pars, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        vars.R_conc_sms[i]  += vars.R_sms[i]/vars.volume[i]
    end

    return nothing
end

"""
    ReactionReservoirConst

A single (vector) constant tracer `R_conc` (constant replacement for a [`ReactionReservoir`](@ref)).

Local name prefix `R` should then be renamed using `variable_links:` in the configuration file.

# Initialisation
Set `:initial_value`, `:initial_delta` on `R_conc` (mol m-3) in the `variable_attributes:` section of the
config file.

TODO salinity normalisation.

# Example configuration in .yaml file
                reservoir_B:  # Constant value for ocean Boron 
                    class: ReactionReservoirConst
                    parameters:
                        field_data: IsotopeLinear
                    variable_links:
                        R*: B
                    variable_attributes:                      
                        R_conc:initial_value:       0.4269239 # contemporary value
                        R_conc:initial_delta:       34.0

# See also
[`ReactionReservoir`](@ref)

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_SETUP)
"""
Base.@kwdef mutable struct ReactionReservoirConst{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionReservoirConst)

    var_R_conc = PB.VarProp("R_conc",   "mol m-3",  "concentration",
        attributes=(:field_data=>rj.pars.field_data[], :specific_light_extinction=>0.0,)
    )

    # specify filterfn to initialize var_R_conc even though it isn't a state variables
    PB.add_method_setup_initialvalue_vars_default!(rj, [var_R_conc], filterfn = v->true)

    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        vars = [
            var_R_conc,
            PB.VarProp("R_delta", "per mil", "isotopic composition"),
        ]
        # use do, not setup, so we handle the case where the value is modified after setup
        PB.add_method_do!(rj, do_reactionreservoirconst, (PB.VarList_namedtuple(vars),) )
    else
        # add a dummy method to block any other reaction from also setting a var_R_conc Property Variable
        PB.add_method_do_nothing!(rj, [var_R_conc])  
    end
    PB.setfrozen!(rj.pars.field_data)

    return nothing
end

function do_reactionreservoirconst(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)

    @inbounds for i in cellrange.indices
        vars.R_delta[i]  = PB.get_delta(vars.R_conc[i])
    end
    
    return nothing
end

"""
    ReactionReservoirForced

A single (vector) constant tracer (constant replacement for a [`ReactionReservoir`](@ref)), with forcing.

Calculates `R_conc = R_conc_initial * R_FORCE`.

Local name prefix `R` should then be renamed using `variable_links:` in the configuration file.

# Initialisation
NB: set `:initial_value`, `:initial_delta` on `R_conc_initial` in the `variable_attributes:`
section of the config file.

TODO salinity normalisation.

# See also
[`ReactionReservoirConst`](@ref), [`ReactionReservoir`](@ref)

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReservoirForced{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
    )
end

function PB.register_methods!(rj::ReactionReservoirForced)

    var_R_conc_initial = PB.VarProp("R_conc_initial",   "mol m-3",  "initial concentration",
        attributes=(:field_data=>rj.pars.field_data[], )
    )
    # specify filterfn to initialize var_R_conc_initial even though it isn't a state variable
    PB.add_method_setup_initialvalue_vars_default!(
        rj, [var_R_conc_initial], 
        filterfn = v->true,
    )
   
    do_vars = [
        var_R_conc_initial,
        PB.VarProp("R_conc",   "mol m-3",  "concentration = initial * forcing",
            attributes=(:field_data=>rj.pars.field_data[], :specific_light_extinction=>0.0,)),
        # PB.VarDep(        "volume",   "m3",       "cell volume"),
        PB.VarDepScalar(  "R_FORCE",  "",         "forcing factor"),
    ]
    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        push!(do_vars, PB.VarProp("R_delta", "per mil", "isotopic composition"))
    end
    PB.setfrozen!(rj.pars.field_data)

    PB.add_method_do!(
        rj,
        do_reactionreservoirconstforced,
        (PB.VarList_namedtuple(do_vars),),
    )

    return nothing
end


function do_reactionreservoirconstforced(m::PB.ReactionMethod, (do_vardata, ), cellrange::PB.AbstractCellRange, deltat)

    force_fac = do_vardata.R_FORCE[]
    @inbounds for i in cellrange.indices
        do_vardata.R_conc[i]  = do_vardata.R_conc_initial[i]*force_fac
   
        if hasfield(typeof(do_vardata), :R_delta)
            do_vardata.R_delta[i]  = PB.get_delta(do_vardata.R_conc_initial[i])
        end
    end

    return nothing
end


"""
    ReactionReservoirWellMixed

A Scalar Reservoir representing a well-mixed tracer in a vector Domain (eg ocean). Provides a scalar state Variable `R` and `R_sms`,
and also vector Variables: `R_conc` (set to uniform concentration), and `R_vec_sms` Target for accumulating fluxes that is added to the scalar
`R_sms`. 

# Initialisation
Initial value is set in the `variable_attributes:` section in the configuration file, using `R:initial_value` and `R:initial_delta`.
NB: may be initialized either from mean concentration or total (set by parameter `initialization_type`)

TODO salinity normalisation.

# See also
[`ReactionReservoir`](@ref) (one value per cell), [`ReactionReservoirScalar`](@ref) (one value for
a reservoir in a 0D Domain eg for COPSE [Bergman2004](@cite)).

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionReservoirWellMixed{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),

        PB.ParString("initialization_type", "conc",
            allowed_values=["conc", "total"],
            description=":initial_value attribute represents conc (mol m-3) or total (mol)"),
    )

    total_norm::Float64 = NaN
end

function PB.register_methods!(rj::ReactionReservoirWellMixed)

    R = PB.VarStateExplicitScalar(  "R",        "mol",      "scalar reservoir",
        attributes=(:field_data=>rj.pars.field_data[],))
    volume_total = PB.VarDepScalar("volume_total", "m^3", "total volume")
    vars = [
        R,
        PB.VarPropScalar(           "R_norm",   "",         "scalar reservoir normalized"),
        PB.VarProp(                 "R_conc",   "mol m-3",  "concentration",           
            attributes=(:field_data=>rj.pars.field_data[], :specific_light_extinction=>0.0,)),
        volume_total,
    ]

    if rj.pars.field_data[] <: PB.AbstractIsotopeScalar
        push!(vars, PB.VarProp("R_delta", "per mil", "isotopic composition"))
    end
    PB.setfrozen!(rj.pars.field_data)

    # callback function to store Variable norm during setup
    function setup_callback(m, attribute_value, v, vdata)
        v.localname == "R" || error("setup_callback unexpected Variable $(PB.fullname(v))")
        if attribute_value == :norm_value
            m.reaction.total_norm = PB.value_ad(PB.get_total(vdata[]))
        end
        return nothing
    end

    if rj.pars.initialization_type[] == "total"
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [R],
            setup_callback=setup_callback
        )
    elseif rj.pars.initialization_type[] == "conc"
        PB.add_method_setup_initialvalue_vars_default!(
            rj, [R],
            convertvars=[volume_total],
            convertfn=((volume_total,), i) -> volume_total[],
            convertinfo=" * total_volume",
            setup_callback=setup_callback,
        )
    end
    PB.setfrozen!(rj.pars.initialization_type)

    threadsafe = get(rj.external_parameters, "threadsafe", false)

    PB.add_method_do!(rj, do_reservoir_well_mixed, (PB.VarList_namedtuple(vars),))

    # add method to sum per-cell vec_sms to scalar sms
    vars_sms = [
        PB.VarDerivScalar(          "R_sms",    "mol yr-1", "scalar reservoir source-sinks",
            attributes=(:field_data=>rj.pars.field_data[],)),
        PB.VarTarget(               "R_vec_sms","mol yr-1", "vector reservoir source-sinks",
            attributes=(:field_data=>rj.pars.field_data[],))
    ]
    PB.add_method_do!(rj, do_reservoir_well_mixed_sms, (PB.VarList_namedtuple(vars_sms),); p=threadsafe)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end


function do_reservoir_well_mixed(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    rj = m.reaction

    if Threads.threadid() == 1 # only set once ! TODO implicitly assuming will always be called from threadid==1
        vars.R_norm[] = PB.get_total(vars.R[])/rj.total_norm
    end

    @inbounds for i in cellrange.indices
        vars.R_conc[i]  = vars.R[]/vars.volume_total[]
    end

    if hasfield(typeof(vars), :R_delta)
        delta = PB.get_delta(vars.R[])
        @inbounds for i in cellrange.indices
            vars.R_delta[i] = delta
        end
    end

    return nothing
end


function do_reservoir_well_mixed_sms(
    m::PB.ReactionMethod,
    (vars, ),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    threadsafe = m.p
    
    # accumulate first into a subtotal,
    # in order to minimise time spent with lock held if this is one tile of a threaded model
    subtotal = zero(vars.R_sms[])
    @inbounds for i in cellrange.indices
        subtotal += vars.R_vec_sms[i]
    end
    if threadsafe
        PB.atomic_add!(vars.R_sms, subtotal)
    else
        vars.R_sms[] += subtotal
    end

    return nothing
end


"""
    ReactionConst, ReactionScalarConst

Create constant Property Variables with names from parameter `constnames`. 

# Initialisation
Constant values set by `:initial_value`, `:initial_delta` attributes in the `variable_attributes:` section of the configuration file.

# Example configuration in .yaml file
        atmfloor:
            reactions:                
                floorstubgasmr:  # Provide mixing-ratio boundary condition for a subset of atmospheric variables
                    class: ReactionConst
                    parameters:
                        constnames: ["O2_mr", "CH4_mr", "CO2_mr", "H2_mr"] #, "H2O_mr"]  
                    variable_attributes:        
                        O2_mr:initial_value:  [0.21]  # mol/mol
                        CH4_mr:initial_value:  [0.7443e-6] # [.NaN] # [1.271e-6] # [1.8e-6]  # mol/mol
                        CO2_mr:initial_value:  [300e-6]  # mol/mol
                        H2_mr:initial_value:   [.NaN]  # mol/mol          

# See also
[`ReactionReservoirConst`](@ref), [`ReactionReservoirScalar`](@ref).  These provide additional variables (eg `R_delta`) to 
allow them to function as a drop-in replacement for a non-constant Reservoir.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_SETUP)
"""
Base.@kwdef mutable struct ReactionConst{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("constnames", ["constvar"],
            description="vector of names for constant Variables. Isotopes use <name>::CIsotope syntax"),
    )

    scalar::Bool = false
end

function PB.register_methods!(rj::ReactionConst)

    vars_const = []

    for varnameisotope in rj.pars.constnames
        varname, IsotopeType = PB.split_nameisotope(varnameisotope, rj.external_parameters)
        if rj.scalar
            constvar = PB.VarPropScalarStateIndep(varname, "unknown", "constant value", attributes=(:field_data=>IsotopeType, ))
        else
            constvar = PB.VarPropStateIndep(varname, "unknown", "constant value", attributes=(:field_data=>IsotopeType, ))
        end
        push!(vars_const, constvar)
    end

    # add a dummy method to block any other reaction from also creating (and modifying!) a Property Variable
    PB.add_method_do_nothing!(rj, vars_const)  

    # specify filterfn to initialize vars_const even though they aren't state variables
    PB.add_method_setup_initialvalue_vars_default!(rj, vars_const, filterfn = v->true)

    return nothing
end

abstract type ReactionScalarConst <: PB.AbstractReaction end
function PB.create_reaction(::Type{ReactionScalarConst}, base::PB.ReactionBase)
    rj = ReactionConst(base=base, scalar=true)
    return rj
end


"""
    ReservoirLinksVector(isotope_data::Dict, reservoirlist) -> (res::Vector, sms::Vector, diag::Vector)

Convenience function to create variables required for a Reaction to link to a list of Reservoir variables.
`res` contains VariableReactions `<reservoir_name>`, `sms` `<reservoir_name>_sms` that link to `State`, `State_sms` variables.
`diag` contains VariableReactions `<reservoir_name>_norm` etc that link to additional properties.

# Arguments
 - `reservoirlist::[(reservoir_name[::Isotope],       units,   description), ...]`: list of Reservoirs
"""
function ReservoirLinksVector(isotope_data::Dict, reservoirlist)
    # create a list of required variables
    reservoirs = Vector{PB.VariableReaction}()
    sms = Vector{PB.VariableReaction}()
    diagnostics = Vector{PB.VariableReaction}()
    for reservoir in reservoirlist
        nameisotopebr, units, description  = reservoir

        (var_optional, nameisotope) = PB.strip_brackets(nameisotopebr)

        name_root, IsotopeType = PB.split_nameisotope(nameisotope, isotope_data)

        if var_optional
            lb="("; rb=")"
        else
            lb=""; rb=""
        end

        push!(
            reservoirs,
            PB.VarDepScalar(lb*name_root*rb, units, description,
                attributes=(:field_data=>IsotopeType,))
        )
        push!(
            sms,
            PB.VarContribScalar(lb*name_root*"_sms"*rb, units*" yr-1", description*" source-sinks",
                attributes=(:field_data=>IsotopeType,))
        )
        push!(
            diagnostics,
            PB.VarDepScalar(lb*name_root*"_norm"*rb, "", "normalized "*description)
        )
        if IsotopeType <: PB.AbstractIsotopeScalar
            push!(
                diagnostics,
                PB.VarDepScalar(lb*name_root*"_delta"*rb, "per mil", "isotope delta "*description)
            )
        end
    end

    return (reservoirs, sms, diagnostics)

end

end # module
