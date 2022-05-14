# -*- coding: utf-8 -*-
module VariableStats

import PALEOboxes as PB
using ..DocStrings

import Infiltrator # Julia debugger


"""
    ReactionSum, ReactionVectorSum

A sum of variables (eg budget).
- If Parameter `component_to_add == 0`, all components of Isotopes are included.
- If Parameter `component_to_add == component_number`, a single component only is included.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionSum{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec( "vars_to_add", ["2*myvar", "myothervar", "-1*mythirdvar"],
            description="vector of variable names to add, eg [2*myvar, myothervar, -1*mythirdvar]"),
        PB.ParString("vars_prefix", "",
            description="optional prefix for vars_to_add"),

        PB.ParInt("component_to_add", 0,
            description="component to add, 0 for all"),

        PB.ParBool("vectorsum", false,
            description="true to accumulate sum into vector Variable, "*
                "false to accumulate sum into scalar (adding vector cells if necessary)"),
    )
end

abstract type ReactionVectorSum <: PB.AbstractReaction end
function PB.create_reaction(::Type{ReactionVectorSum}, base::PB.ReactionBase)
    rj = ReactionSum(base=base)
    PB.setvalueanddefault!(rj.pars.vectorsum, true)
    PB.setfrozen!(rj.pars.vectorsum)
    return rj
end

function PB.register_methods!(rj::ReactionSum)

    vars_to_add = []
    var_multipliers = Float64[]

    for varmultname in rj.pars.vars_to_add.v
        # parse multiplier
        svmn = split(varmultname, ['*', ' '], keepempty=false)
        if length(svmn) == 1
            mult, varname = (1.0, rj.pars.vars_prefix.v*svmn[1])
        elseif length(svmn) == 2
            mult, varname = (parse(Float64, svmn[1]), rj.pars.vars_prefix.v*svmn[2])
        else
            error("reaction ", fullname(rj), "invalid field in vars_to_add ", varmultname)
        end
        @info "reaction $(PB.fullname(rj)) add $mult * $varname"
        push!(var_multipliers, mult)

        # generate new name with domain prefix to disambiguate eg O2 in atm and ocean
        (linkreq_domain, linkreq_subdomain, linkreq_name, link_optional) =
            PB.parse_variablereaction_namestr(varname)
        localname = PB.combine_link_name(linkreq_domain, "", linkreq_name, sep="_")
         # mark all vars_to_add as optional to help diagnose config errors
         # default :field_data=>PB.UndefinedData  to allow Variable to link to any data type (this is checked later)
        push!(vars_to_add, 
            PB.VarDep(localname, "", "", link_namestr="("*varname*")")
        )
    end

    if rj.pars.vectorsum.v
        methodfn = do_vectorsum
        var_sum = PB.VarProp("sum", "", "sum of specified variables")
    else
        methodfn = do_scalarsum
        var_sum = PB.VarPropScalar("sum", "", "sum of specified variables")
    end

    PB.add_method_do!(
        rj,
        methodfn,
        (
            PB.VarList_single(var_sum, components=true),
            PB.VarList_tuple(vars_to_add, components=true)
        ),
        p = (Tuple(var_multipliers), 0:-1) # comprange added later
    )

    return nothing
end

function PB.register_dynamic_methods!(rj::ReactionSum)

    # update method now Variable are linked hence components known
    method_sum = PB.get_method_do(rj, rj.pars.vectorsum.v ? "do_vectorsum" : "do_scalarsum")

    var_sum = PB.get_variable(method_sum, "sum")
    vars_to_add = PB.get_variables(method_sum, filterfn = v->v.localname != "sum")

    # check variable components match and update var_sum.components
    if rj.pars.component_to_add.v == 0
        # add all components of vars_to_add Variables
        # check all Variable have the same data
        firstvar_d = nothing
        for v in vars_to_add
            !isnothing(v.linkvar) ||
                error("Reaction $(PB.fullname(rj)) variable $(v.localname) is not linked: check configuration")            
            linkvar_d = PB.get_attribute(v.linkvar, :field_data)
            if v === first(vars_to_add)
                firstvar_d = linkvar_d
                PB.set_attribute!(var_sum, :field_data,  firstvar_d)
                @info "Reaction $(PB.fullname(rj)) Variable $(PB.fullname(var_sum.linkvar)) adding data=$firstvar_d"
            end
            linkvar_d == firstvar_d ||
                error("Reaction $(PB.fullname(rj)) not all variables to be summed have the same :field_data Type: $(PB.fullname(v.linkvar)) $(linkvar_d) != $(PB.fullname(first(vars_to_add).linkvar)) $(firstvar_d)")
        end        
    else
        # add first component of vars_to_add Variables
        PB.set_attribute!(var_sum, :field_data,  PB.ScalarData)
        @info "Reaction $(PB.fullname(rj)) Variable $(PB.fullname(var_sum.linkvar)) adding single component $(rj.pars.component_to_add.v)"
    end

    if rj.pars.component_to_add.v == 0
        comprange = 1:PB.num_components(PB.get_attribute(var_sum, :field_data))
    else
        comprange = rj.pars.component_to_add.v:rj.pars.component_to_add.v
    end

    # update method_sum
    method_sum.p = (method_sum.p[1], comprange)
end


function do_scalarsum(
    m::PB.ReactionMethod,
    (var_sum_data, vars_to_add_data),
    cellrange::PB.AbstractCellRange,
    deltat
)
    (multipliers, comprange) = m.p

    function _add_var( multiplier, var_to_add)
        for j in comprange
            var_sum_data[j][] += multiplier * sum(var_to_add[j])
        end
        return nothing
    end

    for j in comprange
        var_sum_data[j][] = 0.0
    end

    PB.IteratorUtils.foreach_longtuple(_add_var, multipliers, vars_to_add_data)

    return nothing
end


function do_vectorsum(
    m::PB.ReactionMethod,
    (var_sum_data, vars_to_add_data),
    cellrange::PB.AbstractCellRange,
    deltat
)
    (multipliers, comprange) = m.p

    function _add_var(multiplier, var_to_add)
        @inbounds for j in comprange
            for i in cellrange.indices
                var_sum_data[j][i] += multiplier * var_to_add[j][i]
            end
        end
        return nothing
    end

    for j in comprange
        @inbounds for i in cellrange.indices
            var_sum_data[j][i] = 0.0
        end
    end

    PB.IteratorUtils.foreach_longtuple(_add_var, multipliers, vars_to_add_data)

    return nothing
end


"""
    ReactionWeightedMean

Weighted mean (eg by area or volume) of Variable

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionWeightedMean{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
    )


end


function PB.register_methods!(rj::ReactionWeightedMean)

    vars = [
        PB.VarDep(       "var", "measure-1", "variable to calculate weighted mean from"),
        PB.VarDep(       "measure", "", "cell area or volume"),
        PB.VarDepScalar( "measure_total", "", "total Domain area or volume"),
        PB.VarPropScalar("var_mean", "", "weighted mean over Domain area or volume",
            attributes=(:initialize_to_zero=>true, :atomic=>true)),
    ]

    PB.add_method_do!(rj, do_weighted_mean, (PB.VarList_namedtuple(vars),))

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function do_weighted_mean(m::PB.ReactionMethod, (vars,), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # use a subtotal to minimise time lock held if this is one tile of a threaded model
    subtotal = zero(vars.var_mean[])
    @inbounds for i in cellrange.indices
        subtotal += vars.var[i]*vars.measure[i] # normalisation by measure_total below
    end

    PB.atomic_add!(vars.var_mean, subtotal/vars.measure_total[])
    
    return nothing
end


"""
    ReactionAreaVolumeValInRange

Fraction of Domain area or volume with Variable in a range of values.

# Parameters
$(PARS)

# Methods and Variables
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionAreaVolumeValInRange{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble( "range_min", -Inf, units="mol m-3",
            description="minimum value to include in frac"),
        PB.ParDouble( "range_max", Inf, units="mol m-3",
            description="maximum value to include in frac"),
    )


end


function PB.register_methods!(rj::ReactionAreaVolumeValInRange)

    vars = [
        PB.VarDep(       "rangevar", "mol m-3", "variable to check within range"),
        PB.VarDep(       "measure", "", "cell area or volume"),
        PB.VarDepScalar( "measure_total", "", "total Domain area or volume"),
        PB.VarPropScalar("frac", "", "fraction of Domain area or volume in specified range",
            attributes=(:initialize_to_zero=>true, :atomic=>true)),
    ]

    PB.add_method_do!(rj, do_area_volume_in_range, (PB.VarList_namedtuple(vars),))

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function do_area_volume_in_range(m::PB.ReactionMethod, (vars,), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction

    # use a subtotal to minimise time lock held if this is one tile of a threaded model
    frac_subtotal = zero(vars.frac[])
    @inbounds for i in cellrange.indices
        if vars.rangevar[i] >= rj.pars.range_min.v && vars.rangevar[i] <= rj.pars.range_max.v
            frac_subtotal += vars.measure[i] # normalisation by measure_total below
        end
    end

    PB.atomic_add!(vars.frac, frac_subtotal/vars.measure_total[])

    return nothing
end





end # module
