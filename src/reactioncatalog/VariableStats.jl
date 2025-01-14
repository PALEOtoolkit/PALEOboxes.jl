# -*- coding: utf-8 -*-
module VariableStats

import ...PALEOboxes as PB
using ...PALEOboxes: @public
using ..DocStrings

@public ReactionSum, ReactionVectorSum, ReactionWeightedMean, ReactionAreaVolumeValInRange

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

    var_multipliers::Vector{Float64} = Float64[]
    comprange::UnitRange{Int64} = 0:-1
end

abstract type ReactionVectorSum <: PB.AbstractReaction end
function PB.create_reaction(::Type{ReactionVectorSum}, base::PB.ReactionBase)
    rj = ReactionSum(base=base)
    PB.setvalueanddefault!(rj.pars.vectorsum, true)
    PB.setfrozen!(rj.pars.vectorsum)
    return rj
end

function PB.register_methods!(rj::ReactionSum)

    io = IOBuffer()
    println(io, "register_methods: $(PB.fullname(rj)) $(PB.typename(rj))")

    vars_to_add = []
    empty!(rj.var_multipliers)

    for varmultname in rj.pars.vars_to_add
        mult, varname = PB.parse_number_name(varmultname; io, errmsg="reaction $(PB.fullname(rj)) invalid field in vars_to_add ")
        varname = rj.pars.vars_prefix[]*varname
        println(io, "    add $mult * $varname")
        push!(rj.var_multipliers, mult)

        # generate new name with domain prefix to disambiguate eg O2 in atm and ocean
        (linkreq_domain, linkreq_subdomain, linkreq_name, link_optional) =
            PB.parse_variablereaction_namestr(varname)
        localname = PB.combine_link_name(linkreq_domain, "", linkreq_name, sep="_")

        # mark all vars_to_add as optional to help diagnose config errors
        # default :field_data=>PB.UndefinedData  to allow Variable to link to any data type (this is checked later)
        v =  PB.VarDep(localname => "("*varname*")", "unknown", "")

        # no length check if scalar sum
        if !rj.pars.vectorsum[]
            PB.set_attribute!(v, :check_length, false; allow_create=true)
        end
        push!(vars_to_add, v)
    end

    if rj.pars.vectorsum[]
        methodfn = do_vectorsum
        var_sum = PB.VarProp("sum", "unknown", "sum of specified variables")
    else
        methodfn = do_scalarsum
        var_sum = PB.VarPropScalar("sum", "unknown", "sum of specified variables")
    end

    PB.add_method_do!(
        rj,
        methodfn,
        (
            PB.VarList_single(var_sum, components=true),
            PB.VarList_tuple(vars_to_add, components=true)
        ),
    )

    @info String(take!(io))

    return nothing
end

function PB.register_dynamic_methods!(rj::ReactionSum)

    # update method now Variable are linked hence components known
    method_sum = PB.get_method_do(rj, rj.pars.vectorsum[] ? "do_vectorsum" : "do_scalarsum")

    var_sum = PB.get_variable(method_sum, "sum")
    vars_to_add = PB.get_variables(method_sum, filterfn = v->v.localname != "sum")

    # set units from sum variable
    # (may have been set explicitly in yaml file)
    sum_units = PB.get_attribute(var_sum, :units)
    for v in vars_to_add
        PB.set_attribute!(v, :units,  sum_units)
    end

    # check variable components match and update var_sum.components
    if rj.pars.component_to_add[] == 0
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
        @info "Reaction $(PB.fullname(rj)) Variable $(PB.fullname(var_sum.linkvar)) adding single component $(rj.pars.component_to_add[])"
    end

    if rj.pars.component_to_add[] == 0
        rj.comprange = 1:PB.num_components(PB.get_attribute(var_sum, :field_data))
    else
        rj.comprange = rj.pars.component_to_add[]:rj.pars.component_to_add[]
    end

end


function do_scalarsum(
    m::PB.ReactionMethod,
    (var_sum_data, vars_to_add_data),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    comprange = rj.comprange

    function _add_var(var_to_add, multiplier)
        for j in comprange
            var_sum_data[j][] += multiplier * sum(var_to_add[j])
        end
        return nothing
    end

    for j in comprange
        var_sum_data[j][] = 0.0
    end

    PB.IteratorUtils.foreach_longtuple(_add_var, vars_to_add_data, rj.var_multipliers)

    return nothing
end


function do_vectorsum(
    m::PB.ReactionMethod,
    (var_sum_data, vars_to_add_data),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction
    comprange =  rj.comprange

    function _add_var(var_to_add, multiplier)
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

    PB.IteratorUtils.foreach_longtuple(_add_var, vars_to_add_data, rj.var_multipliers)

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
        PB.ParType(PB.AbstractData, "field_data", PB.ScalarData,
            allowed_values=PB.IsotopeTypes,
            description="disable / enable isotopes and specify isotope type"),
    )


end


function PB.register_methods!(rj::ReactionWeightedMean)

    vars = [
        PB.VarDep(       "var", "unknown", "variable to calculate weighted mean from",
            attributes=(:field_data=>rj.pars.field_data[],)),
        PB.VarDep(       "measure", "unknown", "cell area or volume"),
        PB.VarDepScalar( "measure_total", "unknown", "total Domain area or volume"),
        PB.VarPropScalar("var_mean", "unknown", "weighted mean over Domain area or volume",
            attributes=(:field_data=>rj.pars.field_data[], :initialize_to_zero=>true,)),
    ]

    threadsafe = get(rj.external_parameters, "threadsafe", false)

    PB.add_method_do!(rj, do_weighted_mean, (PB.VarList_namedtuple(vars),); p=threadsafe)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function do_weighted_mean(m::PB.ReactionMethod, (vars,), cellrange::PB.AbstractCellRange, deltat)
    threadsafe = m.p

    # use a subtotal to minimise time lock held if this is one tile of a threaded model
    subtotal = zero(vars.var_mean[])
    @inbounds for i in cellrange.indices
        subtotal += vars.var[i]*vars.measure[i] # normalisation by measure_total below
    end

    var_mean_to_add = subtotal/vars.measure_total[]

    if threadsafe
        PB.atomic_add!(vars.var_mean, var_mean_to_add)
    else
        vars.var_mean[] += var_mean_to_add
    end

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
        PB.VarDep(       "rangevar", "unknown", "variable to check within range";
            attributes=(:field_data=>PB.ScalarData, )), # total only, not isotope
        PB.VarDep(       "measure", "unknown", "cell area or volume"),
        PB.VarDepScalar( "measure_total", "unknown", "total Domain area or volume"),
        PB.VarPropScalar("frac", "", "fraction of Domain area or volume in specified range",
            attributes=(:initialize_to_zero=>true,)),
    ]

    threadsafe = get(rj.external_parameters, "threadsafe", false)
    PB.add_method_do!(rj, do_area_volume_in_range, (PB.VarList_namedtuple(vars),); p=threadsafe)

    PB.add_method_initialize_zero_vars_default!(rj)

    return nothing
end

function do_area_volume_in_range(m::PB.ReactionMethod, pars, (vars,), cellrange::PB.AbstractCellRange, deltat)
    threadsafe = m.p

    # use a subtotal to minimise time lock held if this is one tile of a threaded model
    frac_subtotal = zero(vars.frac[])
    @inbounds for i in cellrange.indices
        if vars.rangevar[i] >= pars.range_min[] && vars.rangevar[i] <= pars.range_max[]
            frac_subtotal += vars.measure[i] # normalisation by measure_total below
        end
    end

    frac_to_add = frac_subtotal/vars.measure_total[]

    if threadsafe
        PB.atomic_add!(vars.frac, frac_to_add)
    else
        vars.frac[] += frac_to_add
    end

    return nothing
end





end # module
