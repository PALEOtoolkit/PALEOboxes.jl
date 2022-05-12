"""
    add_method_do_totals_default!(react::AbstractReaction, variables=PB.get_variables(react))

Create and add a default method to initialize Variables to zero at beginning of each timestep. Defaults to
adding all Variables from `react` with `:initialize_to_zero` attribute `true`.

NB: TODO variables are converted to VarDep (no dependency checking or sorting needed, could define a VarInit or similar?)
"""
function add_method_do_totals_default!(
    reaction::AbstractReaction, total_candidates=get_variables(reaction);
    kwargs...
)
        
    return add_method_do!(
        reaction, 
        create_totals_method(reaction, total_candidates; kwargs...),
    )
end

"""
    create_totals_method(reaction, var; [methodname] [, total_localname] [, operatorID]) -> ReactionMethod 

Create a method to add a total variable (Scalar Property) to `reaction` for `var`. If `total_localname` is not supplied, it will be generated
by appending `_total` to `var.localname`.

NB: total Variables will require initialization to zero using [`add_method_initialize_zero_vars_default!`](@ref)
"""
function create_totals_method(
    reaction::AbstractReaction, var::VariableReaction;
    methodname="totals",
    total_localname=nothing, 
    operatorID=reaction.operatorID
)

    return _create_totals_method(
        reaction, 
        methodname,
        ([var], [create_totalvar(var, total_localname)]),
        operatorID
    )
end

"""
    create_totals_method(reaction, total_candidates; 
        [filterfn] [, methodname] [, total_localnames] [, operatorID]) -> ReactionMethod 

Create a method to add total variables (Scalar Properties), for Variables in `total_candidates` that match `filterfn`
(defaults to those that are Array Variables and have attribute ``:calc_total == true`).

NB: total Variables will require initialization to zero using [`add_method_initialize_zero_vars_default!`](@ref)
"""
function create_totals_method(
    reaction::AbstractReaction, total_candidates;
    filterfn=v->get_attribute(v, :calc_total, false),
    methodname="totals",
    total_localnames=nothing, 
    operatorID=reaction.operatorID
)
     
    if isnothing(total_localnames)
        # remove duplicates (eg from setup and do methods)
        vars_all = Dict(v.localname =>v for v in total_candidates if filterfn(v))        
        vars = [v for (k, v) in vars_all]
    else
        # check no filtering if total_localnames supplied
        vars = [v for v in total_candidates if filterfn(v)]
        length(vars) == length(total_candidates) ||
            error("create_totals_method: configuration error, Reaction $(fullname(reaction)), filterfn not supported with total_localnames")
        length(total_localnames) == length(vars) ||   
            error("create_totals_method: configuration error, Reaction $(fullname(reaction)) length(total_localnames) != number of total Variables")
    end

    var_totals = []
    for (idx, v) in enumerate(vars)
        if isnothing(total_localnames)
            total_localname = nothing
        else
            total_localname = total_localnames[idx]
        end
        push!(var_totals, create_totalvar(v, total_localname))
    end

    return _create_totals_method(
        reaction, 
        methodname,
        (vars, var_totals),
        operatorID
    )     
end

function create_totalvar(var, total_localname=nothing)
 
    if isnothing(total_localname)
        total_localname = var.localname*"_total"
    end
    return VarPropScalar(total_localname, get_attribute(var, :units), "total "*get_attribute(var, :description),
        attributes=(:field_data=>get_attribute(var, :field_data), :initialize_to_zero=>true, :atomic=>true),
    )
end

"""
    do_vartotals(varstats, varstats_data, cellrange::AbstractCellRange)

Calculate totals. `varstats_data` is `accessors` returned by `create_accessors`. Call from the 
`ReactionPhase` supplied to [`add_var`](@ref).
"""
function do_vartotals(m::AbstractReactionMethod, (vars_data, var_totals_data), cellrange::AbstractCellRange, deltat)

    function calc_total(data, total_data, cellrange)
        # accumulate first into a subtotal, and then into total_data[], in order to minimise time spent with lock held
        # if this is one tile of a threaded model
        subtotal = zero(total_data[])
        @inbounds for i in cellrange.indices
            subtotal += data[i]
        end
        atomic_add!(total_data, subtotal)
        return nothing
    end

    length(vars_data) == length(var_totals_data) || 
        error("do_vartotals: components length mismatch $(fullname(m)) $(get_variables_tuple(m)) (check :field_data (ScalarData, IsotopeLinear etc) match)")
    #  if cellrange.operatorID == 0 || cellrange.operatorID in totals_method.operatorID
    for iv in eachindex(vars_data)
        calc_total(vars_data[iv], var_totals_data[iv], cellrange)
    end
    # end
    
    return nothing
end

function _create_totals_method(
    reaction::AbstractReaction, 
    methodname::AbstractString,
    (vars, var_totals),
    operatorID,
) 

    return ReactionMethod(
        do_vartotals,
        reaction,
        methodname,
        (VarList_components(VarDep.(vars)), VarList_components(var_totals)),
        nothing,
        operatorID,
        reaction.domain
    )
    
end
