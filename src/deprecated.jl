
"""
    get_statevar

DEPRECATED - moved to PALEOmodel
"""
function get_statevar end

"""
    get_statevar_norm

DEPRECATED - moved to PALEOmodel
"""
function get_statevar_norm end

"DEPRECATED
Optional: sanitize `values` for storing as model output.
Default implementation is usually OK - only implement for custom types that should be converted to standard types for storage"
get_values_output(values, data_type::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space, mesh) = values

"DEPRECATED sanitized version of `values`, suitable for storing as output"
function get_values_output(field::Field{FieldData, Space, V, N, Mesh}) where {FieldData, Space, V, N, Mesh}
    return get_values_output(field.values, FieldData, field.data_dims, Space, field.mesh)
end


"""
    get_data_output(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int) -> get_values_output(field.values)
    get_data_output(var::VariableDomain, domaindata::AbstractDomainData) -> get_values_output(field.values)

DEPRECATED
Get a sanitized version of Variable `var` data array for storing as output
from [`get_values_output`]@ref)`(`[`Field`](@ref).values`)`
"""
get_data_output(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int=1) = get_values_output(get_field(var, modeldata, arrays_idx))
get_data_output(var::VariableDomain, domaindata::AbstractDomainData) = get_values_output(get_field(var, domaindata))


function add_reaction_factory(ReactionType::Type{<:AbstractReaction})
    Base.depwarn("call to deprecated add_reaction_factory($ReactionType), this does nothing and can be removed", :add_reaction_factory, force=true)
end

"""
    add_par(reaction::AbstractReaction, par::AbstractParameter)
    add_par(reaction::AbstractReaction, objectwithpars)

Add a single parameter or parameters from fields of `objectwithpars` to a new Reaction.

Not usually needed: Parameters in `pars::ParametersTuple`` will be added automatically, only needed if there are additional
Parameters that are not members of `pars`.
"""
function add_par(@nospecialize(reaction::AbstractReaction), par::AbstractParameter)
    
    error("add_par has been removed: all Parameters should be added to field pars::PB.ParametersTuple, please update your code "*
        "reaction $reaction par $par")

    return nothing
end

function add_par(@nospecialize(reaction::AbstractReaction), objectwithpars)
    error("add_par has been removed: all Parameters should be added to field pars::PB.ParametersTuple, please update your code "*
        "reaction $reaction")
end
