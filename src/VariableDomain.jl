"""
    VariableDomPropDep

[`Model`](@ref) ([`Domain`](@ref)) `VariableDomain` linking a single `VariableReaction{VT_ReactProperty}`
to multiple `VariableReaction{VT_ReactDependency}`.
"""
Base.@kwdef mutable struct VariableDomPropDep <: VariableDomain
    ID::Int64
    domain::AbstractDomain
    name::String
    master::VariableReaction

    var_property::Vector{VariableReaction{VT_ReactProperty}}  =
        Vector{VariableReaction{VT_ReactProperty}}() # 0, 1 (or 2 if linked by both setup and do methods)
    var_dependencies::Vector{VariableReaction{VT_ReactDependency}} =
        Vector{VariableReaction{VT_ReactDependency}}()

end

get_var_type(var::VariableDomPropDep) = VT_DomPropDep

"""
    VariableDomContribTarget

[`Model`](@ref) ([`Domain`](@ref)) `VariableDomain` linking a single `VariableReaction{VT_ReactTarget}`
to multiple `VariableReaction{VT_ReactContributor}` and `VariableReaction{VT_ReactDependency}`.
"""
Base.@kwdef mutable struct VariableDomContribTarget <: VariableDomain
    ID::Int64
    domain::AbstractDomain
    name::String
    master::VariableReaction

    var_target::Union{Nothing, VariableReaction{VT_ReactTarget}} = 
        nothing
    var_contributors::Vector{VariableReaction{VT_ReactContributor}} =
        Vector{VariableReaction{VT_ReactContributor}}()
    var_dependencies::Vector{VariableReaction{VT_ReactDependency}} =
        Vector{VariableReaction{VT_ReactDependency}}()

 end

 get_var_type(var::VariableDomContribTarget) = VT_DomContribTarget

####################################
# Properties
###################################

"""
    Base.getproperty(var::VariableDomain, s::Symbol)

Define [`VariableDomain`](@ref) properties.

`var.attributes` is forwarded to the linked [`VariableReaction`](@ref) defined by the `var.master` field.
"""
function Base.getproperty(var::VariableDomain, s::Symbol)
    if s == :attributes
        return var.master.attributes
    else
        return getfield(var, s)
    end
end

is_scalar(var::VariableDomain) = (var.attributes[:space] === ScalarSpace)

"true if data array assigned"
function is_allocated(var::VariableDomain, modeldata::AbstractModelData)
    return !isnothing(get_field(var, modeldata))    
end

"fully qualified name"
function fullname(var::VariableDomain)
    return var.domain.name*"."*var.name
end

function host_dependent(var::VariableDomPropDep)
    return isempty(var.var_property)
end

function host_dependent(var::VariableDomContribTarget)
    return isnothing(var.var_target)
end

################################################
# Field and data access
##################################################

"""
    set_field!(var::VariableDomain, modeldata, field::Field)
 
Set `VariableDomain` data to `field`
"""
function set_field!(var::VariableDomain, modeldata::AbstractModelData, field::Field)

    variable_data = get_domaindata(modeldata, var.domain).variable_data

    variable_data[var.ID] = field

    return nothing
end

"""
    set_data!(var::VariableDomain, modeldata, data)
 
Set [`VariableDomain`](@ref) to a Field containing `data`.

Calls [`wrap_field`](@ref) to create a new [`Field`](@ref).
"""
function set_data!(var::VariableDomain, modeldata::AbstractModelData, data)

    variable_data = get_domaindata(modeldata, var.domain).variable_data

    variable_data[var.ID] =  wrap_field(
        data,
        get_attribute(var, :field_data),
        Tuple(get_data_dimension(domain, dn) for dn in get_attribute(var, :data_dims)),
        get_attribute(var, :datatype),
        get_attribute(var, :space), 
        var.domain.grid,
    )

    return nothing
end


"""
    get_field(var::VariableDomain, modeldata::AbstractModelData) -> field::Field
    get_field(var::VariableDomain, domaindata::AbstractDomainData) -> field::Field

Get Variable `var` `field`
"""
function get_field(var::VariableDomain, modeldata::AbstractModelData)
    domaindata = get_domaindata(modeldata, var.domain)
    return get_field(var, domaindata)
end

function get_field(var::VariableDomain, domaindata::AbstractDomainData)    
    return domaindata.variable_data[var.ID]
end

"""
    get_data(var::VariableDomain, modeldata::AbstractModelData) -> field.values
    get_data(var::VariableDomain, modeldata::AbstractDomainData) -> field.values

Get Variable `var` data array from [`Field`](@ref).values
"""
get_data(var::VariableDomain, domainmodeldata) = get_field(var, domainmodeldata).values

"""
    get_data(var::VariableDomain, modeldata::AbstractModelData) -> get_values_output(field.values)
    get_data(var::VariableDomain, modeldata::AbstractDomainData) -> get_values_output(field.values)

Get a sanitized version of Variable `var` data array for storing as output
from [`get_values_output`]@ref)`(`[`Field`](@ref).values`)`
"""
get_data_output(var::VariableDomain, domainmodeldata) = get_values_output(get_field(var, domainmodeldata))

####################################################################
# Create and add to Domain
##################################################################

"create and add to Domain"
function create_VariableDomPropDep(domain, name, master)
    newvar = VariableDomPropDep(domain=domain, name=name, master=master, ID=_next_variable_ID(domain))
    if haskey(domain.variables, name)
        error("attempt to add duplicate VariableDomPropDep $(domain.name).$(name) to Domain")
    end
    domain.variables[name] = newvar
    _reset_master!(newvar, master)
    return newvar
end

 "create and add to Domain"
 function create_VariableDomContribTarget(domain, name, master)
    newvar = VariableDomContribTarget(domain=domain, name=name, master=master, ID=_next_variable_ID(domain))
    if haskey(domain.variables, name)
        error("attempt to add duplicate VariableDomContribTarget $(domain.name).$(name) to Domain")
    end
    domain.variables[name] = newvar
    _reset_master!(newvar, master)
    return newvar
end

####################################################################
# Manage linked VariableReactions
##################################################################

function add_dependency(vardom::VariableDomPropDep, varreact::VariableReaction{VT_ReactDependency})
    push!(vardom.var_dependencies, varreact)
    if get_attribute(varreact, :vfunction) in (VF_StateExplicit, VF_State) 
        @debug "    Resetting master variable"
        _reset_master!(vardom, varreact)
    end   
    return nothing
end

function add_dependency(vardom::VariableDomContribTarget, varreact::VariableReaction{VT_ReactDependency})
    push!(vardom.var_dependencies, varreact)
    vf = get_attribute(varreact, :vfunction)
    if vf in (VF_StateExplicit, VF_State) 
        error("attempt to link Dependency with :vfunction==$vf to a VariableDomContribTarget "*
            "$(fullname(varreact)) --> $(fullname(vardom))")
    end   
    return nothing
end

function add_contributor(vardom::VariableDomPropDep, varreact::VariableReaction{VT_ReactContributor})
    error("configuration error: attempting to link a VariableReactContrib $(fullname(varreact)) "*
        "(usually a flux) which must link to a Target, to a Property VariableDomPropDep $(fullname(vardom))")
    return nothing
end

function add_contributor(vardom::VariableDomContribTarget, varreact::VariableReaction{VT_ReactContributor})
    push!(vardom.var_contributors, varreact)
    if get_attribute(varreact, :vfunction) in (VF_Deriv, VF_Total, VF_Constraint) 
        @debug "    Resetting master variable"
        _reset_master!(vardom, varreact)
    end
    return nothing
end


# reset master variable
function _reset_master!(var::VariableDomain, master::VariableReaction)
    var.master = master
end

function get_all_links(var::VariableDomPropDep)
    all_links = Vector{VariableReaction}()
    append!(all_links, var.var_property)
    append!(all_links,  var.var_dependencies)
    return all_links
end

function get_all_links(var::VariableDomContribTarget)
    all_links = Vector{VariableReaction}()
    if !isnothing(var.var_target)
        push!(all_links, var.var_target)
    end
    append!(all_links,  var.var_contributors)
    append!(all_links,  var.var_dependencies)
    return all_links
end

function get_modifying_methods(var::VariableDomPropDep)
    methods = AbstractReactionMethod[rv.method for rv in var.var_property]
    return methods
end

function get_modifying_methods(var::VariableDomContribTarget)
    methods = AbstractReactionMethod[rv.method for rv in var.var_contributors]
    return methods
end

function get_dependent_methods(var::VariableDomPropDep)
    methods = AbstractReactionMethod[rv.method for rv in var.var_dependencies]
    return methods
end

function get_dependent_methods(var::VariableDomContribTarget)
    methods = AbstractReactionMethod[rv.method for rv in var.var_dependencies]
    if !isnothing(var.var_target)
        push!(methods, var.var_target.method)
    end    
    return methods
end



##############################
# Pretty printing
############################


"compact display form"
function Base.show(io::IO, var::VariableDomain)
    print(io, 
        typeof(var), 
        "(name='", var.name,"'",
        ", hostdep=", host_dependent(var), 
        ")"
    )
end
"multiline display form"
function Base.show(io::IO, ::MIME"text/plain", var::VariableDomain)
    println(io, typeof(var))
    println(io, "  name='", var.name, "'")
    println(io, "  hostdep=", host_dependent(var))
    println(io, "  attributes=", var.attributes)
end

"""
    show_links(vardom::VariableDomain) 

Display all [`VariableReaction`](@ref)s linked to this [`VariableDomain`](@ref)
"""
show_links(vardom::VariableDomain) = show_links(stdout, vardom)

function show_links(io::IO, vardom::VariableDomPropDep)
    println(io, "\t$(typeof(vardom)) $(fullname(vardom)) links:")
    println(io, "\t\tproperty:\t",  fullname.(vardom.var_property))
    for var in vardom.var_dependencies
        println(io, "\t\tdependency:\t", fullname(var))
    end
    return nothing
end

function show_links(io::IO, vardom::VariableDomContribTarget)
    println(io, "\t$(typeof(vardom)) $(fullname(vardom)) links:")
    println(io, "\t\ttarget:\t", if isnothing(vardom.var_target) "nothing" else fullname(vardom.var_target) end)
    for var in vardom.var_contributors
        println(io, "\t\tcontributor:\t", fullname(var))
    end
    for var in vardom.var_dependencies
        println(io, "\t\tdependency:\t", fullname(var))
    end
    return nothing
end
