"""
    VariableDomain <: VariableBase

Abstract base Type for a ([`Domain`](@ref)) model variable.
    
Defines a named Variable and corresponding data [`Field`](@ref)s that are linked to by [`VariableReaction`](@ref)s.

See [`VariableDomPropDep`](@ref), [`VariableDomContribTarget`](@ref) 
"""
VariableDomain

"""
    VariableDomPropDep <: VariableDomain

[`Model`](@ref) ([`Domain`](@ref)) `VariableDomain` linking a single `VariableReaction{VT_ReactProperty}`
to multiple `VariableReaction{VT_ReactDependency}`.
"""
Base.@kwdef mutable struct VariableDomPropDep <: VariableDomain
    ID::Int64
    domain::AbstractDomain
    name::String
    master::VariableReaction

    var_property::Union{Nothing, VariableReaction{VT_ReactProperty}} =
        nothing
    var_property_setup::Union{Nothing, VariableReaction{VT_ReactProperty}} =
        nothing
    var_dependencies::Vector{VariableReaction{VT_ReactDependency}} =
        Vector{VariableReaction{VT_ReactDependency}}()

end

function get_properties(var::VariableDomPropDep)
    pvars = VariableReaction[]
    isnothing(var.var_property) || push!(pvars, var.var_property)
    isnothing(var.var_property_setup) || push!(pvars, var.var_property_setup)
    return pvars
end

get_var_type(var::VariableDomPropDep) = VT_DomPropDep

"""
    VariableDomContribTarget <: VariableDomain

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
    return isnothing(var.var_property) && isnothing(var.var_property_setup)
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

#############################################################
# Allocate data arrays
#############################################################

"""
    allocate_variables!(
        vars, modeldata; 
        [, eltypemap::Dict{String, DataType}],
        [, default_host_dependent_field_data=nothing])

Allocate memory for [`VariableDomain`](@ref)s `vars`.

Element type of allocated Arrays is determined by `eltype(modeldata)` (the usual case, allowing use of AD types), 
which can be overridden by Variable `:datatype` attribute if present (allowing Variables to ignore AD types).
`:datatype` may be either a Julia `DataType` (eg Float64), or a string to be looked up in `eltypemap`.

Field data type (<:[`AbstractFieldData`](@ref)) is determined by Variable `:field_data` attribute, optionally this can take a
`default_host_dependent_field_data` default for Variables with `host_dependent(v)==true` (these are Variables with no Target or no Property linked,
intended to be external dependencies supplied by the solver).
"""
function allocate_variables!(
    vars, modeldata::AbstractModelData; 
    eltypemap=Dict{String, DataType}(),
    default_host_dependent_field_data=nothing,
)
    
    for v in vars
        check_lengths(v)

        data_dims = Tuple(
            get_data_dimension(v.domain, dimname) 
            for dimname in get_attribute(v, :data_dims)
        )
    
        # eltype usually is eltype(modeldata), but can be overridden by :datatype attribute 
        # (can be used by a Reaction to define a fixed type eg Float64 for a constant Property)
        mdeltype = get_attribute(v, :datatype, eltype(modeldata))
        if mdeltype isa AbstractString
            mdeltype = get(eltypemap, mdeltype, Float64)
        end
        
        thread_safe = false
        if get_attribute(v, :atomic, false) && modeldata.threadsafe
            @info "  $(fullname(v)) allocating Atomic data"
            thread_safe = true
        end

        field_data = get_attribute(v, :field_data)
        space = get_attribute(v, :space)

        if field_data == UndefinedData
            if host_dependent(v) && (get_attribute(v, :vfunction, VF_Undefined) == VF_Undefined) && !isnothing(default_host_dependent_field_data)               
                set_attribute!(v, :field_data, default_host_dependent_field_data)
                field_data = get_attribute(v, :field_data)
                @info "    set :field_data=$field_data for host-dependent Variable $(fullname(v))"
            else
                error("allocate_variables! :field_data=UndefinedData for Variable $(fullname(v)) $v")
            end
        end
        v_field = allocate_field(
            field_data, data_dims, mdeltype, space, v.domain.grid,
            thread_safe=thread_safe, allocatenans=modeldata.allocatenans
        )        
        
        set_field!(v, modeldata, v_field)
      
    end

    return nothing
end

"""
    reallocate_variables!(vars, modeldata, new_eltype) -> [(v, old_eltype), ...]

Reallocate memory for [`VariableDomain`](@ref)s `vars` to `new_eltype`. Returns Vector of
`(reallocated_variable, old_eltype)`.
"""
function reallocate_variables!(vars, modeldata, new_eltype)
    
    reallocated_variables = []
    for v in vars
        v_data = get_data(v, modeldata)        
        if v_data isa AbstractArray
            old_eltype = eltype(v_data)
            if old_eltype != new_eltype            
                set_data!(v, modeldata, similar(v_data, new_eltype))
                push!(reallocated_variables, (v, old_eltype))
            end
        end
    end

    return reallocated_variables
end

"""
    check_lengths(var::VariableDomain)

Check that sizes of all linked Variables match
"""
function check_lengths(var::VariableDomain)

    var_space = get_attribute(var, :space)
    var_scalar = (var_space == ScalarSpace)
    
    for lv in get_all_links(var)
        if get_attribute(lv, :check_length, true)
            var_size = var_scalar ? (1, ) : internal_size(var.domain.grid, lv.linkreq_subdomain)

            link_space = get_attribute(lv, :space)
            link_scalar = (link_space == ScalarSpace)
            link_size = link_scalar ? (1, ) : internal_size(lv.method.domain.grid)

            var_size == link_size || 
                error("check_lengths: VariableDomain $(fullname(var)), :space=$var_space size=$var_size "*
                    "!= $(fullname(lv)), :space=$link_space size=$link_size created by $(typename(lv.method.reaction)) (check size of Domains $(var.domain.name), $(lv.method.domain.name), "*
                    "$(isempty(lv.linkreq_subdomain) ? "" : "subdomain "*lv.linkreq_subdomain*",") and Variables :space)")
        end
    end

    return nothing
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
    append!(all_links,  get_properties(var))
    append!(all_links,  var.var_dependencies)
    return all_links
end

function get_all_links(var::VariableDomContribTarget)
    all_links = Vector{VariableReaction}()
    isnothing(var.var_target) || push!(all_links, var.var_target)
    append!(all_links,  var.var_contributors)
    append!(all_links,  var.var_dependencies)
    return all_links
end

function get_modifying_methods(var::VariableDomPropDep)
    methods = AbstractReactionMethod[rv.method for rv in get_properties(var)]
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
    show_links(io::IO, vardom::VariableDomPropDep)

Display all [`VariableReaction`](@ref)s linked to this [`VariableDomain`](@ref)
"""
show_links(vardom::VariableDomain) = show_links(stdout, vardom)

function show_links(io::IO, vardom::VariableDomPropDep)
    println(io, "\t$(typeof(vardom)) $(fullname(vardom)) links:")
    println(io, "\t\tproperty:\t", if isnothing(vardom.var_property) "nothing" else fullname(vardom.var_property) end)
    println(io, "\t\tproperty_setup:\t", if isnothing(vardom.var_property_setup) "nothing" else fullname(vardom.var_property_setup) end)
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
