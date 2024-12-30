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
function is_allocated(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int)
    return !isnothing(get_field(var, modeldata, arrays_idx))    
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

function get_dimensions(v::VariableDomain; expand_cartesian=false)
    dims = NamedDimension[]
    append!(dims, PB.get_dimensions(v.domain.grid, v.attributes[:space]))
    for ddn in v.attributes[:data_dims]
        append!(get_data_dimension(v.domain, ddn))
    end
    return dims
end

################################################
# Field and data access
##################################################

"""
    set_field!(var::VariableDomain, modeldata, arrays_idx::Int, field::Field)
 
Set `VariableDomain` data to `field`
"""
function set_field!(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int, field::Field)

    variable_data = get_domaindata(modeldata, var.domain, arrays_idx).variable_data

    variable_data[var.ID] = field

    return nothing
end

"""
    set_data!(var::VariableDomain, modeldata, arrays_idx::Int, data)
 
Set [`VariableDomain`](@ref) to a newly-created [`Field`](@ref) containing `data`.
"""
function set_data!(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int, data)

    variable_data = get_domaindata(modeldata, var.domain, arrays_idx).variable_data

    variable_data[var.ID] =  Field(
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
    get_field(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int) -> field::Field
    get_field(var::VariableDomain, domaindata::AbstractDomainData) -> field::Field

Get Variable `var` `field`
"""
function get_field(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int=1)
    domaindata = get_domaindata(modeldata, var.domain, arrays_idx)
    return get_field(var, domaindata)
end

function get_field(var::VariableDomain, domaindata::AbstractDomainData)
    return domaindata.variable_data[var.ID]
end

"""
    get_data(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int=1) -> field.values
    get_data(var::VariableDomain, domaindata::AbstractDomainData) -> field.values

Get Variable `var` data array from [`Field`](@ref).values
"""
get_data(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int=1) = get_field(var, modeldata, arrays_idx).values
get_data(var::VariableDomain, domaindata::AbstractDomainData) = get_field(var, domaindata).values

"""
    get_data_output(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int) -> get_values_output(field.values)
    get_data_output(var::VariableDomain, domaindata::AbstractDomainData) -> get_values_output(field.values)

Get a sanitized version of Variable `var` data array for storing as output
from [`get_values_output`]@ref)`(`[`Field`](@ref).values`)`
"""
get_data_output(var::VariableDomain, modeldata::AbstractModelData, arrays_idx::Int=1) = get_values_output(get_field(var, modeldata, arrays_idx))
get_data_output(var::VariableDomain, domaindata::AbstractDomainData) = get_values_output(get_field(var, domaindata))

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
        vars, modeldata, arrays_idx; 
        [, eltypemap::Dict{String, DataType}],
        [, default_host_dependent_field_data=nothing],
        [, allow_base_link=true],
        [. use_base_transfer_jacobian=true],
        [, use_base_vars=String[]],
        [, check_units_opt=:no])

Allocate or link memory for [`VariableDomain`](@ref)s `vars` in `modeldata` array set `arrays_idx`

Element type of allocated Arrays is determined by `eltype(modeldata, arrays_idx)` (the usual case, allowing use of AD types), 
which can be overridden by Variable `:datatype` attribute if present (allowing Variables to ignore AD types).
`:datatype` may be either a Julia `DataType` (eg Float64), or a string to be looked up in `eltypemap`.

If `allow_base_link==true`, and any of the following are true a link is made to the base array (`arrays_idx=1`), 
instead of allocating a new array in array set `arrays_idx`:
    - Variable element type matches `modeldata` base eltype (arrays_idx=1)
    - `use_base_transfer_jacobian=true` and Variable `:transfer_jacobian` attribute is set
    - Variable full name is in `use_base_vars`

Field data type is determined by Variable `:field_data` attribute, optionally this can take a
`default_host_dependent_field_data` default for Variables with `host_dependent(v)==true` (these are Variables with no Target or no Property linked,
intended to be external dependencies supplied by the solver).

If `check_units_opt != :no` then the `:units` field of linked variable is checked, resulting in either a warning (if `check_units_opt=:warn`)
or error (if `check_units_opt=:error`).
"""
function allocate_variables!(
    vars, modeldata::AbstractModelData, arrays_idx::Int; 
    eltypemap=Dict{String, DataType}(),
    default_host_dependent_field_data=ScalarData,
    allow_base_link=true,
    use_base_transfer_jacobian=true,
    use_base_vars=String[],
    check_units_opt=:no,
)
    
    for v in vars
        check_lengths(v)

        check_units(v; check_units_opt)

        data_dims = Tuple(
            get_dimension(v.domain, dimname) 
            for dimname in get_attribute(v, :data_dims)
        )
    
        # eltype usually is eltype(modeldata, arrays_idx), but can be overridden by :datatype attribute 
        # (can be used by a Reaction to define a fixed type eg Float64 for a constant Property)
        mdeltype = get_attribute(v, :datatype, eltype(modeldata, arrays_idx))
        if mdeltype isa AbstractString
            mdeltype_str = mdeltype
            mdeltype = get(eltypemap, mdeltype_str, Float64)
            @debug "Variable $(fullname(v)) mdeltype $mdeltype_str -> $mdeltype"
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

        # allocate or link
        if (arrays_idx != 1 && allow_base_link) && 
                (
                    mdeltype == eltype(modeldata, 1) ||
                    (use_base_transfer_jacobian && get_attribute(v, :transfer_jacobian, false)) ||
                    fullname(v) in use_base_vars
                )
            # link to existing array
            v_field = get_field(v, modeldata, 1)
        else
            # allocate new field array
            v_field = Field(
                field_data, data_dims, mdeltype, space, v.domain.grid;
                allocatenans=modeldata.allocatenans
            )
        end
        
        set_field!(v, modeldata, arrays_idx, v_field)
      
    end

    return nothing
end

"""
    reallocate_variables!(vars, modeldata, arrays_idx, new_eltype) -> [(v, old_eltype), ...]

Reallocate memory for [`VariableDomain`](@ref)s `vars` to `new_eltype`. Returns Vector of
`(reallocated_variable, old_eltype)`.
"""
function reallocate_variables!(vars, modeldata::AbstractModelData, arrays_idx::Int, new_eltype)
    
    reallocated_variables = []
    for v in vars
        v_data = get_data(v, modeldata)        
        if v_data isa AbstractArray
            old_eltype = eltype(v_data)
            if old_eltype != new_eltype            
                set_data!(v, modeldata, arrays_idx, similar(v_data, new_eltype))
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
    
    for lv in get_all_links(var)
        if get_attribute(lv, :check_length, true)
            try
                var_size = internal_size(var_space, var.domain.grid, lv.linkreq_subdomain)

                link_space = get_attribute(lv, :space)
                link_size = internal_size(link_space, lv.method.domain.grid)

                var_size == link_size || 
                    error("check_lengths: VariableDomain $(fullname(var)), :space=$var_space size=$var_size "*
                        "!= $(fullname(lv)), :space=$link_space size=$link_size created by $(typename(lv.method.reaction)) (check size of Domains $(var.domain.name), $(lv.method.domain.name), "*
                        "$(isempty(lv.linkreq_subdomain) ? "" : "subdomain "*lv.linkreq_subdomain*",") and Variables :space)")
            catch
                error("check_lengths: exception VariableDomain $(fullname(var)), :space=$var_space "*
                        "linked by $(fullname(lv)), created by $(typename(lv.method.reaction)) (check size of Domains $(var.domain.name), $(lv.method.domain.name), "*
                        "$(isempty(lv.linkreq_subdomain) ? "" : "subdomain "*lv.linkreq_subdomain*",") and Variables :space)")
            end
        end
    end

    return nothing
end

"""
    check_units(var::VariableDomain; check_units_opt=:warn)

Check that units of all linked Variables match
"""
function check_units(var::VariableDomain; check_units_opt=:warn)

    check_units_opt in (:no, :warn, :error) ||
        error("check_units(): unsupported option check_units_opt=$check_units_opt (allowed values are :no, :warn, :error)")

    check_units_opt == :no && return

    var_units = get_attribute(var, :units)

    num_errors = 0
    
    for lv in get_all_links(var)
        lv_units = get_attribute(lv, :units)
        
        if !_compare_units(var_units, lv_units)
            num_errors += 1 
            @warn "check_units: VariableDomain $(fullname(var)), :units=\"$var_units\" (from master $(typename(var.master.method.reaction)) $(fullname(var.master)))"*
                  " != $(fullname(lv)), :units=\"$lv_units\" (created by $(typename(lv.method.reaction)) $(fullname(lv.method.reaction)))"
        end
    end

    if check_units_opt == :error && !iszero(num_errors)
        error("check_units: VariableDomain $(fullname(var)), :units=$var_units units of linked variables do not match")
    end

    return num_errors
end

function _compare_units(units1, units2)
    # very crude regularization of unit strings:  m^3 and m3 are both accepted
    units1 = replace(units1, "^"=>"")
    units2 = replace(units2, "^"=>"") 
    return (units1 == units2) || (units1 == "unknown") || (units2 == "unknown")
end

####################################################################
# Manage linked VariableReactions
##################################################################

function add_dependency(vardom::VariableDomPropDep, varreact::VariableReaction{VT_ReactDependency})
    push!(vardom.var_dependencies, varreact)
    if get_attribute(varreact, :vfunction) in (VF_StateExplicit, VF_StateTotal, VF_State) 
        @debug "    Resetting master variable"
        _reset_master!(vardom, varreact)
    end   
    return nothing
end

function add_dependency(vardom::VariableDomContribTarget, varreact::VariableReaction{VT_ReactDependency})
    push!(vardom.var_dependencies, varreact)
    vf = get_attribute(varreact, :vfunction)
    if vf in (VF_StateExplicit, VF_StateTotal, VF_State) 
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
    show_links(model::Model, varnamefull::AbstractString)
    show_links(io::IO, vardom::VariableDomain)
    show_links(io::IO, model::Model, varnamefull::AbstractString)

Display all [`VariableReaction`](@ref)s linked to this [`VariableDomain`](@ref)

`varnamefull` should be of form "<domain name>.<variable name>"

Linked variables are shown as "<domain name>.<reaction name>.<method name>.<local name>"
"""
function show_links end

show_links(vardom::VariableDomain) = show_links(stdout, vardom)

# implementation of show_links(model::Model, varnamefull::AbstractString) is in Model.jl 

function show_links(io::IO, vardom::VariableDomPropDep)
    println(io, "\t$(typeof(vardom)) \"$(fullname(vardom))\" links:")
    println(io, "\t\tproperty:\t", if isnothing(vardom.var_property) "nothing" else "\""*fullname(vardom.var_property)*"\"" end)
    if !isnothing(vardom.var_property_setup)
        println(io, "\t\tproperty_setup:\t", "\""*fullname(vardom.var_property_setup)*"\"")
    end
    println(io, "\t\tdependencies:\t", String[fullname(var) for var in vardom.var_dependencies])
    return nothing
end

function show_links(io::IO, vardom::VariableDomContribTarget)
    println(io, "\t$(typeof(vardom)) \"$(fullname(vardom))\" links:")
    println(io, "\t\ttarget:\t\t", if isnothing(vardom.var_target) "nothing" else "\""*fullname(vardom.var_target)*"\"" end)
    println(io, "\t\tcontributors:\t", String[fullname(var) for var in vardom.var_contributors])
    println(io, "\t\tdependencies:\t", String[fullname(var) for var in vardom.var_dependencies])
    return nothing
end
