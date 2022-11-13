import Infiltrator
# import StaticArrays
"""
    VariableReaction{T}

Reaction view on a model variable.

Reactions define [`AbstractVarList`](@ref)s of `VariableReaction`s when creating a [`ReactionMethod`](@ref). These
are used to create to [`VariableDomain`](@ref) variables when the model is initialised, 
linking together `VariableReaction`s with the same name. 
`views` on Domain data Arrays are then passed to the [`ReactionMethod`](@ref) function at each timestep.

The Type parameter `T` is one of `VT_ReactProperty`, `VT_ReactDependency`, `VT_ReactContributor`, `VT_ReactTarget`.

`localname` identifies the `VariableReaction` within the `Reaction`, and can be used to set `variable_attributes:`
and `variable_links:` in the .yaml configuration file.

`linkreq_domain.linkreq_subdomain.linkreq_name` defines the Domain, Subdomain and name for run-time linking
to [`VariableDomain`](@ref) variables. 

Not created directly - see [`CreateVariableReaction`](@ref)
"""
Base.@kwdef mutable struct VariableReaction{T} <: VariableBase
    method::Union{Nothing, AbstractReactionMethod} = nothing
    localname::String

    attributes::Dict{Symbol, Any}       = default_variable_attributes()

    "Requested (ie set by config file) named domain-hosted variable to link to"
    linkreq_domain::String              = ""
    linkreq_subdomain::String           = ""
    linkreq_name::String                = ""
    link_optional::Bool                 = false  # VT_ReactDependency, VT_ReactContributor only

    "Link variable"
    linkvar::Union{Nothing, VariableDomain} = nothing
end


get_var_type(var::VariableReaction{T}) where T = T

const VarPropT        = VariableReaction{VT_ReactProperty}
const VarDepT         = VariableReaction{VT_ReactDependency}
const VarTargetT      = VariableReaction{VT_ReactTarget}
const VarContribT     = VariableReaction{VT_ReactContributor}

function Base.copy(v::VariableReaction{T}) where T
    vcopy = VariableReaction{T}(
        method = v.method,
        localname = v.localname,
        attributes = copy(v.attributes), # NB: no deepcopy
        # attributes = Dict{Symbol, Any}(k=>copy(v) for (k, v) in v.attributes), # 
        linkreq_domain = v.linkreq_domain,
        linkreq_subdomain = v.linkreq_subdomain,
        linkreq_name = v.linkreq_name,
        link_optional = v.link_optional,
        linkvar = v.linkvar,
    )
    return vcopy
end

"""
    get_domvar_attribute(var::VariableReaction, name::Symbol, missing_value=missing) -> value

Get the 'master' Variable attribute from the VariableDomain a VariableReaction is linked to.

TODO: this is almost always what is wanted, as only the 'master' VariableReaction will be updated
by the configuration file.
"""
function get_domvar_attribute(var::VariableReaction, name::Symbol, missing_value=missing)
    domvar = var.linkvar
    !isnothing(domvar) || error("get_domvar_attribute: VariableReaction $(fullname(var)) not linked")
        
    return get_attribute(domvar, name, missing_value)
end


###################################################
# Constructors and creation methods
##################################################

"""
    CreateVariableReaction(T, localname => link_namestr, units, description; attributes=Tuple()) -> VariableReaction{T}
    CreateVariableReaction(T, linklocal_namestr, units, description; attributes=Tuple()) -> VariableReaction{T}    
    [deprecated] CreateVariableReaction(T, localname, units, description; link_namestr, attributes=Tuple()) -> VariableReaction{T}

Create a [`VariableReaction`](@ref).

Not called directly: use short names `VarProp`, `VarDep`, `VarContrib`, `VarTarget`, 
`VarPropScalar`, `VarDepScalar`, `VarContribScalar`, `VarTargetScalar`.

# Arguments
- `var_type::VariableType`:  one of `VT_ReactProperty`, `VT_ReactDependency`, `VT_ReactContributor`, `VT_ReactTarget`
- `localname::AbstractString`: Reaction-local Variable name
- `link_namestr::AbstractString`: `<linkreq_domain>.[linkreq_subdomain.]linkreq_name`. Parsed by [`parse_variablereaction_namestr`](@ref)
  to define the requested linking to `Domain` Variable.
- `linklocal_namestr::AbstractString`:  `<linkreq_domain>.[linkreq_subdomain.]localname`. Convenience form to define both `localname`
  and requested linking to Domain Variable, for the common case where `linkreq_name == localname`.
- `units::AbstractString`: units ("" if not applicable)
- `description::AbstractString`: text describing the variable
# Keywords
- `attributes::Tuple(:attrb1name=>attrb1value, :attrb2name=>attrb2value, ...)`: 
  variable attributes, see [`StandardAttributes`](@ref), [`set_attribute!`](@ref), [`get_attribute`](@ref)
"""
function CreateVariableReaction end

CreateVariableReaction(
    var_type::VariableType,
    localname_linkname::Pair{<:AbstractString, <:AbstractString},
    units::AbstractString,
    description::AbstractString;                                      
    attributes::Tuple{Vararg{Pair}}=(),  # (:atrb1=>value1, :atrb2=>value2)
) = CreateVariableReaction(
    var_type, first(localname_linkname), units, description;
    link_namestr=last(localname_linkname), attributes=attributes
)

function CreateVariableReaction(
    var_type::VariableType,
    linklocal_namestr::AbstractString,
    units::AbstractString,
    description::AbstractString;                                      
    attributes::Tuple{Vararg{Pair}}=(),  # (:atrb1=>value1, :atrb2=>value2)
    link_namestr = linklocal_namestr
)

    (_, _, localname, _) = parse_variablereaction_namestr(linklocal_namestr)
    localname = sub_variablereaction_linkreq_name(localname, "")  # strip %reaction%
    if linklocal_namestr != link_namestr
        linklocal_namestr == localname ||
            error("CreateVariableReaction: invalid combination of explicit link_namestr=$link_namestr and localname=$linklocal_namestr")
    end

    (linkreq_domain, linkreq_subdomain, linkreq_name, link_optional) = parse_variablereaction_namestr(link_namestr)

    var_type in (VT_ReactProperty, VT_ReactDependency, VT_ReactContributor, VT_ReactTarget) || 
        error("VariableReaction $name invalid var_type=$var_type")

    newvar = VariableReaction{var_type}(
        localname=localname,                                 
        linkreq_domain=linkreq_domain,
        linkreq_subdomain=linkreq_subdomain,
        linkreq_name=linkreq_name,
        link_optional=link_optional,
    )

    # default initialize_to_zero
    if (get_var_type(newvar) == VT_ReactTarget) 
        set_attribute!(newvar, :initialize_to_zero, true)
    end
    if (get_var_type(newvar) == VT_ReactContributor)
        # do this first before applying attributes to newvar with set_attribute!,
        # so that :initialize_to_zero can be overridden by an entry in attributes
        for (namesymbol, value) in  attributes
            if namesymbol == :vfunction && value in (VF_Deriv, VF_Constraint, VF_Total)                
                # :vfunction host-dependent, so there will be no VT_ReactTarget within the model to initialize_to_zero
                # so set :initialize_to_zero on the VT_ReactContributor
                set_attribute!(newvar, :initialize_to_zero, true)
            end
        end
    end
    
    if ((get_var_type(newvar) == VT_ReactProperty) || (get_var_type(newvar) == VT_ReactTarget)) && link_optional
        error("VariableReaction $name invalid link_optional=$link_optional")
    end

    set_attribute!(newvar, :units, units)
    set_attribute!(newvar, :description, description)

    for (namesymbol, value) in  attributes
        set_attribute!(newvar, namesymbol, value, allow_create=true)
    end
    return newvar
end


"""
    reset_link_namestr!(var, link_namestr)

Reset `link_namestr` after Variable creation, but before Model link_variables
"""
function reset_link_namestr!(var, link_namestr)
    (linkreq_domain, linkreq_subdomain, linkreq_name, link_optional) = 
        parse_variablereaction_namestr(link_namestr)

    var.linkreq_domain = linkreq_domain
    var.linkreq_subdomain = linkreq_subdomain
    var.linkreq_name = linkreq_name
    var.link_optional = link_optional
    return nothing
end


VarPropScalar(localname, units, description; attributes::Tuple=(), kwargs...) =
    VarProp(localname, units, description; attributes=(:space=>ScalarSpace, attributes...), kwargs...)
VarProp(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    CreateVariableReaction(VT_ReactProperty, localname, units, description; attributes=(:field_data=>ScalarData, attributes...), kwargs...)
            
VarPropScalarStateIndep(localname, units, description; attributes::Tuple=(), kwargs... ) =
    VarPropScalar(localname, units, description; attributes=(attributes..., :datatype=>Float64), kwargs...)
VarPropStateIndep(localname, units, description; attributes::Tuple=(), kwargs...) =
    VarProp(localname, units, description; attributes=(attributes..., :datatype=>Float64),  kwargs...)

VarDepScalar(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarDep(localname, units, description; attributes=(:space=>ScalarSpace, attributes...), kwargs...)
VarDepColumn(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarDep(localname, units, description; attributes=(:space=>ColumnSpace, attributes...), kwargs...)
VarDep(localname, units, description; kwargs... ) = 
    CreateVariableReaction(VT_ReactDependency, localname, units, description; kwargs...)
    
VarDepScalarStateIndep(localname, units, description; attributes::Tuple=(), kwargs... ) =
    VarDepScalar(localname, units, description; attributes=(attributes..., :datatype=>Float64), kwargs...)
VarDepColumnStateIndep(localname, units, description; attributes::Tuple=(), kwargs... ) =
    VarDepColumn(localname, units, description; attributes=(attributes..., :datatype=>Float64), kwargs...)
VarDepStateIndep(localname, units, description; attributes::Tuple=(), kwargs...) =
    VarDep(localname, units, description; attributes=(attributes..., :datatype=>Float64),  kwargs...)


# create a VarDep suitable for linking to a VarProp or VarTarget
VarDep(v::VarDepT) = v
function VarDep(v::Union{VarPropT, VarTargetT})
    vdep = VarDepT(        
        localname = v.localname,  
        attributes = copy(v.attributes), # NB: no deepcopy
        linkreq_domain = v.linkreq_domain,
        linkreq_subdomain = v.linkreq_subdomain,
        linkreq_name = v.linkreq_name
    )

    # remove as v will handle this
    has_attribute(vdep, :initialize_to_zero) && set_attribute!(vdep, :initialize_to_zero, false) 
    has_attribute(vdep, :calc_total) && set_attribute!(vdep, :calc_total, false)
    return vdep
end

VarTargetScalar(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarTarget(localname, units, description; attributes=(:space=>ScalarSpace, attributes...), kwargs...)
VarTarget(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    CreateVariableReaction(VT_ReactTarget, localname, units, description; attributes=(:field_data=>ScalarData, attributes...), kwargs...)
        
VarContribScalar(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarContrib(localname, units, description; attributes=(:space=>ScalarSpace, attributes...), kwargs...)
VarContribColumn(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarContrib(localname, units, description; attributes=(:space=>ColumnSpace, attributes...), kwargs...)
VarContrib(localname, units, description; kwargs... ) = 
    CreateVariableReaction(VT_ReactContributor, localname, units, description; kwargs...)

VarContrib(v::VarContribT) = v
VarContrib(v::VarTargetT) = VarContribT(        
    localname = v.localname,  
    attributes = copy(v.attributes), # NB: no deepcopy
    linkreq_domain = v.linkreq_domain,
    linkreq_subdomain = v.linkreq_subdomain,
    linkreq_name = v.linkreq_name
)

VarStateExplicitScalar(localname, units, description; attributes::Tuple=(), kwargs...) =
    VarDepScalar(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_StateExplicit), kwargs...)
VarStateExplicit(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarDep(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_StateExplicit), kwargs...)
    
VarTotalScalar(localname, units, description; attributes::Tuple=(), kwargs...) =
    VarContribScalar(localname, units, description;
        components, attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_Total), kwargs...)
VarTotal(localname, units, description; attributes::Tuple=(), kwargs... ) = 
        VarContrib(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_Total), kwargs...)
         
VarDerivScalar(localname, units, description; attributes::Tuple=(), kwargs... ) =
    VarContribScalar(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_Deriv), kwargs...)
VarDeriv(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarContrib(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_Deriv),  kwargs...)
         
VarConstraintScalar(localname, units, description; attributes::Tuple=(), kwargs... ) =
    VarContribScalar(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_Constraint), kwargs...)
VarConstraint(localname, units, description; attributes::Tuple=(), kwargs... ) = 
        VarContrib(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_Constraint), kwargs...)
         
VarStateScalar(localname, units, description; attributes::Tuple=(), kwargs...) =
    VarDepScalar(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_State), kwargs...)
VarState(localname, units, description; attributes::Tuple=(), kwargs... ) = 
    VarDep(localname, units, description; attributes=(:field_data=>ScalarData, attributes..., :vfunction=>VF_State), kwargs...)
   
# TODO: define a VarInit Type. Currently (ab)using VarDep   
VarInit(v::Union{VarPropT, VarTargetT, VarContribT}) = VarDepT(        
    localname = v.localname,
    attributes = copy(v.attributes), # NB: no deepcopy
    linkreq_domain = v.linkreq_domain,
    linkreq_subdomain = v.linkreq_subdomain,
    linkreq_name = v.linkreq_name
)

function VarInit(v::VarDepT)
    error("VarInit for VarDepT $(fullname(v))")
end

"""
    VarVector(variable_ctorfn, variables_list) -> Vector{VariableReaction{T}}

Create and return a `Vector{VariableReaction{T}}` defined by specified `variable_ctorfn`
from `variables_list = [(linklocal_namestr, units, description), ...]``
"""
VarVector(variable_ctorfn, variables_list) = 
    [
        variable_ctorfn(linklocal_namestr, units, description) 
        for (linklocal_namestr, units, description) in variables_list
    ]

#########################################################################
# VarLists
#########################################################################

"""
    AbstractVarList

Variables required by a [`ReactionMethod`](@ref) `methodfn` are specified by 
a Tuple of VarList_xxx <: AbstractVarList, each containing a collection of [`VariableReaction`](@ref).

These are then converted (by the [`create_accessors`](@ref) method) to a corresponding
Tuple of collections of views on Domain data arrays , which are then be passed to the [`ReactionMethod`](@ref) `methodfn`.

# Implementation
Subtypes of `AbstractVarList` should implement:
- a constructor that takes a collection of VariableReactions
- [`create_accessors`](@ref), returning the views on Domain data arrays in a subtype-specific collection.
- [`get_variables`](@ref), returning the collection of VariableReactions (as a flat list).
"""
abstract type AbstractVarList
end

"""
    create_accessors(varlist::AbstractVarList, modeldata::AbstractModelData, arrays_idx::Int) -> vardata

Return a collection `vardata` of views on Domain data arrays for VariableReactions in `varlist`.
Collection and view are determined by `varlist` Type.
"""
function create_accessors(va::AbstractVarList, modeldata::AbstractModelData, arrays_idx::Int) end

"""
    get_variables(varlist::AbstractVarList) -> Vector{VariableReaction}

Return VariableReaction in `varlist` as a flat Vector.
"""
function get_variables(va::AbstractVarList) end


"""
    VarList_single(var; components=false) -> VarList_single

Create a `VarList_single` describing a single `VariableReaction`,
`create_accessors` will then return a single accessor.
"""
struct VarList_single <: AbstractVarList
    var::VariableReaction
    components::Bool
    VarList_single(var; components=false) = new(copy(var), components)
end

get_variables(vl::VarList_single) = [vl.var]

create_accessors(vl::VarList_single, modeldata::AbstractModelData, arrays_idx::Int) = 
    create_accessor(vl.var, modeldata, arrays_idx, vl.components)


"""
    VarList_components(varcollection) -> VarList_components

Create a `VarList_components` describing a collection of `VariableReaction`s,
`create_accessors` will then return a Vector with data array components concatenated (flattened).
"""
struct VarList_components <:  AbstractVarList
    vars::Vector{VariableReaction}
    allow_unlinked::Bool
    VarList_components(varcollection; allow_unlinked=false) = 
        new([copy(v) for v in varcollection], allow_unlinked)
end

get_variables(vl::VarList_components) = vl.vars

function create_accessors(vl::VarList_components, modeldata::AbstractModelData, arrays_idx::Int)
    accessors_generic = []
    for v in vl.vars
        a = create_accessor(v, modeldata, arrays_idx, true)
        if isnothing(a) || isempty(a)
            if vl.allow_unlinked
                append!(accessors_generic, fill(nothing, num_components(v)))
            else
                error("create_accessors(::VarList_components, ) - unlinked variable $(fullname(v)) $v")
            end
        else
            append!(accessors_generic, a)
        end
    end
    
    accessors_typed = [a for a in accessors_generic]  # narrow Type
    return accessors_typed
end

# _dynamic_svector(av) = _dynamic_svector(Val(length(av)), av)
# _dynamic_svector(::Val{N}, av) where {N} = StaticArrays.SVector{N}(av)

struct VarList_namedtuple <: AbstractVarList
    vars::Vector{VariableReaction}
    keys::Vector{Symbol}
    components::Bool
end

"""
    VarList_namedtuple(varcollection; components=false) -> VarList_namedtuple

Create a `VarList_namedtuple` describing a collection of `VariableReaction`s,
`create_accessors` will then return a NamedTuple with field names = `VariableReaction.localname`
and field values = corresponding data arrays.

If `components = true`, each NamedTuple field will be a Vector of data array components.
"""
function VarList_namedtuple(varcollection; components=false)
    keys = Symbol.([v.localname for v in varcollection])
    vars = [copy(v) for v in varcollection]
 
    return VarList_namedtuple(vars, keys, components)
end

"""
    VarList_namedtuple_fields(objectwithvars; components=false) -> VarList_namedtuple

Create a `VarList_namedtuple` describing `VariableReaction` fields in `objectwithvars`,
`create_accessors` will then return a NamedTuple with field names = field names in `objectwithvars` and
field values = corresponding data arrays.

If `components = true`, each NamedTuple field will be a Vector of data array components.
"""
function VarList_namedtuple_fields(objectwithvars; components=false)
    # find field names and VariableReactions
    fieldnamesvars = [(f, getproperty(objectwithvars, f)) for f in propertynames(objectwithvars) 
                        if getproperty(objectwithvars, f) isa VariableReaction]

    keys      = [f for (f,v) in fieldnamesvars]
    vars      = [copy(v) for (f,v) in fieldnamesvars]
    return VarList_namedtuple(vars, keys, components)
end


get_variables(vl::VarList_namedtuple) = vl.vars

create_accessors(vl::VarList_namedtuple, modeldata::AbstractModelData, arrays_idx::Int) =
    NamedTuple{Tuple(vl.keys)}(create_accessor(v, modeldata, arrays_idx, vl.components) for v in vl.vars)

"""
    VarList_tuple(varcollection; components=false) -> VarList_tuple

Create a `VarList_tuple` describing a collection of `VariableReaction`s,
`create_accessors` will then return a Tuple of data arrays.

If `components = true`, each Tuple field will be a Vector of data array components.
"""
struct VarList_tuple <: AbstractVarList
    vars::Vector{VariableReaction}
    components::Bool
    VarList_tuple(varcollection; components=false) = new([copy(v) for v in varcollection], components)
end

get_variables(vl::VarList_tuple) = vl.vars

create_accessors(vl::VarList_tuple, modeldata::AbstractModelData, arrays_idx::Int) =
    Tuple(create_accessor(v, modeldata, arrays_idx, vl.components) for v in vl.vars)

"""
    VarList_vector(varcollection; components=false, forceview=false) -> VarList_vector

Create a `VarList_vector` describing a collection of `VariableReaction`s,
`create_accessors` will then return a Vector of data arrays.

If `components = true`, each Vector element will be a Vector of data array components.

If `forceview = true`, each accessor will be a 1-D `view` to help type stability,
even if this is redundant (ie no `view` required, v::Vector -> view(v, 1:length(v)))
"""
struct VarList_vector <: AbstractVarList
    vars::Vector{VariableReaction}
    components::Bool
    forceview::Bool
    VarList_vector(varcollection; components=false, forceview=false) =
        new([copy(v) for v in varcollection], components, forceview)
end

get_variables(vl::VarList_vector) = vl.vars

create_accessors(vl::VarList_vector, modeldata::AbstractModelData, arrays_idx::Int) = 
    [create_accessor(v, modeldata, vl.components, arrays_idx, forceview=vl.forceview) for v in vl.vars]


"""
    VarList_vvector(Vector{Vector{VariableReaction}}::vars; components=false) -> VarList_vvector

Create a `VarList_vvector` describing a Vector of Vectors of `VariableReaction`s,
`create_accessors` will then return a Vector of Vectors of data arrays.

If `components = true`, each Vector of Vectors element will be a Vector of data array components.
"""
struct VarList_vvector <: AbstractVarList
    vars::Vector{Vector{VariableReaction}}
    components::Bool
    VarList_vvector(vars; components=false) = new([[copy(v) for v in vv] for vv in vars], components)
end

get_variables(vl::VarList_vvector) = vcat(vl.vars...)

create_accessors(vl::VarList_vvector, modeldata::AbstractModelData, arrays_idx::Int) = 
    [[create_accessor(v, modeldata, arrays_idx, vl.components) for v in vv] for vv in vl.vars]

"""
    VarList_nothing() -> VarList_nothing

Create a placeholder for a missing/unavailable VariableReaction.
`create_accessors` will then return `nothing`.
"""   
struct VarList_nothing <: AbstractVarList
end

get_variables(vl::VarList_nothing) = VariableReaction[] 
create_accessors(vl::VarList_nothing, modeldata::AbstractModelData, arrays_idx::Int) = nothing

"""
    VarList_tuple_nothing(nvar) -> VarList_tuple_nothing

Create a placeholder for `nvar` missing/unavailable VariableReactions.
`create_accessors` will then return an `NTuple{nvar, Nothing}`.
"""   
struct VarList_tuple_nothing <: AbstractVarList
    nvar::Int
end

get_variables(vl::VarList_tuple_nothing) = VariableReaction[] 
create_accessors(vl::VarList_tuple_nothing, modeldata::AbstractModelData, arrays_idx::Int) = ntuple(x->nothing, vl.nvar)

"""
    VarList_fields(varcollection) -> VarList_fields

Create a `VarList_fields` describing a collection of `VariableReaction`s,
`create_accessors` will then return a Tuple with unprocessed Fields.
"""
struct VarList_fields <: AbstractVarList
    vars::Vector{VariableReaction}
    VarList_fields(varcollection) = new([copy(v) for v in varcollection])
end

get_variables(vl::VarList_fields) = vl.vars

function create_accessors(vl::VarList_fields, modeldata::AbstractModelData, arrays_idx::Int)
    accessors = []
    for v in vl.vars
        isempty(v.linkreq_subdomain) || error("variable $v subdomains not supported")
        if isnothing(v.linkvar)
            v.link_optional || error("unlinked variable $v")
            push!(accessors, nothing)
        else
            push!(accessors, get_field(v.linkvar, modeldata, arrays_idx))    
        end
    end
    return Tuple(accessors)
end

# convert Tuple of VarLists to Tuple of data array views
create_accessors(method::AbstractReactionMethod, modeldata::AbstractModelData, arrays_idx::Int) =
    Tuple(create_accessors(vl, modeldata, arrays_idx) for vl in method.varlists)

 """
     create_accessor(var::VariableReaction, modeldata, arrays_idx, components, [,forceview=false])
        -> accessor or (accessor, subdomain_indices)

Creates a `view` on a (single) [`VariableDomain`](@ref) data array linked by `var::`[`VariableReaction`](@ref).
Called by an [`AbstractVarList`](@ref) [`create_accessors`](@ref)
implementation to generate a collection of `views` for multiple [`VariableReaction`](@ref)s.

Returns:
- if `var` is linked, an `accessor` or Tuple `(accessor, subdomain_indices)` that provides a `view` on variable data.
- if `var` is not linked, `nothing` if `var` is optional, or errors and doesn't return if `var` is non-optional.

Mapping of indices for Subdomains <--> Domains:
- if no Subdomain, returns unmodified indices (if `forceview=false`), 
    or an equivalent view (if `forceview=true`, this is to help type stability) 
- if Variable is a Domain interior and Subdomain is a Domain boundary,
    `accessor` is a view with a subset of Domain indices.
- if Variable is a Domain boundary and Subdomain is the Domain interior,
    returns a Tuple with `subdomain_indices`,  length(Subdomain size), with `missing` for interior points.

Mapping of multi-component (Isotope) Variables:
- If `components=false`:
    -  map multi-component Variable to `accessor::IsotopeArray`
    - return a single-component Variable as a `accessor::AbstractArray`. 
- If `components=true`:
    - variable data as a `accessor::Vector{Array}`, length=number of components
"""
function create_accessor(
    var::VariableReaction, modeldata::AbstractModelData, arrays_idx::Int, components; 
    forceview=false
)
    errstring = "create_accessor: VariableReaction $(fullname(var))  ->  $(isnothing(var.linkvar) ? nothing : fullname(var.linkvar))"
    @debug errstring
    if isnothing(var.linkvar)
        var.link_optional || error("$errstring:  unlinked variable")
        return nothing
    else
        linkvar_field = get_field(var.linkvar, modeldata, arrays_idx)

        # find subdomain (if requested)
        if !isempty(var.linkreq_subdomain)
            !isnothing(var.linkvar.domain.grid) || error("$errstring:   linkreq_subdomain='$(var.linkreq_subdomain)' but Domain has no grid defined") 
            linksubdomain = Grids.get_subdomain(var.linkvar.domain.grid, var.linkreq_subdomain)
            !isnothing(linksubdomain) || error("$errstring:  linkreq_subdomain='$(var.linkreq_subdomain)' not found") 
            subdomain_indices = Grids.subdomain_indices(linksubdomain) # nothing, unless accessing boundary from interior
        else
            linksubdomain = nothing
            subdomain_indices = nothing
        end

        try
            var_accessor = create_accessor(
                get_attribute(var, :field_data),
                # TODO - check space, data_dims,
                linkvar_field, linksubdomain,
                forceview=forceview, components=components,
            )

            !isnothing(var_accessor) || error("$errstring:   create_accessor returned $nothing")
    
            if isnothing(subdomain_indices) 
                return var_accessor  # usual case: no subdomain, or accessing interior from boundary as a view
            else
                return (var_accessor, subdomain_indices) # accessing boundary from interior
            end
        catch ex
            if isa(ex, MethodError)
                @warn "$errstring:  invalid :field_data combination - attempting to link a Variable "*
                    "with :field_data=$(get_attribute(var, :field_data)) to a Variable with :field_data=$(field_data(linkvar_field))"
            end
            rethrow(ex)
        end        
    end
end



###################################################################
# Get / set / modify variable data arrays
# (optional, provide convenient ways of handling unlinked variables)
####################################################################

@Base.propagate_inbounds function add_if_available(vardata, idx, flux)
    vardata[idx] += flux
    return nothing
end
add_if_available(vardata::Nothing, idx, flux) = nothing

@Base.propagate_inbounds function add_if_available(vardata, flux)
    vardata[] += flux
    return nothing
end
add_if_available(vardata::Nothing, flux) = nothing

@Base.propagate_inbounds function set_if_available(vardata, idx, val)
    vardata[idx] = val
    return nothing
end
set_if_available(vardata::Nothing, idx, val) = nothing

@Base.propagate_inbounds function set_if_available(vardata, val)
    vardata[] = val
    return nothing
end
set_if_available(vardata::Nothing, val) = nothing

@Base.propagate_inbounds get_if_available(vardata, defaultval) = vardata[]
get_if_available(vardata::Nothing, defaultval) = defaultval
@Base.propagate_inbounds get_if_available(vardata, idx, defaultval) = vardata[idx]
get_if_available(vardata::Nothing, idx, defaultval) = defaultval


###########################################################
# Parsing of variable names into domain, subdomain etc 
########################################################
"Split a name into component parts.
    name = domain.subdomain.name
"
function split_link_name(fullname::AbstractString)
    splitname = split(fullname, ".")
    if length(splitname) == 1
        (domain, subdomain, name) = ("",            "",         splitname[1])        
    elseif length(splitname) == 2
        (domain, subdomain, name) = (splitname[1],  "",         splitname[2]) 
    elseif length(splitname) == 3
        (domain, subdomain, name) = (splitname[1], splitname[2], splitname[3])          
    else
        error("invalid format $fullname")
    end

    return (domain, subdomain, name)
end

"Combine component parts to form name"
function combine_link_name(domain::AbstractString, subdomain::AbstractString, name::AbstractString; sep=".")
    fullname = name
    if !isempty(subdomain)
        fullname = subdomain*sep*fullname
    end
    if !isempty(domain)
        fullname = domain*sep*fullname
    end
    return fullname
end

function combine_link_name(var::VariableReaction)
    return combine_link_name(var.linkreq_domain, var.linkreq_subdomain, var.linkreq_name)
end

function sub_variablereaction_linkreq_name(linkreq_name::AbstractString, reactionname::AbstractString)
    return replace(linkreq_name, "%reaction%" => reactionname)
end


function strip_brackets(linkstr::AbstractString)
    if linkstr[1]=='(' && linkstr[end] == ')'
        hasbrackets = true
        linkstr_core = linkstr[2:end-1]
    else
        hasbrackets = false
        linkstr_core = linkstr
    end

    return (hasbrackets, linkstr_core)
end

"""
    parse_variablereaction_namestr(linkstr) 
        -> (linkreq_domain, linkreq_subdomain, linkreq_name, link_optional)

Parse a linkstr into component parts.
 
`linkstr` is of format: `[(][<linkreq_domain>.][<linkreq_subdomain>.]<linkreq_name>[)]`
  * Optional brackets `( ... )` set `link_optional=true`
  * `linkreq_name` may contain `%reaction%` which will later be substituted with `<Reaction name>/`
 
# Examples:
```jldoctest; setup = :(import PALEOboxes)
julia> PALEOboxes.parse_variablereaction_namestr("foo")  # Common case eg for a property that should be public
("", "", "foo", false)

julia> PALEOboxes.parse_variablereaction_namestr("%reaction%foo")  # Reaction-private by default
("", "", "%reaction%foo", false)

julia> PALEOboxes.parse_variablereaction_namestr("ocean.foo")  # Request link to variable of same name in ocean Domain
("ocean", "", "foo", false)

julia> PALEOboxes.parse_variablereaction_namestr("(ocean.oceansurface.goo)") # Full syntax
("ocean", "oceansurface", "goo", true)
```
"""
function parse_variablereaction_namestr(linkstr::AbstractString)

    any(occursin.([' ', '>', '<'], linkstr)) && 
        error("invalid character in linkstr='$linkstr'")

    (link_optional, slink) = strip_brackets(linkstr)
    
    (linkreq_domain, linkreq_subdomain, linkreq_name) = 
        split_link_name(slink)
    
    return (linkreq_domain, linkreq_subdomain, linkreq_name, link_optional)
end


#####################################
# Pretty printing
########################################

"fully qualified name (domainname.reactionname.methodname.varname)"
function fullname(var::VariableReaction)
    if isnothing(var.method)
        return "<no method>."*var.localname
    else
        return fullname(var.method)*"."*var.localname
    end
end

"link name + link status"
function show_link_name(var::VariableReaction)
    link_name = combine_link_name(var)
    if isnothing(var.linkvar)
        link_name = link_name*"[unlinked]"
    end
    return link_name
end

show_links(var::VariableReaction) = show_links(stdout, var)

function show_links(io::IO, var::VariableReaction)
    if isnothing(var.linkvar)
        println(io, "\t$(typeof(var)) $(fullname(var)) not linked")
    else
        println(io, "\t$(typeof(var)) $(fullname(var)) --> $(fullname(var.linkvar))")
        show_links(var.linkvar)
    end
    return nothing
end

show_links(vars::AbstractVector) = foreach(show_links, vars)

"compact form"
function Base.show(io::IO, var::VariableReaction)
    print(io, typeof(var), "(localname='$(var.localname)', link_name='$(show_link_name(var))')")
end
"multiline form"
function Base.show(io::IO, ::MIME"text/plain", var::VariableReaction)
    println(io, typeof(var))
    println(io, "  localname='$(var.localname)'")
    println(io, "  link_name='$(show_link_name(var))'")
    println(io, "  attributes=", var.attributes)
end





