# Get scalar value from variable x (discarding any AD derivatives)
value_ad(x) = x
# Model code should implement this for any AD types used, eg
# value_ad(x::SparsityTracing.ADval) = SparsityTracing.value(x)
# value_ad(x::ForwardDiff.Dual) = ForwardDiff.value(x)

# get scalar or ad from variable x, as specified by first argument
value_ad(::Type{T}, x::T) where {T} = x        # pass through AD
value_ad(::Type{T}, x::Float64) where {T} = x  # pass through Float64
value_ad(::Type{Float64}, x)  = value_ad(x)    # strip AD
value_ad(::Type{Float64}, x::Float64)  = x     # avoid method ambiguity

# TODO there doesn't seem to be an easy way of parameterising ModelData by an Array type ?
const PaleoArrayType = Array

"""
    AbstractParameter

Base Type for Parameters
    
See also: [`Parameter`](@ref), [`VecParameter`](@ref)
"""
abstract type AbstractParameter
end

#################################
# Reactions
###################################


"""
    AbstractReactionMethod

Base Type for Reaction Methods.

"""
abstract type AbstractReactionMethod
end


"""
    AbstractReaction

Base Type for Reactions.

# Implementation
Derived types should include a field `base::`[`ReactionBase`](@ref)
"""
abstract type AbstractReaction
end


##############################################
# Variables
###############################################

"""
    VariableBase
 
A `Model` biogeochemical `Variable`. `Reaction`s access `Variable`s using derived Types [`VariableReaction`](@ref)
which are links to [`VariableDomain`](@ref)s.
"""
abstract type VariableBase
end

"""
    show_variables(obj, ...) -> Table

Show all Variables attached to PALEO object `obj`
""" 
function show_variables end

"""
    get_variable(obj, varname, ...) -> variable

Get `variable <: VariableBase` by name from PALEO object
""" 
function get_variable end

"""
    has_variable(obj, varname, ...) -> Bool

True if PALEO object containts `varname`
""" 
function has_variable end


"""
    @enum VariableType

Enumeration of `VariableBase` subtypes. Allowed values:
- `VariableReaction`: `VT_ReactProperty`, `VT_ReactDependency`, `VT_ReactContributor`, `VT_ReactTarget`
- `VariableDomain` : `VT_DomPropDep`, `VT_DomContribTarget`
"""
@enum VariableType::Cint begin
    VT_Undefined = 0
    VT_ReactProperty
    VT_ReactDependency
    VT_ReactContributor
    VT_ReactTarget
    VT_DomPropDep
    VT_DomContribTarget
end

"""
    VariableDomain

Host ([`Domain`](@ref)) model variable.

See also: [`VariableDomPropDep`](@ref), [`VariableDomContribTarget`](@ref)
"""
abstract type VariableDomain <: VariableBase
end

"""
    AbstractDomain

A model region containing Fields and Reactions that act on them. 
"""
abstract type AbstractDomain
end

abstract type AbstractSubdomain
end

abstract type AbstractMesh
end

"""
    function get_mesh(obj, ...)

Return an [`AbstractMesh`](@ref) for PALEO object `obj`
"""
function get_mesh end

"""
    AbstractCellRange

Defines a range of cells within a [`Domain`](@ref)
"""
abstract type AbstractCellRange
end

abstract type AbstractSpace
end

abstract type AbstractData
end

"parse eg \"ScalarData\" as ScalarData"
function Base.parse(::Type{AbstractData}, str::AbstractString)   
    dtype = getproperty(@__MODULE__, Symbol(str))

    dtype <: AbstractData || 
        throw(ArgumentError("$str is not a subtype of AbstractData"))
    return dtype
end


"""
    AbstractField

Defines a Field in a Domain
"""
abstract type AbstractField
end

"`model` data arrays etc"
abstract type AbstractModelData
end

"struct to hold Domain Field data"
abstract type AbstractDomainData
end

"""
    function get_table(obj, ...)

Return a Tables.jl data table view for PALEO object `obj`
"""
function get_table end

#############################################
# ReactionMethodDispatchList
#############################################

"""
    ReactionMethodDispatchList

Defines a list of [`ReactionMethod`](@ref) with corresponding [`CellRange`](@ref)
and views on Variable data (sub)arrays.
"""
struct ReactionMethodDispatchList{MF <:Tuple, M <:Tuple, V <:Tuple, C <: Tuple}
    methodfns::MF
    methods::M
    vardatas::V
    cellranges::C
end

# See https://discourse.julialang.org/t/pretty-print-of-type/19555
# Customize typeof function, as full type name is too verbose (as Tuples are of length ~ number of ReactionMethods to call)
Base.show(io::IO, ::Type{PALEOboxes.ReactionMethodDispatchList{MF, M, V, C}}) where {MF, M, V, C} = 
    print(io, "PALEOboxes.ReactionMethodDispatchList{MF::Tuple, M::Tuple, V::Tuple, C::Tuple each of length=$(fieldcount(MF))}")

"compact form"
function Base.show(io::IO, dispatchlist::ReactionMethodDispatchList)
    print(io, typeof(dispatchlist))
end

"multiline form"
function Base.show(io::IO, ::MIME"text/plain", dl::ReactionMethodDispatchList)
    println(io, typeof(dl))
    for i in eachindex(dl.methodfns)
        println(io, "  ", dl.methods[i], ", ", dl.cellranges[i])
    end 
end
