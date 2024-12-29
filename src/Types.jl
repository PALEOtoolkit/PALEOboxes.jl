# TODO there doesn't seem to be an easy way of parameterising ModelData by an Array type ?
const PaleoArrayType = Array

################################
# CommonDataModel adaptor
################################

"""
    CDModel(x)

Create a CommonDataModel adaptor for PALEO object x
"""
function CDModel end

################################
# Parameters
###############################

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


abstract type VariableDomain <: VariableBase
end

abstract type AbstractModel
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

const AbstractMeshOrNothing = Union{AbstractMesh, Nothing}

"""
    function get_mesh(obj, ...)

Return an [`AbstractMesh`](@ref) for PALEO object `obj`
"""
function get_mesh end

abstract type AbstractCellRange
end

abstract type AbstractSpace
end

abstract type AbstractData
end

"parse eg \"ScalarData\" or \"PALEOboxes.ScalarData\" as ScalarData"
function Base.parse(::Type{AbstractData}, str::AbstractString)   
    dtype = tryparse(AbstractData, str)

    !isnothing(dtype) || 
        throw(ArgumentError("$str is not a subtype of AbstractData"))
    return dtype
end

function Base.tryparse(::Type{AbstractData}, str::AbstractString)   
    dtype = getproperty(@__MODULE__, Symbol(replace(str, "PALEOboxes."=>"")))

    return (dtype <: AbstractData) ? dtype : nothing
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
struct ReactionMethodDispatchList{M <:Tuple, V <:Tuple, C <: Tuple}
    methods::M
    vardatas::V
    cellranges::C
end

ReactionMethodDispatchList(methods::Vector, vardatas::Vector, cellranges::Vector) = 
    ReactionMethodDispatchList(Tuple(methods), Tuple(vardatas), Tuple(cellranges))

struct ReactionMethodDispatchListNoGen
    methods::Vector
    vardatas::Vector
    cellranges::Vector
end

# See https://discourse.julialang.org/t/pretty-print-of-type/19555
# Customize typeof function, as full type name is too verbose (as Tuples are of length ~ number of ReactionMethods to call)
Base.show(io::IO, ::Type{PALEOboxes.ReactionMethodDispatchList{M, V, C}}) where {M, V, C} = 
    print(io, "PALEOboxes.ReactionMethodDispatchList{M::Tuple, V::Tuple, C::Tuple each of length=$(fieldcount(M))}")

function Base.show(io::IO, dispatchlist::ReactionMethodDispatchList)
    print(io, typeof(dispatchlist))
end

function Base.show(io::IO, ::MIME"text/plain", dl::ReactionMethodDispatchList)
    println(io, typeof(dl))
    for i in eachindex(dl.methodfns)
        println(io, "  ", dl.methods[i], ", ", dl.cellranges[i])
    end 
end

"""
    infoerror(io::IOBuffer, message::AbstractString)

Output accumulated log messages in io, then raise `ErrorException` with message
"""
function infoerror(io::IOBuffer, message::AbstractString)
    s = String(take!(io))
    isempty(s) || @info s
    error(message)
end