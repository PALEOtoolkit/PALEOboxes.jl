
import StructArrays

#############################
# AbstractIsotopeScalar
################################

"""
    AbstractIsotopeScalar <: AbstractData

An IsotopeScalar represents a quantity or flux with isotopic composition.
Can be added, subtracted, multiplied by a scalar, and decomposed into components with the
same (bulk) transport properties.

# Implementation
Each IsotopeScalar should be added to `IsotopeTypes`, and implement:
- get_total(is::IsotopeScalar) -> total
- arithmetic operations +/- (ie IsotopeScalars can be added and subtracted) and *number, /number 
  (IsotopeScalars can be scaled by a real number).
- methods of the [`AbstractData`](@ref) interface

and optional isotope-specific functions, eg for a single isotope:
- isotope_totaldelta(::Type{<: IsotopeScalar}, total, delta) -> IsotopeScalar()
- get_delta(is::IsotopeScalar) -> delta

"""
abstract type AbstractIsotopeScalar <: AbstractData
end

num_components(is::AbstractIsotopeScalar)       = num_components(typeof(is))
num_components(iv::AbstractArray)               = num_components(eltype(iv))

Base.length(is::AbstractIsotopeScalar)          = 1
Base.size(is::AbstractIsotopeScalar)          = ()


"""
    @isotope_totaldelta(IsotopeType::Type, total, delta) -> Type

Convenience macro to create either an `IsotopeXXX` of `IsotopeType::Type`, 
or a Float64 if `IsotopeType==ScalarData` (in which case `delta` argument is ignored).

This generates code equivalent to:
    IsotopeType == ScalarData ? total : isotope_totaldelta(IsotopeType, total, delta)
"""
macro isotope_totaldelta(IsotopeType, total, delta)
    return :($(esc(IsotopeType)) == ScalarData ? $(esc(total)) : isotope_totaldelta($(esc(IsotopeType)), $(esc(total)), $(esc(delta)) ) )
end

"generic get_total for non-isotope variable"
function get_total(scalar)
    return scalar
end

"""
    get_property(data; propertyname::Union{Symbol, Nothing}=nothing, default_total=false) -> data

Generic wrapper for test scripts etc to provide access to field `propertyname` (which can be `nothing`)
from struct-valued `data`, or to always apply [`get_total`](@ref).

This is intended for use in test scripts etc, where it is useful to have a function that can take arguments
to control access either to a property or to the raw `data`.  This is not efficient (it allocates), not
recommended for use in application code.
"""
function get_property(
    data;
    propertyname::Union{Symbol, Nothing}=nothing,
    default_total=false
)

    if !isnothing(propertyname)
        if eltype(data) <: AbstractVector
            data = [getproperty.(d, propertyname) for d in data]
        else
            data = getproperty.(data, propertyname)
        end
    elseif default_total
        if eltype(data) <: AbstractVector
            data = [get_total.(d) for d in data]
        else
            data = get_total.(data)
        end
    end

    return data
end

#####################################################
# IsotopeLinear
#####################################################
"""
    IsotopeLinear <: AbstractIsotopeScalar

Linearized representation of isotopic composition, where 
`moldelta` = `total` * `delta`.
"""
struct IsotopeLinear{T1, T2} <: AbstractIsotopeScalar
    v::T1
    v_moldelta::T2
end

function Base.getproperty(obj::IsotopeLinear, sym::Symbol)
    if sym === :v_delta
        return get_delta(obj)
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

function Base.propertynames(obj::IsotopeLinear, private::Bool=false)
    return (:v, :v_moldelta, :v_delta)
end

"compact form"
function Base.show(io::IO, il::IsotopeLinear)
    print(io, "(v=", il.v,", v_moldelta=", il.v_moldelta,", ‰=", il.v_delta, ")")
end


"""
    isotope_totaldelta(::Type{IsotopeLinear}, total, delta) -> IsotopeLinear

Create an `IsotopeLinear` from `total` and `delta`

# Examples:
```jldoctest; setup = :(import PALEOboxes as PB)
julia> a = PB.isotope_totaldelta(PB.IsotopeLinear, 10.0, -2.0)
(v=10.0, v_moldelta=-20.0, ‰=-2.0)

julia> b = a*2
(v=20.0, v_moldelta=-40.0, ‰=-2.0)

julia> c = a + b
(v=30.0, v_moldelta=-60.0, ‰=-2.0)

julia> PB.get_total(c)
30.0

julia> PB.get_delta(c)
-2.0
```
"""
isotope_totaldelta(::Type{IsotopeLinear}, total, delta) = IsotopeLinear(total, total*delta)
isotope_totaldelta(::Type{IsotopeLinear{T, T}}, total::T, delta::T) where {T} = IsotopeLinear(total, total*delta)

"""
    get_total(is::IsotopeLinear)

Get isotope total
"""
function get_total(is::IsotopeLinear)
    return is.v
end

"""
    get_delta(is::IsotopeLinear)

Get isotope delta (per mil)
"""
function get_delta(is::IsotopeLinear, add_eps_denom=0.0)
    delta =  is.v_moldelta / (is.v + add_eps_denom) # (is.v + 1e-3*sign(is.v)*abs(is.v_moldelta)) 
    return delta # min(max(delta, -100.0), 100.0)
end

function get_delta_limit(is::IsotopeLinear, add_eps_denom, delta_limit)
    if is.v < add_eps_denom
        delta = 0.0
    else
        delta =  is.v_moldelta / is.v # (is.v + 1e-3*sign(is.v)*abs(is.v_moldelta)) 
    end
    return min(max(delta, -delta_limit), delta_limit)
end

Base.zero(x::IsotopeLinear{T1, T2}) where{T1,T2}= zero(typeof(x))
Base.zero(::Type{IsotopeLinear{T1, T2}}) where{T1,T2} = IsotopeLinear(zero(T1), zero(T2))
Base.eltype(::Type{IsotopeLinear{T1, T2}}) where {T1, T2} = IsotopeLinear{T1, T2}
Base.:-(a::IsotopeLinear)                       = IsotopeLinear(-a.v, -a.v_moldelta)
Base.:-(a::IsotopeLinear,  b::IsotopeLinear)    = IsotopeLinear(a.v-b.v, a.v_moldelta-b.v_moldelta)
Base.:+(a::IsotopeLinear,  b::IsotopeLinear)    = IsotopeLinear(a.v+b.v, a.v_moldelta+b.v_moldelta)
Base.:*(a::IsotopeLinear,  b::Real)             = IsotopeLinear(a.v*b, a.v_moldelta*b)
Base.:*(a::Real, b::IsotopeLinear)              = IsotopeLinear(a*b.v, a*b.v_moldelta)
Base.:/(a::IsotopeLinear,  b::Real)             = IsotopeLinear(a.v/b, a.v_moldelta/b)

# convert methods for cases (eg from Float64 to an AD type) where scalar components v, v_moldelta can be converted
IsotopeLinear{T1, T2}(x::IsotopeLinear) where {T1, T2} = IsotopeLinear{T1, T2}(x.v, x.v_moldelta)
Base.convert(::Type{IsotopeLinear{T1, T2}}, x::IsotopeLinear) where {T1, T2} = IsotopeLinear{T1, T2}(x)

"""
    IsotopeTypes::Tuple

Available [`AbstractIsotopeScalar`](@ref) types
"""
const IsotopeTypes = (ScalarData, IsotopeLinear)


"""
    split_nameisotope(
        nameisotope::AbstractString, isotope_data::Dict=Dict(); 
        default=UndefinedData,
    ) -> (name::AbstractString, field_data::AbstractData)

Split `nameisotope` string of form  <name>[::XIsotope],  and look up `field_data` from key `XIsotope` in supplied Dict `isotope_data`.

Returns ScalarData if ::XIsotope not present in `nameisotope`, `field_data = default` if `XIsotope` key is not present in `isotope_data`.
    
# Examples
```jldoctest; setup = :(import PALEOboxes as PB)
julia> PB.split_nameisotope("myflux::CIsotope", Dict("CIsotope"=>PB.IsotopeLinear))
("myflux", PALEOboxes.IsotopeLinear)

julia> PB.split_nameisotope("myflux::PALEOboxes.IsotopeLinear")
("myflux", PALEOboxes.IsotopeLinear)

julia> PB.split_nameisotope("myflux", Dict("CIsotope"=>PB.IsotopeLinear))
("myflux", PALEOboxes.ScalarData)

julia> PB.split_nameisotope("::CIsotope", Dict("CIsotope"=>PB.IsotopeLinear))
("", PALEOboxes.IsotopeLinear)
```
"""
function split_nameisotope(
    nameisotope::AbstractString, isotope_data::Dict=Dict();
    default=UndefinedData
)    
    convertd(v_data::AbstractString)       = parse(AbstractData, v_data)
    convertd(v_data::Type{D}) where {D <: AbstractData}   = D
    convertd(v_data) = error("split_nameisotope: $v_data isotope_data must contain String or <:AbstractData")

    spfn = split(nameisotope, "::")
    name = spfn[1]
    if length(spfn) == 1            
        field_data = ScalarData
    elseif length(spfn) == 2
        if haskey(isotope_data, spfn[2]) 
            field_data = convertd(isotope_data[spfn[2]])
        elseif !isnothing(tryparse(AbstractData, spfn[2]))
            field_data = parse(AbstractData, spfn[2])
        elseif default != UndefinedData
            field_data = default
        else
            error("flux '$nameisotope' isotope not defined, isotope_data=$isotope_data")
        end
    else
        error("invalid flux '$nameisotope'")
    end

    return (name, field_data)
end


#####################################################
# IsotopeLinear AbstractData implementation
#####################################################

function allocate_values(
    data::Type{IsotopeLinear}, data_dims::Tuple{}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}};
    allocatenans,
)
    v = allocate_values(ScalarData, data_dims, data_type, space, spatial_size; allocatenans)
    v_moldelta = allocate_values(ScalarData, data_dims, data_type, space, spatial_size; allocatenans)

    return StructArrays.StructArray{IsotopeLinear{data_type, data_type}}((v, v_moldelta))
end

function check_values(
    existing_values, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}}
)
    eltype(existing_values) <: IsotopeLinear ||
        throw(ArgumentError("type mismatch: supplied eltype $(eltype(existing_values))"))

    check_data_type(existing_values, data_type)

    existing_size = size(existing_values) 
    existing_size == spatial_size ||
        throw(ArgumentError("size mismatch: supplied $existing_size, require spatial_size=$spatial_size"))

    return nothing
end

function init_values!(
    values, data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{S},
    init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange,
    info::NTuple{3, String}
) where {S <: Union{ScalarSpace, CellSpace}}
    varinfo, convertinfo, trsfrinfo = info

    initial_value = get_attribute(attribv, init_value)

    if init_value == :norm_value
        initial_delta = 1.0
        @info "init_values!     :$(rpad(init_value, 20)) $(rpad(varinfo, 30)) "*
            "= $initial_value$convertinfo, delta=$initial_delta (fixed delta value to calculate norm)$trsfrinfo"
    else
        initial_delta = get_attribute(attribv, :initial_delta)
        @info "init_values!     :$(rpad(init_value, 20)) $(rpad(varinfo, 30)) "*
            "= $initial_value$convertinfo, delta=$initial_delta$trsfrinfo"
    end
  
    # . as initial_value may be a Vector
    i_initial_value = isotope_totaldelta.(IsotopeLinear, initial_value, initial_delta)
    
    if S === ScalarSpace
        values[] = i_initial_value*convertfn(convertvalues, nothing)
    else
        if initial_value isa Real # ie not a Vector
            for i in cellrange.indices
                values[i] = i_initial_value*convertfn(convertvalues,  i)
            end
        elseif length(initial_value) == length(values)
            for i in cellrange.indices
                values[i] = i_initial_value[i]*convertfn(convertvalues,  i)
            end
        else
            error("init_values!: configuration error, Vector Variable $varinfo and Vector :initial_value lengths mismatch ($(length(values)) != $(length(initial_value)))")
        end
    end
    
    return nothing
end

function zero_values!(values, data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{ScalarSpace}, cellrange)
    values[] = zero(eltype(values))
    return nothing
end

function zero_values!(values, data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange)
    @inbounds for i in cellrange.indices
        values[i] = zero(eltype(values))
    end
    return nothing
end

dof_values(data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{ScalarSpace}, mesh, cellrange) = 2
dof_values(data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{CellSpace}, mesh, cellrange) = 2*length(cellrange.indices)
dof_values(data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{CellSpace}, mesh, cellrange::Nothing) = 2*mesh.ncells

function copyfieldto!(dest, doff, values, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{<:AbstractSpace}, cellrange)
    n1 = copyfieldto!(dest, doff, values.v, ScalarData, data_dims, space, cellrange)
    n2 = copyfieldto!(dest, doff+n1, values.v_moldelta, ScalarData, data_dims, space, cellrange)
   
    return n1 + n2
end

function copytofield!(values, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{<:AbstractSpace}, cellrange, src, soff)
    n1 = copytofield!(values.v, ScalarData, data_dims, space, cellrange, src, soff)
    n2 = copytofield!(values.v_moldelta, ScalarData, data_dims, space, cellrange, src, soff+n1)
   
    return n1 + n2
end

add_field!(dest, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{<:AbstractSpace}, a, cellrange, src) =
    add_field!(dest, ScalarData, data_dims::Tuple{}, space, a, cellrange, src)

function add_field_vec!(dest, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{ScalarSpace}, a, cellrange, src, soff)
    dest[].v += a*IsotopeLinear(src[soff], src[soff+1])
    return 2
end

function add_field_vec!(dest, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{CellSpace}, a, cellrange, src, soff)
    # assume src layout matches that of copytofield!
    @inbounds for i in cellrange.indices
        dest[i] += a*IsotopeLinear(src[soff], src[soff+length(cellrange.indices)])
        soff += 1
    end
    return 2*length(cellrange.indices)
end

add_field_vec!(dest, field_data::Type{IsotopeLinear}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange::Nothing, src, soff) =
    add_field_vec!(dest, field_data, data_dims, space, a, (indices=1:length(dest),), src, soff)


num_components(field_data::Type{IsotopeLinear})         = 2
num_components(field_data::Type{IsotopeLinear{T1, T2}})  where {T1,T2}      = num_components(IsotopeLinear)

get_components(values, data::Type{IsotopeLinear}) = [values.v, values.v_moldelta]


# provide a create_accessor method for ScalarData --> IsotopeLinear that takes just the first (total) component
function create_accessor(
    output_data::Type{ScalarData},
    # TODO - check space, data_dims,
    linkvar_field::Field{IsotopeLinear, S, V, W, M}, linksubdomain::Union{Nothing, AbstractSubdomain};
    forceview, components,
) where {S, V, W, M}

    # create accessor
    if isnothing(linksubdomain)
        if forceview
            var_accessor = view(linkvar_field.values.v, 1:length(linkvar_field.values.v))
        else
            var_accessor = linkvar_field.values.v
        end
    else              
        var_accessor = Grids.subdomain_view(linkvar_field.values.v, linksubdomain)
    end

    if components
        return [var_accessor]
    else
        return var_accessor
    end
end

#############################################################
# Special cases for atomic addition
##########################################################

atomic_add!(ias::StructArrays.StructVector{<:IsotopeLinear}, iv::IsotopeLinear) = (atomic_add!(ias.v, iv.v); atomic_add!(ias.v_moldelta, iv.v_moldelta);)

