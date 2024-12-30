

############################################################
# (Function) Spaces
###########################################################
"""
    AbstractSpace

Defines a function Space within a Domain, on a mesh defined by a Grid 
"""
AbstractSpace

"""
    ScalarSpace <: AbstractSpace

A Domain position-independent quantity
"""
struct ScalarSpace <: AbstractSpace
end

"""
    CellSpace <: AbstractSpace

A per-cell quantity. Use as Variable attribute :space to create a Variable with data array dimensions from Grid cells.
"""
struct CellSpace <: AbstractSpace
end

"""
    ColumnSpace <: AbstractSpace

A per-column quantity. Use as Variable attribute :space to create a Variable with data array dimensions from Grid columns.
"""
struct ColumnSpace <: AbstractSpace
end


"""
    Face1DColumnSpace <: AbstractSpace

A quantity defined on upper and lower faces of a cell in a 1D column
"""
struct Face1DColumnSpace <: AbstractSpace
end

"parse eg \"CellSpace\" as CellSpace"
function Base.parse(::Type{AbstractSpace}, str::AbstractString)   
    dtype = getproperty(@__MODULE__, Symbol(str))
    dtype <: AbstractSpace || 
        throw(ArgumentError("$str is not a subtype of AbstractSpace"))
    return dtype
end

################################################################
# AbstractMesh Spaces and sizes
#
# Concrete types (UnstructuredVectorGrid, CartesianLinearGrid, ...) 
# should implement internal_size, and optionally cartesian_size
#################################################################

"""
    has_internal_cartesian(mesh::AbstractMesh, Space::Type{<:AbstractSpace})::Bool

`true` if `Mesh` uses different internal and external (cartesian) spatial array representations for `Space`
"""
function has_internal_cartesian end


"""
    internal_size(Space::Type{<:AbstractSpace}, mesh::AbstractMesh; [subdomain=""] [space=:cell]) -> NTuple{ndims, Int}

Array size to use for model Variables.

All `AbstractMesh` concrete subtypes (UnstructuredVectorGrid, CartesianLinearGrid, ...) should implement this method.

# Optional Keyword Arguments
- `subdomain::String=""`: a named subdomain  
"""
function internal_size end

"""
    cartesian_size(Space::Type{<:AbstractSpace}, mesh::AbstractMesh) -> NTuple{ndims, Int}

Optional (only regular Cartesian grids should implement this method): Array size of Cartesian Domain.

NB: this may be different from `internal_size` if the `mesh` implements a mapping 
eg to a Vector for internal model Variables.
"""
function cartesian_size end


"""
    spatial_size(::Type{<:AbstractSpace}, mesh::AbstractMesh) -> NTuple{ndims, Int}

Array size for given Space and mesh.
"""
spatial_size(Space::Type{<:AbstractSpace}, mesh) = internal_size(Space, mesh)


################################################################
# AbstractData interface 
#
# Concrete types (ScalarData, ArrayScalarData, IsotopeData) should implement these methods
#################################################################
"""
    AbstractData

Defines a Data type that can be composed with an [`AbstractSpace`](@ref) to form a Field

Concrete subtypes should implement:

[`allocate_values`](@ref), [`check_values`](@ref), [`zero_values!`](@ref), [`dof_values`](@ref),
[`get_values_output`](@ref)

If the subtype needs to provide values for a numerical solver (eg as a state variable), it also needs to implement:

[`init_values!`](@ref), [`copyfieldto!`](@ref), [`copytofield!`](@ref), [`add_field!`](@ref), [`add_field_vec!`](@ref)

If the subtype has a representation as components, it should implement:
[`num_components`](@ref), [`get_components`](@ref)

If the subtype needs to provide a thread-safe atomic addition operation eg to provide scalar accumulator variables for Domain totals with a tiled model,
it should implement [`atomic_add!`](@ref) for the field values.
"""
AbstractData


"""
    allocate_values(
        FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, data_type, Space::Type{<:AbstractSpace}, spatial_size::Tuple;
        allocatenans::Bool,
    ) -> values::V

allocate `Field.values::V` (eg an Array of elements with Type `data_type`) for `FieldData` with dimensions defined by `spatial_size` and `data_dims`
"""
function allocate_values(
    FieldData::Type{<:AbstractData}, 
    data_dims::Tuple{Vararg{NamedDimension}}, 
    data_type, 
    Space::Type{<:AbstractSpace}, 
    spatial_size::Tuple{Integer, Vararg{Integer}}, # an NTuple with at least one element
    allocatenans::Bool,
) 
end

"""
    check_values(
        existing_values, 
        FieldData::Type{<:AbstractData},
        data_dims::Tuple{Vararg{NamedDimension}}, 
        data_type, 
        Space::Type{<:AbstractSpace}, 
        spatial_size::Tuple{Integer, Vararg{Integer}}
    )

Check `existing_values` is of suitable type, size etc for use as `Field.values`, throw exception if not.
"""
function check_values(
    existing_values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, data_type, Space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}},
)
end

"""
    check_data_type(existing_values, data_type::Union{Missing, Type})

Helper function for [`check_values`](@ref) implementations: check `existing_values` are consistent
with `data_type` (if supplied), throw exception if not.
"""
check_data_type(existing_values, data_type::Missing) = nothing

function check_data_type(existing_values, data_type::Type)
    existing_eltype = eltype(existing_values)
    existing_eltype === data_type ||
        throw(ArgumentError("data_type mismatch: supplied $existing_eltype, require $data_type"))
    return nothing
end

"""
    init_values!(
        values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace},
        init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange, info::NTuple{3, String}
    ) 
    
Initialize `values` at model start to `init_value` over region `cellrange` using information from
Variable `attribv` attributes, scaled by `convertfn` and `convertvalues`. 
 
Optional: only required if this `FieldData` type is used for a model (state) Variable that requires initialisation.

Arguments:
- `values`: data to be zeroed
- `init_value::Symbol`: one of :initial_value, :norm_value, requesting type of initial value required
- `attribv::VariableBase`: Variable with attributes to use for initialisation
- `convertfn::Function`:  apply multiplier `convertfn(convertvalues, i)` to initialisation value for cell i.
  Typically this is used to convert units eg concentration to mol.
- `convertvalues`: parameters (if any) required by `convertfn`, eg a volume measure.
- `cellrange`: range of cells to initialise
- `info::::NTuple{3, String}`: Tuple (varinfo, convertinfo, trsfrinfo) of identifier strings to use for log messages
  
"""
function init_values!(
    values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace},
    init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange, info::NTuple{3, String}
)
end

"""
    zero_values!(values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, cellrange)

Set `values` over spatial region `cellrange` to zero at start of main loop 
"""
function zero_values!(values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, cellrange)
end

"""
    field_single_element(FieldData::Type{<:AbstractData}, N) -> Bool

Return true if `FieldData` with length(data_dims) = N is represented with 
a single value that can be accessed as [] (used to optimize FieldRecord storage).

Default is probably OK, unless `FieldData` uses a Vector of values per element.
"""
function field_single_element(FieldData::Type{<:AbstractData}, N)
    if N == 0
        return true
    else
        return false
    end
end

"""
    dof_values(
        FieldData::Type{<:AbstractData}, 
        data_dims::Tuple{Vararg{NamedDimension}}, 
        Space::Type{<:AbstractSpace}, 
        mesh, 
        cellrange
    ) -> dof::Int

Return degrees-of-freedom for `FieldData` over spatial region `cellrange`.
"""
function dof_values(FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, mesh, cellrange) 
end

"""
    copyfieldto!(
        dest,
        doff, 
        values, 
        FieldData::Type{<:AbstractData}, 
        data_dims::Tuple{Vararg{NamedDimension}}, 
        Space::Type{<:AbstractSpace}, 
        cellrange
    ) -> num_copied::Int

Copy Field.values `values` from spatial region defined by `cellrange`, to `dest` Array starting at index `doff`.

Number of values over whole Domain should equal degrees-of-freedom returned by [`dof_values`](@ref)

Required if this `FieldData` type needs to provide values for a numerical solver.
"""
function copyfieldto!(dest, doff, values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, cellrange)
end

"""
    copytofield!(
        values,
        FieldData::Type{<:AbstractData},
        data_dims::Tuple{Vararg{NamedDimension}},
        Space::Type{<:AbstractSpace},
        cellrange, 
        src, 
        soff
    ) -> num_copied::Int

Copy from `src` Array starting at index `soff` to Field.values `values` for spatial region defined by `cellrange`.

Number of values over whole Domain should equal degrees-of-freedom returned by [`dof_values`](@ref)

Required if this `FieldData` type needs to provide values for a numerical solver.
"""
function copytofield!(values, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, cellrange, src, soff)
end

"""
    add_field!(dest, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, a, cellrange, src)

Implement `dest += a*src` where `dest`, `src` are Field.values, `a` is a number, over region defined by `cellrange`
"""
function add_field!(dest, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, a, cellrange, src)
end

"""
    add_field_vec!(
        dest, 
        FieldData::Type{<:AbstractData}, 
        data_dims::Tuple{Vararg{NamedDimension}}, 
        Space::Type{<:AbstractSpace}, 
        a, 
        cellrange,
        src, 
        soff
    ) -> num_added::Int

Implement `dest += a*src` where `dest` is a Field.values, `src` is an Array, `a` is a number, over region defined by `cellrange`,
starting at index `soff` in `src`.

Returns number of elements of `src` used.

See  [`copytofield!`](@ref), [`copyfieldto!`](@ref) for the relationship between Array `src` and Field values `dest`.
"""
function add_field_vec!(dest, FieldData::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space::Type{<:AbstractSpace}, a, cellrange, src, soff)
end

"""
    num_components(FieldData::Type{<:AbstractData}) -> Int

get number of components (optional - implement if `FieldData` has a representation as components)
"""
function num_components(FieldData::Type{<:AbstractData}) end

"""
    get_components(values, FieldData::Type{<:AbstractData}) -> Vector

Convert Field `values` to a Vector of components
(optional - implement if `FieldData` has a representation as components)
"""
function get_components(values, FieldData::Type{<:AbstractData}) end

"Optional: sanitize `values` for storing as model output.
Default implementation is usually OK - only implement for custom types that should be converted to standard types for storage"
get_values_output(values, data_type::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, Space, mesh) = values



################################################################
# UndefinedData
#################################################################


"""
    UndefinedData <: AbstractData

Undefined data type (no methods implemented). Used to indicate that a Variable can link to any data type.
"""
struct UndefinedData <: AbstractData
end


#############################################################
# Field
#############################################################


"""
    Field{FieldData <: AbstractData, Space <: AbstractSpace, V, N, Mesh <: AbstractMeshOrNothing}
    Field(FieldData, data_dims::NTuple{N, NamedDimensions}, data_type, Space, mesh::AbstractMeshOrNothing; allocatenans)
    Field(existing_values, FieldData, data_dims::NTuple{N, NamedDimension}, data_type::Union{DataType, Missing}, Space, mesh::AbstractMeshOrNothing)

A Field of `values::V` of type `FieldData <: AbstractData` defined on function space `Space` over grid of type `Mesh`
and (optionally) with `N` `data_dims::NTuple{N, NamedDimensions}`.

`values::V` is usually Array or similar of elements with Type `data_type` and dimensions defined by `Space`, `Mesh` and `data_dims`,
and is either created by [`allocate_values`](@ref) or supplied by `existing_values`.
"""
struct Field{FieldData <: AbstractData, Space <: AbstractSpace, V, N, Mesh <: AbstractMeshOrNothing}
    values::V
    data_dims::NTuple{N, NamedDimension}
    mesh::Mesh
end

# create a new Field, allocating `values` data arrays
function Field(
    FieldData::Type, data_dims::NTuple{N, NamedDimension}, data_type::Type, Space::Type{<:AbstractSpace}, mesh::AbstractMeshOrNothing;
    allocatenans
) where {N}
    v = allocate_values(
        FieldData, data_dims, data_type, Space, spatial_size(Space, mesh);
        allocatenans,
    )

    return Field{FieldData, Space, typeof(v), N, typeof(mesh)}(v, data_dims, mesh)
end

# create a new Field, containing supplied `existing_values` data arrays
function Field(
    existing_values, FieldData::Type, data_dims::NTuple{N, NamedDimension}, data_type::Union{DataType, Missing}, Space::Type{<:AbstractSpace}, mesh::AbstractMeshOrNothing
) where {N}
    check_values(
        existing_values, FieldData, data_dims, data_type, Space, spatial_size(Space, mesh), 
    )
    return Field{FieldData, Space, typeof(existing_values), N, typeof(mesh)}(existing_values, data_dims, mesh)
end

field_data(field::Field{FieldData, Space, V, N, Mesh}) where {FieldData, Space, V, N, Mesh} = FieldData
space(field::Field{FieldData, Space, V, N, Mesh}) where {FieldData, Space, V, N, Mesh} = Space

function get_dimensions(f::Field; expand_cartesian=false)
    dims = NamedDimension[nd for nd in get_dimensions(f.mesh, space(f); expand_cartesian)]
    append!(dims, f.data_dims)
    return dims
end

# use default compact form Base.show

# multiline form
function Base.show(io::IO, ::MIME"text/plain", f::Field{FieldData, Space, V, N, Mesh}) where {FieldData, Space, V, N, Mesh}
    println(io, "Field{", FieldData, ", ", Space, ", ", V, ", ", N, ", ", Mesh, "}")
    println(io, "  values: ", f.values)
    println(io, "  data_dims: ", f.data_dims)
    println(io, "  mesh: ", f.mesh)
    println(io, "  dimensions:")
    for nd in get_dimensions(f)
        println(io, "    ", nd)
    end
    return nothing
end

"""
    get_field(obj, ...) -> Field

Get Field from PALEO object `obj`
"""
function get_field end

"""
    add_field!(obj, f::Field ...)

Add Field or Field to PALEO object `obj`
"""
function add_field! end



"zero out `field::Field` over region defined by `cellrange`"
function zero_field!(field::Field{FieldData, Space, V, N, Mesh}, cellrange) where {FieldData, Space, V, N, Mesh}
    zero_values!(field.values, FieldData, field.data_dims, Space, cellrange)
end

"initialize `field::Field` to `init_value` (`:initial_value` or `:norm_value`)
from Variable `attribv` attributes, over region defined by `cellrange`.
Optionally calculate transformed initial values from supplied `convertfn` and  `convertvalues`"
function init_field!(
    field::Field{FieldData, Space, V, N, Mesh}, init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange, info,
) where {FieldData, Space, V, N, Mesh}
    init_values!(field.values, FieldData, field.data_dims, Space, init_value, attribv, convertfn, convertvalues, cellrange, info)
end


"calculate number of degrees-of-freedom for `field::Field` over region defined by `cellrange`"
function dof_field(field::Field{FieldData, Space, V, N, Mesh}, cellrange) where {FieldData, Space, V, N, Mesh}
    return dof_values(FieldData, field.data_dims, Space, field.mesh, cellrange)
end

"copy `src::Field`` to `dest::Vector`, optionally restricting to region defined by `cellrange`"
function Base.copyto!(dest, doff, src::Field{FieldData, Space, V, N, Mesh}, cellrange) where {FieldData, Space, V, N, Mesh}
    return copyfieldto!(dest, doff, src.values, FieldData, src.data_dims, Space, cellrange)
end

"copy `src::Vector`` to `dest::Field`, optionally restricting to region defined by `cellrange`"
function Base.copyto!(dest::Field{FieldData, Space, V, N, Mesh}, cellrange, src, soff) where {FieldData, Space, V, N, Mesh}
    return copytofield!(dest.values, FieldData, dest.data_dims, Space, cellrange, src, soff)
end

"Calculate `dest::Field = dest::Field + a * src::Field`, optionally restricting to region defined by `cellrange`"
function add_field!(dest::Field{FieldData, Space, V, N, Mesh}, a, cellrange, src::Field{FieldData, Space, V, N, Mesh}) where {FieldData, Space, V, N, Mesh}
    return add_field!(dest.values, FieldData, dest.data_dims, Space, a, cellrange, src.values)
end

function add_field_vec!(dest::Field{FieldData, Space, V, N, Mesh}, a, cellrange, srcvalues::AbstractVector, soff) where {FieldData, Space, V, N, Mesh}
    return add_field_vec!(dest.values, FieldData, dest.data_dims, Space, a, cellrange, srcvalues, soff)
end

"sanitized version of `values`, suitable for storing as output"
function get_values_output(field::Field{FieldData, Space, V, N, Mesh}) where {FieldData, Space, V, N, Mesh}
    return get_values_output(field.values, FieldData, field.data_dims, Space, field.mesh)
end


# get values from `linkvar_field`, optionally applying view defined by `linksubdomain`
function create_accessor(
    output_data::Union{Type{FieldData}, Type{UndefinedData}},
    # TODO - check space, data_dims,
    linkvar_field::Field{FieldData, Space, V, N, Mesh}, linksubdomain::Union{Nothing, AbstractSubdomain};
    forceview, components,
) where {FieldData, Space, V, N, Mesh}

    # create accessor
    if isnothing(linksubdomain)
        if forceview
            if components
                var_accessor = [view(vc, 1:length(vc)) for vc in get_components(linkvar_field.values, FieldData)]
            else
                var_accessor = view(linkvar_field.values, 1:length(linkvar_field.values))
            end
        else
            if components
                var_accessor = get_components(linkvar_field.values, FieldData)
            else
                var_accessor = linkvar_field.values
            end
        end
    else              
        if components
            var_accessor = [Grids.subdomain_view(vc, linksubdomain) for vc in get_components(linkvar_field.values, FieldData)]
        else
            var_accessor = Grids.subdomain_view(linkvar_field.values, linksubdomain)
        end
    end

    return var_accessor
end
