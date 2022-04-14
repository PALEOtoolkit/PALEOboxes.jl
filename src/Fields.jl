
import Infiltrator
############################################################
# (Function) Spaces
###########################################################

"""
    ScalarSpace

A Domain position-independent quantity
"""
struct ScalarSpace <: AbstractSpace
end

"""
    CellSpace

A per-cell quantity. Use as Variable attribute :space to create a Variable with data array dimensions from Grid
"""
struct CellSpace <: AbstractSpace
end


"""
    Face1DColumnSpace

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
    internal_size(mesh::AbstractMesh) -> NTuple{ndims, Int}

Array size to use for per-cell model Variables.

All `AbstractMesh` concrete subtypes (UnstructuredVectorGrid, CartesianLinearGrid, ...) should implement this method.
"""
function internal_size end

"""
    cartesian_size(mesh::AbstractMesh) -> NTuple{ndims, Int}

Optional (only regular Cartesian grids should implement this method): Array size of Cartesian Domain.

NB: this may be different from `internal_size` if the `mesh` implements a mapping 
eg to a Vector for internal model Variables.
"""
function cartesian_size end

"""
    spatial_size(::Type{<:AbstractSpace}, mesh::AbstractMesh) -> NTuple{ndims, Int}

Array size for given Space and mesh.
"""
spatial_size(::Type{ScalarSpace}, mesh) = (1, )
spatial_size(::Type{CellSpace}, mesh) = internal_size(mesh)


################################################################
# AbstractData interface 
#
# Concrete types (ScalarData, ArrayScalarData, IsotopeData) should implement these methods
#################################################################

"""
    allocate_values(
        field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple;
        thread_safe::Bool, allocatenans::Bool,
    ) -> values

allocate `Field.values` (eg an Array) for `field_data` with dimensions defined by `spatial_size` and `data_dims`
"""
function allocate_values(
    field_data::Type{<:AbstractData}, 
    data_dims::Tuple{Vararg{NamedDimension}}, 
    data_type, 
    space::Type{<:AbstractSpace}, 
    spatial_size::Tuple{Integer, Vararg{Integer}}, # an NTuple with at least one element
    thread_safe::Bool, allocatenans::Bool,
) 
end

"""
    check_values(
        existing_values, 
        field_data::Type{<:AbstractData},
        data_dims::Tuple{Vararg{NamedDimension}}, 
        data_type, 
        space::Type{<:AbstractSpace}, 
        spatial_size::Tuple{Integer, Vararg{Integer}}
    )

Check `existing_values` is of suitable type, size etc for use as `Field.values`, throw exception if not.
"""
function check_values(
    existing_values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}},
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
        values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace},
        init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange, info::NTuple{3, String}
    ) 
    
Initialize `values` at model start to `init_value` over region `cellrange` using information from
Variable `attribv` attributes, scaled by `convertfn` and `convertvalues`. 
 
Optional: only required if this `field_data` type is used for a model (state) Variable that requires initialisation.

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
    values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace},
    init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange, info::NTuple{3, String}
)
end

"""
    zero_values!(values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, cellrange)

Set `values` over spatial region `cellrange` to zero at start of main loop 
"""
function zero_values!(values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, cellrange)
end

"""
    field_single_element(field_data::Type{<:AbstractData}, N) -> Bool

Return true if `field_data` with length(data_dims) = N is represented with 
a single value that can be accessed as [] (used to optimize FieldRecord storage).

Default is probably OK, unless `field_data` uses a Vector of values per element.
"""
function field_single_element(field_data::Type{<:AbstractData}, N)
    if N == 0
        return true
    else
        return false
    end
end

"""
    dof_values(
        field_data::Type{<:AbstractData}, 
        data_dims::Tuple{Vararg{NamedDimension}}, 
        space::Type{<:AbstractSpace}, 
        mesh, 
        cellrange
    ) -> dof::Int

Return degrees-of-freedom for `field_data` over spatial region `cellrange`.
"""
function dof_values(field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, mesh, cellrange) 
end

"""
    copyfieldto!(
        dest,
        doff, 
        values, 
        field_data::Type{<:AbstractData}, 
        data_dims::Tuple{Vararg{NamedDimension}}, 
        space::Type{<:AbstractSpace}, 
        cellrange
    ) -> num_copied::Int

Copy Field.values `values` from spatial region defined by `cellrange`, to `dest` Array starting at index `doff`.

Number of values over whole Domain should equal degrees-of-freedom returned by [`dof_values`](@ref)

Required if this `field_data` type needs to provide values for a numerical solver.
"""
function copyfieldto!(dest, doff, values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, cellrange)
end

"""
    copytofield!(
        values,
        field_data::Type{<:AbstractData},
        data_dims::Tuple{Vararg{NamedDimension}},
        space::Type{<:AbstractSpace},
        cellrange, 
        src, 
        soff
    ) -> num_copied::Int

Copy from `src` Array starting at index `soff` to Field.values `values` for spatial region defined by `cellrange`.

Number of values over whole Domain should equal degrees-of-freedom returned by [`dof_values`](@ref)

Required if this `field_data` type needs to provide values for a numerical solver.
"""
function copytofield!(values, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, cellrange, src, soff)
end

"""
    add_field!(dest, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, a, cellrange, src)

Implement `dest += a*src` where `dest`, `src` are Field.values, `a` is a number, over region defined by `cellrange`
"""
function add_field!(dest, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, a, cellrange, src)
end

"""
    add_field_vec!(
        dest, 
        field_data::Type{<:AbstractData}, 
        data_dims::Tuple{Vararg{NamedDimension}}, 
        space::Type{<:AbstractSpace}, 
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
function add_field_vec!(dest, field_data::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space::Type{<:AbstractSpace}, a, cellrange, src, soff)
end

"""
    num_components(field_data::Type{<:AbstractData}) -> Int

get number of components (optional - implement if `field_data` has a representation as components)
"""
function num_components(field_data::Type{<:AbstractData}) end

"""
    get_components(values, field_data::Type{<:AbstractData}) -> Vector

Convert Field `values` to a Vector of components
(optional - implement if `field_data` has a representation as components)
"""
function get_components(values, field_data::Type{<:AbstractData}) end

"Optional: sanitize `values` for storing as model output.
Default implementation is usually OK - only implement eg for Atomic types that should be converted to standard types for storage"
get_values_output(values, data_type::Type{<:AbstractData}, data_dims::Tuple{Vararg{NamedDimension}}, space, mesh) = values



################################################################
# UndefinedData
#################################################################


"""
    UndefinedData

Undefined data type (no methods implemented). Used to indicate that a Variable can link to any data type.
"""
struct UndefinedData <: AbstractData
end


#############################################################
# Field
#############################################################


"""
    Field{D <: AbstractData, S <: AbstractSpace, V, N, M}

A Field of data type `D` defined on function space `S` over `mesh` and (optionally) `data_dims`
"""
struct Field{D <: AbstractData, S <: AbstractSpace, V, N, M}
    values::V
    data_dims::NTuple{N, NamedDimension}
    mesh::M
end

"""
    get_field(obj, ...) -> Field

Get Field or FieldRecord from PALEO object `obj`
"""
function get_field end

"""
    add_field!(obj, f::Field ...)
    add_field!(obj, fr::FieldRecord ...)

Add Field or FieldRecord to PALEO object `obj`
"""
function add_field! end

"create a new Field, allocating `values` data arrays"
function allocate_field(
    field_data::Type, data_dims::NTuple{N, NamedDimension}, data_type::Type, space::Type{<:AbstractSpace}, mesh;
    thread_safe::Bool, allocatenans
) where {N}
    v = allocate_values(
        field_data, data_dims, data_type, space, spatial_size(space, mesh), 
        thread_safe=thread_safe, allocatenans=allocatenans,
    )

    return Field{field_data, space, typeof(v), N, typeof(mesh)}(v, data_dims, mesh)
end

"create a new Field, containing supplied `existing_values` data arrays"
function wrap_field(
    existing_values, field_data::Type, data_dims::NTuple{N, NamedDimension}, data_type::Union{DataType, Missing}, space::Type{<:AbstractSpace}, mesh
) where {N}
    check_values(
        existing_values, field_data, data_dims, data_type, space, spatial_size(space, mesh), 
    )
    return Field{field_data, space, typeof(existing_values), N, typeof(mesh)}(existing_values, data_dims, mesh)
end

"zero out `field::Field` over region defined by `cellrange`"
function zero_field!(field::Field{D, S, V, N, M}, cellrange) where {D, S, V, N, M}
    zero_values!(field.values, D, field.data_dims, S, cellrange)
end

"initialize `field::Field` to `init_value` (`:initial_value` or `:norm_value`)
from Variable `attribv` attributes, over region defined by `cellrange`.
Optionally calculate transformed initial values from supplied `convertfn` and  `convertvalues`"
function init_field!(
    field::Field{D, S, V, N, M}, init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange, info,
) where {D, S, V, N, M}
    init_values!(field.values, D, field.data_dims, S, init_value, attribv, convertfn, convertvalues, cellrange, info)
end


"calculate number of degrees-of-freedom for `field::Field` over region defined by `cellrange`"
function dof_field(field::Field{D, S, V, N, M}, cellrange) where {D, S, V, N, M}
    return dof_values(D, field.data_dims, S, field.mesh, cellrange)
end

"copy `src::Field`` to `dest::Vector`, optionally restricting to region defined by `cellrange`"
function Base.copyto!(dest, doff, src::Field{D, S, V, N, M}, cellrange) where {D, S, V, N, M}
    return copyfieldto!(dest, doff, src.values, D, src.data_dims, S, cellrange)
end

"copy `src::Vector`` to `dest::Field`, optionally restricting to region defined by `cellrange`"
function Base.copyto!(dest::Field{D, S, V, N, M}, cellrange, src, soff) where {D, S, V, N, M}
    return copytofield!(dest.values, D, dest.data_dims, S, cellrange, src, soff)
end

"Calculate `dest::Field = dest::Field + a * src::Field`, optionally restricting to region defined by `cellrange`"
function add_field!(dest::Field{D, S, V, N, M}, a, cellrange, src::Field{D, S, V, N, M}) where {D, S, V, N, M}
    return add_field!(dest.values, D, dest.data_dims, S, a, cellrange, src.values)
end

function add_field_vec!(dest::Field{D, S, V, N, M}, a, cellrange, srcvalues::AbstractVector, soff) where {D, S, V, N, M}
    return add_field_vec!(dest.values, D, dest.data_dims, S, a, cellrange, srcvalues, soff)
end

"sanitized version of `values`, suitable for storing as output"
function get_values_output(field::Field{D, S, V, N, M}) where {D, S, V, N, M}
    return get_values_output(field.values, D, field.data_dims, S, field.mesh)
end


"get values from `linkvar_field`, optionally applying view defined by `linksubdomain`"
function create_accessor(
    output_data::Union{Type{D}, Type{UndefinedData}},
    # TODO - check space, data_dims,
    linkvar_field::Field{D, S, V, N, M}, linksubdomain::Union{Nothing, AbstractSubdomain};
    forceview, components,
) where {D, S, V, N, M}

    #=
    if !(outputcomponents in (VC_Single, VC_ValueOnly, VC_Undefined))
        @warn "create_accessor: invalid combination of outputcomponents=$outputcomponents and Vector data"
        return nothing
    end
    =#

    # create accessor
    if isnothing(linksubdomain)
        if forceview
            if components
                var_accessor = [view(vc, 1:length(vc)) for vc in get_components(linkvar_field.values, D)]
            else
                var_accessor = view(linkvar_field.values, 1:length(linkvar_field.values))
            end
        else
            if components
                var_accessor = get_components(linkvar_field.values, D)
            else
                var_accessor = linkvar_field.values
            end
        end
    else              
        if components
            var_accessor = [Grids.subdomain_view(vc, linksubdomain) for vc in get_components(linkvar_field.values, D)]
        else
            var_accessor = Grids.subdomain_view(linkvar_field.values, linksubdomain)
        end
    end

    return var_accessor
end

"""
    get_array(f::Field; [attributes=nothing], [selectargs...]) -> FieldArray

Return a [`FieldArray`](@ref) containing `f::Field` data values and
any attached coordinates, for the spatial region defined by `selectargs`.

Available `selectargs` depend on the grid `f.mesh`, and 
are passed to `PB.Grids.get_region`.

`attributes` (if present) are added to `FieldArray`
"""
function get_array(
    f::Field{D, ScalarSpace, V, N, M};
    attributes=nothing,
) where {D, V, N, M}

    return FieldArray("", f.values, f.data_dims, attributes)
end

function get_array(
    f::Field{D, CellSpace, V, 0, M};
    attributes=nothing,
    selectargs...
) where {D, V, M}

    values, dims = Grids.get_region(f.mesh, f.values; selectargs...)

    return FieldArray("", values, dims, attributes)
end

# single data dimension
# TODO generalize this to arbitrary data dimensions
function get_array(
    f::Field{D, CellSpace, V, 1, M};
    attributes=nothing,
    selectargs...
) where {D, V, M}

    f_space_dims_colons = ntuple(i->Colon(), ndims(f.values) - 1)
    f_size_datadim = size(f.values)[end]

    dvalues, dims = Grids.get_region(f.mesh, f.values[f_space_dims_colons..., 1]; selectargs...)
   
    d = (size(dvalues)..., f_size_datadim)
    values = Array{eltype(dvalues), length(d)}(undef, d...)
  
    if length(d) == 1
        # single cell - space dimension squeezed out
        for i in 1:f_size_datadim
            values[i], dims = Grids.get_region(f.mesh, f.values[f_space_dims_colons..., i]; selectargs...)
        end
    else
        dvalues_colons = ntuple(i->Colon(), ndims(dvalues))
        for i in 1:f_size_datadim
            dvalues, dims = Grids.get_region(f.mesh, f.values[f_space_dims_colons..., i]; selectargs...)
            values[dvalues_colons..., i] .= dvalues
        end
    end

    return FieldArray("", values, (dims..., f.data_dims...), attributes)
end

###########################################################################
# FieldRecord
##########################################################################

"""
    FieldRecord{D <: AbstractData, S <: AbstractSpace, ...}

A series of records each containing a [`Field`](@ref).

# Implementation
Fields with array values are stored in `records` as a Vector of arrays.
Fields with single values (field_single_element true) are stored as a Vector of eltype(Field.values). 
"""
struct FieldRecord{D <: AbstractData, S <: AbstractSpace, V, N, M, R}
    records::Vector{R}
    data_dims::NTuple{N, NamedDimension}
    mesh::M
    attributes::Dict{Symbol, Any}
    coords_record::Vector{FixedCoord} # coordinates attached to record dimension
end

function Base.show(io::IO, fr::FieldRecord)
    print(io, 
        "FieldRecord(eltype=", eltype(fr),", length=", length(fr), 
        ", attributes=", fr.attributes, ", coords_record=", fr.coords_record, ")"
    )
end

function Base.show(io::IO, ::MIME"text/plain", fr::FieldRecord)
    println(io, "FieldRecord(eltype=", eltype(fr),", length=", length(fr), ")") 
    println(io, "  data_dims: ", fr.data_dims)
    println(io, "  mesh: ", fr.mesh)
    println(io, "  attributes: ", fr.attributes)
    println(io, "  coords_record: ", fr.coords_record)
    return nothing
end

"test whether Field contains single elements"
function field_single_element(::Type{field_data}, N, ::Type{S}, ::Type{M}) where {field_data <: AbstractData, S <: AbstractSpace, M}
    if field_single_element(field_data, N) && (S == ScalarSpace || (S == CellSpace && M == Nothing))
        return true
    else
        return false
    end
end

field_single_element(::Type{Field{D, S, V, N, M}}) where {D, S, V, N, M} = field_single_element(D, N, S, M)
field_single_element(::Type{FR}) where {FR <: FieldRecord} = field_single_element(eltype(FR))
field_single_element(f::T) where {T} = field_single_element(T)

"create empty FieldRecord"
function FieldRecord(
    f::Field{D, S, V, N, M}, attributes;
    coords_record, 
    sizehint::Union{Nothing, Int}=nothing
) where {D, S, V, N, M}
    if field_single_element(f)
        # if Field contains single elements, store as a Vector of elements
        records = Vector{eltype(f.values)}()
    else
        # if Field contains something else, store as a Vector of those things
        records = Vector{typeof(f.values)}()
    end
    if !isnothing(sizehint)
        sizehint!(records, sizehint)
    end
    return FieldRecord{D, S, V, N, M, eltype(records)}(records, f.data_dims, f.mesh, attributes, coords_record)
end

"create a new FieldRecord, containing supplied `existing_values` data arrays"
function wrap_fieldrecord(
    existing_values, 
    field_data::Type, 
    data_dims::NTuple{N, NamedDimension},
    data_type::Union{DataType, Missing},
    space::Type{<:AbstractSpace}, 
    mesh::M,
    attributes;
    coords_record
) where {N, M}
    # check_values(
    #    existing_values, field_data, data_dims, data_type, space, spatial_size(space, mesh), 
    # )
    if field_single_element(field_data, N, space, M)
        # assume existing_values is a Vector, with element a Field value, to be stored in Field as a length-1 Vector
        # eltype(existing_values) <: Real || 
        #    throw(ArgumentError("unexpected eltype(existing_values) $(eltype(existing_values)) for ScalarSpace"))
        V = Array{eltype(existing_values), 1}
    else
        # assume existing_values is a Vector of Field values
        V = eltype(existing_values)
    end

    return FieldRecord{field_data, space, V, N, typeof(mesh), eltype(existing_values)}(
        existing_values, data_dims, mesh, attributes, coords_record
    )
end

function Base.push!(fr::FieldRecord{D, S, V, N, M, R}, f::Field{D, S, V, N, M}) where {D, S, V, N, M, R}
    if field_single_element(fr)
        # if Field contains single elements, store as a Vector of elements
        push!(fr.records, f.values[])
    else
        # if Field contains something else, store as a Vector of those things
        push!(fr.records, copy(f.values))
    end
    return fr
end

Base.length(fr::FieldRecord) = length(fr.records)

Base.eltype(::Type{FieldRecord{D, S, V, N, M, R}}) where {D, S, V, N, M, R} = Field{D, S, V, N, M}

function Base.getindex(fr::FieldRecord{D, S, V, N, M, R}, i::Int) where {D, S, V, N, M, R}

    if field_single_element(fr)
        # if Field contains single elements, FieldRecord stores as a Vector of elements
        return wrap_field([fr.records[i]], D, fr.data_dims, missing, S, fr.mesh)
    else
        # if Field contains something else, FieldRecord stores as a Vector of those things
        return wrap_field(fr.records[i], D, fr.data_dims, missing, S, fr.mesh)       
    end
end

Base.lastindex(fr::FieldRecord) = lastindex(fr.records)

function Base.copy(fr::FieldRecord{D, S, V, N, M, R}) where {D, S, V, N, M, R}
    return FieldRecord{D, S, V, N, M, R}(
        deepcopy(fr.records),
        fr.data_dims,
        fr.mesh,
        deepcopy(fr.attributes),
        copy(fr.coords_record),
    )
end

"""
    get_array(fr::FieldRecord; [recordarg] [,selectargs...]) -> FieldArray

Return a [`FieldArray`](@ref) containing `fr::FieldRecord` data values and
any attached coordinates, for records defined by `recordarg` and the 
spatial region defined by `selectargs`.

`recordarg` can be one of:
- `records=r::Int` to select a single record, or `records = first:last` to select a range.
- `<record coord>`: (where eg <record coord>=tmodel), `<record_coord>=t::Float64 to select a single record
  with the nearest value of `fr.coords_record`, or `<record_coord>=(first::Float64, last::Float64)` (a Tuple) to select a range
  starting at the record with the nearest value of `fr.coords_record` before `first` and ending at the nearest record after
  `last`.

Available `selectargs` depend on the grid `fr.mesh`, and 
are passed to `PB.Grids.get_region`.
"""
function get_array(
    fr::FieldRecord{D, S, V, N, M, R}; selectargs...
) where {D, S, V, N, M, R}

    # select records to use and create NamedDimension
    ridx = 1:length(fr)
    selectargs_region = Dict()
    selectargs_records = NamedTuple()
    for (k, v) in selectargs
        if k==:records
            if v isa Integer
                ridx = [v]                
            else
                ridx = v
            end
            selectargs_records = (records=v,)
        elseif String(k) in getfield.(fr.coords_record, :name)
            # find ridx corresponding to a coordinate
            for cr in fr.coords_record
                if String(k) == cr.name
                    ridx, cvalue = find_indices(cr.values, v)
                    selectargs_records=NamedTuple((k=>cvalue,))
                end
            end            
        else
            selectargs_region[k] = v
        end
    end
    records_dim = NamedDimension("records", length(ridx), get_region(fr.coords_record, ridx))

    # add attributes for selection used
    attributes = copy(fr.attributes)
    attributes[:filter_records] = selectargs_records
    attributes[:filter_region] = NamedTuple(selectargs_region)

    # Generate name from attributes
    name = get(attributes, :domain_name, "")
    name *= isempty(name) ? "" : "."
    name *= get(attributes, :var_name, "")
    if !isempty(selectargs_region) || !isempty(selectargs_records)
        name *= "(" * join(["$k=$v" for (k, v) in merge(Dict(pairs(selectargs_records)), selectargs_region)], ", ") * ")"
    end

    # Select selectargs_region 
    if field_single_element(fr)
        # create FieldArray directly from FieldRecord
        isempty(selectargs_region) || 
            throw(ArgumentError("invalid index $selectargs_region"))
        if length(ridx) > 1
            return FieldArray(
                name,
                fr.records[ridx], 
                (fr.data_dims..., records_dim),
                attributes
            )
        else
            # squeeze out records dimension
            return FieldArray(
                name,
                fr.records[first(ridx)], 
                fr.data_dims,
                attributes
            )
        end
    else
        # pass through to Field

        # get FieldArray from first Field record and use this to figure out array shapes etc
        far = get_array(fr[first(ridx)]; selectargs_region...)        
        if length(ridx) > 1
            # add additional (last) dimension for records
            favalues = Array{eltype(far.values), length(far.dims)+1}(undef, size(far.values)..., length(ridx))
            fa = FieldArray(
                name,
                favalues, 
                (far.dims..., records_dim), 
                attributes,
            )
            # copy values for first record
            if isempty(far.dims)
                fa.values[1] = far.values
            else
                fa.values[fill(Colon(), length(far.dims))..., 1] .= far.values
            end
        else
            # squeeze out record dimension so this is just fa with attributes added
            fa = FieldArray(
                name, 
                far.values,
                far.dims,
                attributes
            ) 
        end
        
        # fill with values from FieldArrays for Fields for remaining records
        if length(ridx) > 1
            for (i, r) in enumerate(ridx[2:end])
                far = get_array(fr[r]; selectargs_region...)
                if isempty(far.dims)
                    fa.values[i+1] = far.values
                else
                    fa.values[fill(Colon(), length(far.dims))..., i+1] .= far.values
                end
            end
        end

        return fa
    end

end

