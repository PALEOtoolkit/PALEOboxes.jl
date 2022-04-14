###############################################
# ArrayScalarData
################################################

"""
    ArrayScalarData

An Array of Field scalars (eg a intensity on a wavelength grid).

NB: only implements minimal methods sufficient to allow use as a model internal variable, not as a state variable.
"""
struct ArrayScalarData <: AbstractData
end


function allocate_values(
    field_data::Type{ArrayScalarData}, data_dims::Tuple{NamedDimension, Vararg{NamedDimension}}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}};
    thread_safe, allocatenans,
)
    !thread_safe || throw(ArgumentError("ArrayScalarData with thread_safe==true"))

    s = _size_arrayscalardata(data_dims, space, spatial_size)
    d = Array{data_type, length(s)}(undef, s...)

    if allocatenans
        fill!(d, NaN)
    end
       
    return d
end

function check_values(
    existing_values, field_data::Type{ArrayScalarData}, data_dims::Tuple{NamedDimension, Vararg{NamedDimension}}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}},
) where {N}

    check_data_type(existing_values, data_type)

    existing_size = size(existing_values)
    req_size = _size_arrayscalardata(data_dims, space, spatial_size)
    existing_size == req_size ||
        throw(ArgumentError("size mismatch: supplied $existing_size, require $req_size"))
    return nothing
end


function zero_values!(values, field_data::Type{ArrayScalarData}, data_dims::Tuple{NamedDimension, Vararg{NamedDimension}}, space::Type{ScalarSpace}, cellrange) where {N}
    values .= 0.0
    return nothing
end

function zero_values!(values, field_data::Type{ArrayScalarData}, data_dims::Tuple{NamedDimension, Vararg{NamedDimension}}, space::Type{CellSpace}, cellrange)  where {N}
    data_dims_colons = ntuple(i->Colon(), length(data_dims))
    @inbounds for i in cellrange.indices
        values[i, data_dims_colons...] .= 0.0
    end
    return nothing
end

# NB: Tuple{Int, Vararg{Int}} is an NTuple with at least one element
function _size_arrayscalardata(data_dims::Tuple{NamedDimension, Vararg{NamedDimension}}, space, spatial_size::Tuple{Integer, Vararg{Integer}})
    size_data_dims = Tuple(d.size for d in data_dims)
    if space === ScalarSpace
        # Don't include a first dimension length 1 for ScalarSpace
        return size_data_dims
    else
        return (spatial_size..., size_data_dims...)
    end
end
