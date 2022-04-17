
################################################################
# ScalarData
#################################################################


"""
    ScalarData <: AbstractData

A Field scalar (eg a biogeochemical concentration)
"""
struct ScalarData <: AbstractData
end

num_components(field_data::Type{ScalarData})         = 1

get_components(values, field_data::Type{ScalarData}) = [values]

function allocate_values(
    field_data::Type{ScalarData}, data_dims::Tuple{}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}};
    thread_safe, allocatenans,
)

    if thread_safe
        spatial_size == (1, ) ||
            throw(ArgumentError("thread_safe with spatial_size != (1, )"))
        d = AtomicScalar{data_type}()
    else
        d = Array{data_type, length(spatial_size)}(undef, spatial_size...)
    end

    if allocatenans
        fill!(d, NaN)
    end
    
    return d
end

function get_values_output(values::AtomicScalar, data_type::Type{ScalarData}, data_dims::Tuple{}, space, mesh)
    return [values[]]
end

function check_values(
    existing_values, field_data::Type{ScalarData}, data_dims::Tuple{}, data_type, space::Type{<:AbstractSpace}, spatial_size::Tuple{Integer, Vararg{Integer}}
)
    check_data_type(existing_values, data_type)

    existing_size = size(existing_values) 
    existing_size == spatial_size ||
        throw(ArgumentError("size mismatch: supplied $existing_size, require spatial_size=$spatial_size"))

    return nothing
end

function init_values!(
    values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{S},
    init_value::Symbol, attribv::VariableBase, convertfn, convertvalues, cellrange,
    info::NTuple{3, String}
) where {S <: Union{ScalarSpace, CellSpace}}
    varinfo, convertinfo, trsfrinfo = info

    initial_value = get_attribute(attribv, init_value)

    @info "init_values!     :$(rpad(init_value, 20)) $(rpad(varinfo, 30)) "*
                "= $initial_value$convertinfo$trsfrinfo"
    
    if S === ScalarSpace
        values[] = initial_value*convertfn(convertvalues, nothing)
    else
        if initial_value isa Real # ie not a Vector
            for i in cellrange.indices
                values[i] = initial_value*convertfn(convertvalues,  i)
            end
        elseif length(initial_value) == length(values)
            for i in cellrange.indices
                values[i] = initial_value[i]*convertfn(convertvalues,  i)
            end
        else
            error("init_values!: configuration error, Vector Variable $varinfo and Vector :initial_value lengths mismatch ($(length(values)) != $(length(initial_value)))")
        end
    end
    
    return nothing
end

function zero_values!(values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{ScalarSpace}, cellrange)
    values[] = 0.0
    return nothing
end

function zero_values!(values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange)
    @inbounds for i in cellrange.indices
        values[i] = 0.0
    end
    return nothing
end

dof_values(field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{ScalarSpace}, mesh, cellrange) = 1
dof_values(field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, mesh, cellrange) = length(cellrange.indices)
dof_values(field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, mesh, cellrange::Nothing) = mesh.ncells

function copyfieldto!(dest, doff, values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{ScalarSpace}, cellrange)
    dest[doff] = values[]
    return 1
end

function copyfieldto!(dest, doff, values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange)
    @inbounds for i in cellrange.indices
        dest[doff] = values[i]
        doff += 1
    end
    return length(cellrange.indices)
end

function copyfieldto!(dest, doff, values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange::Nothing)
    @inbounds copyto!(dest, doff, values, 1, length(values))
    return length(values)
end

function copytofield!(values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{ScalarSpace}, cellrange, src, soff)
    values[] = src[soff]
    return 1
end

function copytofield!(values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange, src, soff)
    @inbounds for i in cellrange.indices
        values[i] = src[soff]
        soff += 1
    end
    return length(cellrange.indices)
end

function copytofield!(values, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange::Nothing, src, soff)
    @inbounds copyto!(values, 1, src, soff, length(values))
    return length(values)
end

function add_field!(dest, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{ScalarSpace}, a, cellrange, src)
    dest[] += a*src[]
    return nothing
end

function add_field!(dest, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, a, cellrange, src)
    @inbounds for i in cellrange.indices
        dest[i] += a*src[i]
    end
    return nothing
end

function add_field!(dest, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, a, cellrange::Nothing, src)
    dest .+= a .* src
    return nothing
end

function add_field_vec!(dest, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{ScalarSpace}, a, cellrange, src, soff)
    dest[] += a*src[soff]
    return 1
end

function add_field_vec!(dest, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, a, cellrange, src, soff)
    @inbounds for i in cellrange.indices
        dest[i] += a*src[soff]
        soff += 1
    end
    return length(cellrange.indices)
end

add_field_vec!(dest, field_data::Type{ScalarData}, data_dims::Tuple{}, space::Type{CellSpace}, cellrange::Nothing, src, soff) =
    add_field_vec!(dest, field_data, data_dims, space, a, (indices=1:length(dest),), src, soff)


