
################################
# Coordinates
#################################

"""
    FixedCoord

A fixed (state independent) coordinate
"""
mutable struct FixedCoord
    name::String
    values::Vector{Float64}
    attributes::Dict{Symbol, Any}
end

"""
    append_units(name::AbstractString, attributes) -> "name (units)"

Utility function to append variable units string to a variable name for display.
"""
function append_units(name::AbstractString, attributes::Dict{Symbol, Any})
    units = get(attributes, :units, "")
    if isempty(units)
        return name
    else
        return name*" ($units)"
    end
end

append_units(name::AbstractString, attributes::Nothing) = name

"""
    build_coords_edges(coords_vec::Vector{FixedCoord}) -> Vector{Float64}

Build a vector of coordinate edges (length `n+1``) from `coords_vec`, assuming the PALEO
convention that `coords_vec` contains three elements with 
cell midpoints, lower edges, upper edges each of length `n`, in that order.

Falls back to just returning the first entry in `coords_vec` for other cases.
"""
function build_coords_edges(coords_vec::Vector{FixedCoord})
    
    if length(coords_vec) == 1 || length(coords_vec) > 3
        # 1 coordinate or something we don't understand - take first
        co = first(coords_vec)
        co_values = co.values
        co_label = append_units(co.name, co.attributes)
    elseif length(coords_vec) in (2, 3)
        # 2 coordinates assume lower, upper edges
        # 3 coordinates assume mid, lower, upper
        co_lower = coords_vec[end-1]
        co_upper = coords_vec[end]
        co_label = append_units(co_lower.name*", "*co_upper.name, co_lower.attributes)
        first(co_lower.values) < first(co_upper.values) ||
            @warn "build_coords_edges: $co_label co_lower is > co_upper - check model grid"            
        if co_lower.values[end] > co_lower.values[1] # ascending order
            co_lower.values[2:end] == co_upper.values[1:end-1] ||
                @warn "build_coords_edges: $co_label lower and upper edges don't match"
            co_values = [co_lower.values; co_upper.values[end]]
        else # descending order
            co_lower.values[1:end-1] == co_upper.values[2:end] ||
                @warn "build_coords_edges: $co_label lower and upper edges don't match"
            co_values = [co_upper.values[1]; co_lower.values]
        end

    end

    return co_values, co_label
end

"guess coordinate edges from midpoints, assuming uniform spacing"
function guess_coords_edges(x_midpoints)
    first_x = x_midpoints[1] - 0.5*(x_midpoints[2] - x_midpoints[1])
    last_x = x_midpoints[end] + 0.5*(x_midpoints[end] - x_midpoints[end-1])
    return [first_x; 0.5.*(x_midpoints[1:end-1] .+ x_midpoints[2:end]); last_x]
end


function get_region(fc::FixedCoord, indices::AbstractVector)
    return FixedCoord(fc.name, fc.values[indices], fc.attributes)
end

function get_region(fcv::Vector{FixedCoord}, indices::AbstractVector)
    return [FixedCoord(fc.name, fc.values[indices], fc.attributes) for fc in fcv]
end


"find indices of coord from first before range[1] to first after range[2]"
function find_indices(coord::AbstractVector, range)
    length(range) == 2 ||
        throw(ArgumentError("find_indices: length(range) != 2  $range"))

    idxstart = findlast(t -> t<=range[1], coord)
    isnothing(idxstart) && (idxstart = 1)

    idxend = findfirst(t -> t>=range[2], coord)
    isnothing(idxend) && (idxend = length(coord))

    return idxstart:idxend, (coord[idxstart], coord[idxend])
end

"find indices of coord nearest val"
function find_indices(coord::AbstractVector, val::Real)
    idx = 1
    for i in 1:length(coord)
        if abs(coord[i] - val) < abs(coord[idx] - val)
            idx = i
        end
    end

    return [idx], coord[idx]
end

#################################################
# Dimensions
#####################################################

"""
    NamedDimension

A named dimension, with optional preferred coordinates `coords`

PALEO convention is that where possible `coords` contains three elements, for cell
midpoints, lower edges, upper edges, in that order.
"""
mutable struct NamedDimension
    name::String
    size::Int64    
    coords::Vector{String} # may be empty
end

if false
"create from size only (no coords)"
function NamedDimension(name, size::Integer, coords=String[])
    return NamedDimension(
        name, 
        size, 
        String[],
    )
end

"create from coord mid-points"
function NamedDimension(name, coord::AbstractVector)
    return NamedDimension(
        name, 
        length(coord), 
        [
            FixedCoord(name, coord, Dict{Symbol, Any}()),
        ]
    )
end

"create from coord mid-points and edges"
function NamedDimension(name, coord::AbstractVector, coord_edges::AbstractVector)
    if coord[end] > coord[1]
        # ascending order
        coord_lower = coord_edges[1:end-1]
        coord_upper = coord_edges[2:end]
    else
        # descending order
        coord_lower = coord_edges[2:end]
        coord_upper = coord_edges[1:end-1]
    end
    return NamedDimension(
        name, 
        length(coord), 
        [
            FixedCoord(name, coord, Dict{Symbol, Any}()),
            FixedCoord(name*"_lower", coord_lower, Dict{Symbol, Any}()),
            FixedCoord(name*"_upper", coord_upper, Dict{Symbol, Any}()),
        ]
    )
end

function get_region(nd::NamedDimension, indices::AbstractVector)
    return NamedDimension(nd.name, length(indices), get_region(nd.coords, indices))
end

"""
    build_coords_edges(nd::NamedDimension) -> Vector{Float64}

Call [`build_coords_edges`](@ref)(nd.coords), or fallback to just returning indices
if no coords present.
"""
function build_coords_edges(nd::NamedDimension)
    if !isempty(nd.coords)
        return build_coords_edges(nd.coords)
    else
        @warn "no coords for NamedDimension $(nd.name), returning indices"
        return collect(1:nd.size), nd.name*" (indices)"
    end
end

end

function Base.show(io::IO, nd::NamedDimension)
    print(io, "NamedDimension(name=", nd.name, ", size=", nd.size, ", coords=", nd.coords, ")")
    return nothing
end
