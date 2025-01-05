module Grids

import NCDatasets
import PALEOboxes as PB

import Infiltrator # Julia debugger

###########################
# Subdomains
###########################
"""
    AbstractSubdomain

Defines the relationship between two [`PB.Domain`](@ref)s by mapping indices in one Domain to related indices in the other,
eg interior cells adjacent to a boundary.

Concrete subtypes should implement:

[`subdomain_view`](@ref)

[`subdomain_indices`](@ref)
"""
PB.AbstractSubdomain

"""
    BoundarySubdomain <: PB.AbstractSubdomain

A 2D subdomain corresponding to the 2D boundary Domain associated with a 3D interior Domain:
- `indices[ibnd]` is the index of the 3D interior cell corresponding to a 2D boundary cell `ibnd`.
"""
struct BoundarySubdomain <: PB.AbstractSubdomain
    indices::Vector{Int}
end

"""
    subdomain_view(values, subdomain::BoundarySubdomain) -> view

Create a `view` on `values` in an interior `Domain` to access cells corresponding to indices
in a boundary `Domain`.
"""
function subdomain_view(values, subdomain::BoundarySubdomain)
    return view(values, subdomain.indices)
end

"""
    subdomain_indices(subdomain::BoundarySubdomain) -> nothing

No additional indices required to access Variables in an interior `Domain` from a boundary `Domain`
(view created by [`subdomain_view`](@ref) is sufficient).
"""
function subdomain_indices(subdomain::BoundarySubdomain)
    return nothing
end

"""
    InteriorSubdomain <: PB.AbstractSubdomain

A 3D subdomain corresponding to the 3D interior Domain associated with a 2D boundary Domain: `indices[iint]` is either:
- `missing` if `iint` is the index of an interior cell in the 3D Domain, or
- the index of the 2D boundary Domain cell corresponding to the boundary-adjacent 3D interior Domain cell `iint`
"""
struct InteriorSubdomain <: PB.AbstractSubdomain
    indices::Vector{Union{Missing, Int}}
end

"Create InteriorSubdomain given size of interior Domain and 
 boundaryindices[ibnd] index of interior point corresponding to boundary point ibnd"
function InteriorSubdomain(ninterior, boundaryindices)
    indices = Vector{Union{Missing, eltype(boundaryindices)}}(undef, ninterior)
    for i in 1:length(boundaryindices)
        indices[boundaryindices[i]] = i 
    end
    return InteriorSubdomain(indices)
end


"""
    subdomain_view(values, subdomain::InteriorSubdomain) -> var

Return unmodified `values` from a boundary `Domain`, `indices` to access 
from interior supplied by [`subdomain_indices`](@ref)
"""
function subdomain_view(values, subdomain::InteriorSubdomain)
    return values
end

"""
    subdomain_indices(subdomain::InteriorSubdomain) -> indices

Return `indices` to access Variables in a boundary `Domain` from interior `Domain`
(will have `missing` entries where interior cells do not correspond to boundary)
"""
function subdomain_indices(subdomain::InteriorSubdomain)
    return subdomain.indices
end

function is_boundary(subdomain::InteriorSubdomain, i)
    return !ismissing(subdomain.indices[i])
end

"""
    subdomain_view(values, subdomain::Nothing) -> values

Fallback when `subdomain == nothing`
"""
function subdomain_view(values, subdomain::Nothing)
    return values
end

"""
    subdomain_indices(subdomain::Nothing) -> nothing

fallback when `subdomain == nothing`
"""
function subdomain_indices(subdomain::Nothing)
    return nothing
end



#########################################
# Grids
##########################################

"""
    AbstractMesh

Defines additional geometric and topological information for [`PB.Domain`](@ref)

Concrete subtypes should implement methods:

[`available_spaces`](@ref) tuple of [`PB.AbstractSpace`](@ref) types supported. 

[`PB.has_internal_cartesian`](@ref) `true` if subtype uses an internal spatial representation that differs from 
the external (cartesian) representation.

[`PB.get_dimensions`](@ref)

[`PB.internal_size`](@ref), optionally [`PB.cartesian_size`](@ref)

[`PB.Grids.set_subdomain!`](@ref), [`PB.Grids.get_subdomain`](@ref)

[`PB.Grids.create_default_cellrange`](@ref)

# Internal and cartesian representations
If a subtype uses an internal representation that differs from the external (cartesian) representation,
it should define [`PB.has_internal_cartesian`](@ref) `true` and implement [`cartesian_to_internal`](@ref),
[`internal_to_cartesian`](@ref), and [`PB.cartesian_size`](@ref).
"""
PB.AbstractMesh

"parse eg \"CartesianLinearGrid\" as CartesianLinearGrid"
function Base.parse(::Type{PB.AbstractMesh}, str::AbstractString)   
    dtype = getproperty(@__MODULE__, Symbol(str))

    dtype <: PB.AbstractMesh || 
        throw(ArgumentError("$str is not a subtype of AbstractMesh"))
    return dtype
end

"""
    create_default_cellrange(domain, grid, [; operatorID=0]) -> CellRange

Create a CellRange for entire `domain` and supplied `operatorID`
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::PB.AbstractMeshOrNothing) end


"""
    set_subdomain!(grid::PB.AbstractMesh, subdomainname::AbstractString, subdom::PB.AbstractSubdomain, allowcreate::Bool=false)

Set Subdomain
"""
function set_subdomain!(grid::PB.AbstractMesh, subdomainname::AbstractString, subdom::PB.AbstractSubdomain, allowcreate::Bool=false)
    if !haskey(grid.subdomains, subdomainname) && !allowcreate
        error("attempt to create new subdomain name='$subdomainname' failed (allowcreate=false)")
    end
    grid.subdomains[subdomainname] = subdom
    return nothing
end

"""
    get_subdomain(grid::PB.AbstractMesh, subdomainname::AbstractString) -> PB.AbstractSubdomain

Get Subdomain
"""
function get_subdomain(grid::PB.AbstractMesh, subdomainname::AbstractString)
    subdomain = get(grid.subdomains, subdomainname, nothing)
    !isnothing(subdomain) || error("get_subdomain: no subdomain $subdomainname")
    
    return subdomain
end

"""
    available_spaces(grid::PB.AbstractMesh) -> NTuple{PB.Space}

Tuple of Spaces supported by this grid
"""
function available_spaces end


# generic handler when subdomainname present
function PB.internal_size(space::Type{<:PB.AbstractSpace}, grid::PB.AbstractMeshOrNothing, subdomainname::AbstractString)
    if isempty(subdomainname)
        return PB.internal_size(space, grid)
    else
        if space === PB.CellSpace
            subdomain = get_subdomain(grid, subdomainname)
            return (length(subdomain.indices),)
        else
            error("internal_size: cannot specify subdomain with space=$space (subdomainname=$subdomainname)")
        end
    end
end

# generic handler when subdomain present
function PB.get_dimensions(
    grid::PB.AbstractMeshOrNothing, space::Type{<:PB.AbstractSpace}, subdomainname::AbstractString,
    expand_cartesian=false,
)
    if isempty(subdomainname)
        return PB.get_dimensions(grid, space; expand_cartesian)
    else
        if space === PB.CellSpace
            subdomain = get_subdomain(grid, subdomainname)
            return (PB.NamedDimension("subdomain_"*subdomainname, length(subdomain.indices)),)
        else
            error("get_dimensions: cannot specify subdomain with space=$space (subdomainname=$subdomainname)")
        end
    end
end

# scalar space always supported
PB.internal_size(space::Type{PB.ScalarSpace}, grid::PB.AbstractMeshOrNothing) = (1, )
PB.get_dimensions(grid::PB.AbstractMeshOrNothing, space::Type{PB.ScalarSpace}; expand_cartesian=false) = ()
# ignore subdomain 
PB.internal_size(space::Type{PB.ScalarSpace}, grid::PB.AbstractMeshOrNothing, subdomainname::AbstractString) = (1, )

# fallbacks for cartesian 
PB.cartesian_size(space, grid::PB.AbstractMeshOrNothing) = PB.internal_size(space, grid)
PB.has_internal_cartesian(grid::PB.AbstractMeshOrNothing, Space::Type{<:PB.AbstractSpace}) = false

"""
    cartesian_to_internal(grid::PB.AbstractMeshOrNothing, griddata::AbstractArray)

Map a cartesian (external) array to internal representation
"""
cartesian_to_internal(grid::PB.AbstractMeshOrNothing, griddata::AbstractArray) = griddata

"""
    cartesian_to_internal(grid::PB.AbstractMeshOrNothing, griddata::AbstractArray)

Map an internal array to cartesian (external) representation
"""
internal_to_cartesian(grid::PB.AbstractMeshOrNothing, griddata::AbstractArray) = griddata

# generic get all dimensions
function PB.get_dimensions(grid::PB.AbstractMesh)
    dims = PB.NamedDimension[]
    for space in available_spaces(grid)
        append!(dims, PB.get_dimensions(grid, space))
        if PB.has_internal_cartesian(grid, space)
            append!(dims, PB.get_dimensions(grid, space; expand_cartesian=true))
        end
    end
    for sdn in keys(grid.subdomains)
        append!(dims, PB.get_dimensions(grid, PB.CellSpace, sdn))
    end
    append!(dims, _extra_dimensions(grid))
    return unique(d->d.name, dims)
end

# implement if subtype has defines additional dimensions not included in list by Space and subdomain
_extra_dimensions(grid::PB.AbstractMesh) = ()

# generic coordinate set / get 
function PB.set_coordinates!(grid::PB.AbstractMesh, dimname::AbstractString, coordinates::Vector{<:AbstractString})
    dims = PB.get_dimensions(grid)
    idx = findfirst(nd -> nd.name==dimname, dims)
    !isnothing(idx) || error("grid does not contain dimension $dimname")
    delete!(grid.coordinates, dimname)
    if !isempty(coordinates)
        grid.coordinates[dimname] = coordinates
    end
end

function PB.get_coordinates(grid::PB.AbstractMesh, dimname::AbstractString)    
    return get(grid.coordinates, dimname, String[])
end

# helper function
function _show_subdomains_dimensions_coordinates(io, grid::PB.AbstractMesh)
    println(io, "  spaces: ", available_spaces(grid))
    for space in available_spaces(grid)
        println(io, "      ", space)
        if PB.has_internal_cartesian(grid, space)
            println(io, "        dimensions (internal): ", PB.get_dimensions(grid, space; expand_cartesian=false))
            println(io, "        dimensions (cartesian): ", PB.get_dimensions(grid, space; expand_cartesian=true))
        else
            println(io, "        dimensions: ", PB.get_dimensions(grid, space))
        end
    end
    println(io, "  subdomains: ", keys(grid.subdomains))
    for subdomain in keys(grid.subdomains)
        println(io, "        dimensions: ", PB.get_dimensions(grid, PB.CellSpace, subdomain))
    end
    
    extra_dimensions = _extra_dimensions(grid)
    if !isempty(extra_dimensions)
        println(io, "  additional dimensions: ", extra_dimensions)
    end

    if !isempty(grid.coordinates)
        println(io, "  coordinates:")
        for (dimname, coordnames) in grid.coordinates
            println(io, "      ", dimname, " => ", coordnames)
        end
    end

    return nothing
end

# generic fallback for named cells 
substitute_cell_names(grid::PB.AbstractMeshOrNothing, cells) = cells
substitute_cell_names(grid::PB.AbstractMeshOrNothing, cells::Union{Number, Symbol}) =
    substitute_cell_names(grid, [cells])[]

# generic fallback for column indices
function column_indices(grid::PB.AbstractMeshOrNothing, column)
    throw(ArgumentError("grid $grid does not support column selection"))
end



###################################
# Fallbacks for Domain with no grid
#####################################

# allow Vector variables length 1 in a 0D Domain without a grid
available_spaces(grid::Nothing) = (PB.ScalarSpace, PB.CellSpace)
PB.internal_size(space::Type{PB.CellSpace}, grid::Nothing) = (1, )

PB.get_dimensions(grid::Nothing, space::Type{PB.CellSpace}; expand_cartesian=false) = (PB.NamedDimension("cells", 1), )

get_subdomain(grid::Nothing, subdomainname::AbstractString) = error("get_subdomain: no subdomain $subdomainname")

"""
    create_default_cellrange(domain, grid::Nothing [; operatorID=0]) -> CellRange

Create a CellRange for entire `domain`. Fallback for a domain with no grid
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::Nothing; operatorID=0)
    return PB.CellRange(domain=domain, indices=1:PB.get_length(domain), operatorID=operatorID)
end


##################################
# UnstructuredVectorGrid
##################################

"""
    UnstructuredVectorGrid <: PB.AbstractMesh

Minimal Grid for a Vector Domain, defines only some named cells for plotting
"""
Base.@kwdef mutable struct UnstructuredVectorGrid <: PB.AbstractMesh
    ncells::Int64

    "Define some named cells for plotting (only)"
    cellnames::Dict{Symbol,Int} = Dict{Symbol,Int}()

    subdomains::Dict{String, PB.AbstractSubdomain} = Dict{String, PB.AbstractSubdomain}()

    coordinates::Dict{String, Vector{String}} = Dict{String, Vector{String}}()
end

function Base.show(io::IO, grid::UnstructuredVectorGrid)
    print(io, "UnstructuredVectorGrid(ncells=", grid.ncells, 
        ", cellnames=", grid.cellnames, ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", grid::UnstructuredVectorGrid)
    println(io, "UnstructuredVectorGrid")
    println(io, "  ncells=", grid.ncells)
    println(io, "  cellnames=", grid.cellnames)
    _show_subdomains_dimensions_coordinates(io, grid)
    return nothing
end

available_spaces(grid::UnstructuredVectorGrid) = (PB.ScalarSpace, PB.CellSpace)
PB.internal_size(space::Type{PB.CellSpace}, grid::UnstructuredVectorGrid) = (grid.ncells, )

PB.get_dimensions(grid::UnstructuredVectorGrid, space::Type{PB.CellSpace}; expand_cartesian=false) = 
    (PB.NamedDimension("cells", grid.ncells), )

"""
    create_default_cellrange(domain, grid::UnstructuredVectorGrid [; operatorID=0]) -> CellRange

Create a CellRange for entire `domain` - use linear index.
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::UnstructuredVectorGrid; operatorID=0)
    return PB.CellRange(domain=domain, indices=1:grid.ncells, operatorID=operatorID)
end

substitute_cell_names(grid::UnstructuredVectorGrid, cells::Union{Number, Symbol}) =
    substitute_cell_names(grid, [cells])[]
function substitute_cell_names(grid::UnstructuredVectorGrid, cells)
    cell_indices = Int[]

    for cell in cells
        if cell isa Int
            push!(cell_indices, cell)
        else
            idx = get(grid.cellnames, cell, nothing)
            !isnothing(idx) ||
                throw(ArgumentError("cell ':$cell' not present in  grid.cellnames=$(grid.cellnames)"))
            push!(cell_indices, idx)
        end
    end

    return cell_indices
end

############################################
# UnstructuredColumnGrid
###########################################

"""
    UnstructuredColumnGrid <: PB.AbstractMesh

Minimal Grid for a Vector Domain composed of columns (not necessarily forming a 2-D array).

# Fields
- `ncells::Int` total number of cells in this Domain
- `Icolumns::Vector{Vector{Int}}`: Icolumns[n] should be the indices of column n, in order from surface to floor, where n is also the index of 
any associated boundary surface.
- `columnnames::Vector{Symbol}`: optional column names
- `coordinates::Dict{String, Vector{String}}`: optional attached coordinates
"""
Base.@kwdef mutable struct UnstructuredColumnGrid <: PB.AbstractMesh
    ncells::Int64    
    Icolumns::Vector{Vector{Int}}    

    "Define optional column names"
    columnnames::Vector{Symbol} = Symbol[]

    subdomains::Dict{String, PB.AbstractSubdomain} = Dict{String, PB.AbstractSubdomain}()

    coordinates::Dict{String, Vector{String}} = Dict{String, Vector{String}}()
end

function Base.show(io::IO, grid::UnstructuredColumnGrid)
    print(io, "UnstructuredColumnGrid(ncells=", grid.ncells, ", columns: ", length(grid.Icolumns),
        ", columnnames=", grid.columnnames, ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", grid::UnstructuredColumnGrid)
    println(io, "UnstructuredColumnGrid")
    println(io, "  ncells=", grid.ncells)
    println(io, "  columnames=", grid.columnnames)
    _show_subdomains_dimensions_coordinates(io, grid)
    return nothing
end

available_spaces(grid::UnstructuredColumnGrid) = (PB.ScalarSpace, PB.ColumnSpace, PB.CellSpace)
PB.internal_size(space::Type{PB.CellSpace}, grid::UnstructuredColumnGrid) = (grid.ncells, )
PB.internal_size(space::Type{PB.ColumnSpace}, grid::UnstructuredColumnGrid) = (length(grid.Icolumns), )

PB.get_dimensions(grid::UnstructuredColumnGrid, space::Type{PB.CellSpace}; expand_cartesian=false) = 
    (PB.NamedDimension("cells", grid.ncells), )

PB.get_dimensions(grid::UnstructuredColumnGrid, space::Type{PB.ColumnSpace}; expand_cartesian=false) = 
    (PB.NamedDimension("columns", length(grid.Icolumns)), )

"""
    create_default_cellrange(domain, grid::UnstructuredColumnGrid [; operatorID=0]) -> CellRangeColumns

Create a CellRange for entire `domain`. Return a [`PB.CellRangeColumns`](@ref) with iterators for columns and cells.
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::UnstructuredColumnGrid; operatorID=0)
    return PB.CellRangeColumns(
        domain=domain, 
        indices=1:grid.ncells, 
        columns=[
            (isurf, PB.replace_contiguous_range(grid.Icolumns[isurf])) 
            for isurf in eachindex(grid.Icolumns)
        ],
        operatorID=operatorID
    )
end

function column_indices(grid::UnstructuredColumnGrid, column)
    if column isa Int
        column in 1:length(grid.Icolumns) ||
            throw(ArgumentError("column index $column out of range"))
        colidx = column
    else
        colidx = findfirst(isequal(column), grid.columnnames)
        !isnothing(colidx) || 
            throw(ArgumentError("columnname '$column' not present in  grid.columnnames=$(grid.columnnames)"))
    end

    return grid.Icolumns[colidx]
end

######################################
# CartesianLinearGrid
#######################################

"""
    CartesianLinearGrid <: PB.AbstractMesh

nD grid with netcdf CF1.0 coordinates, using Vectors for PALEO internal representation of Variables,
with a mapping linear indices <--> some subset of grid indices.

The linear indices mapping should be set with `set_linear_index`. Conversion of Field values between
the PALEO internal representation (a Vector with a linear index) and a Cartesian Array with multiple dimensions
(for import and export of model output) is then implemented by `cartesian_to_internal` and `internal_to_cartesian`.

# Fields
- `ncells::Int64`: number of cells in Domain = `length(linear_index)` (may be subset of total points in `prod(dims)`)
- `dimensions::Vector{PB.NamedDimensions}`: names and sizes of cartesian spatial dimensions (ordered list)
- `coordinates::Dict{String, Vector{String}}`: optional attached coordinates
"""
Base.@kwdef mutable struct CartesianLinearGrid{N} <: PB.AbstractMesh
    
    ncells::Int64                           = -1 
    ncolumns::Int64                         = -1

    dimensions::Vector{PB.NamedDimension}   = Vector{PB.NamedDimension}(undef, N) # netcdf dimension names
    dimensions_extra::Vector{PB.NamedDimension} = [PB.NamedDimension("bnds", 2)] # additional dimension names

    "index of longitude dimension (if any)"
    londim::Int                            = 0
    "index of latitude dimension (if any)"
    latdim::Int                             = 0
    "index of z dimension (if any)"
    zdim::Int                               = 0
    "index of surface in z coordinate (if any) (1 or length(z dim))" 
    zidxsurface::Int                       = 0
    "multiplier to use for display (eg -1.0 to convert depth to height)"
    display_mult::Vector{Float64}          = ones(N)

    subdomains::Dict{String, PB.AbstractSubdomain}   = Dict{String, PB.AbstractSubdomain}()

    coordinates::Dict{String, Vector{String}} = Dict{String, Vector{String}}()

    "cartesian -> linear index mapping (initialized to `missing` ie linear index contains 0 cells)"
    linear_index::Array{Union{Missing, Int32}, N} = Array{Union{Missing, Int32},N}(undef, zeros(Int,N)...) 
    "linear -> cartesian index mapping (initialized to empty Vector ie linear index contains 0 cells)"
    cartesian_index::Vector{CartesianIndex{N}}  = Vector{CartesianIndex{N}}()
end

function Base.show(io::IO, grid::CartesianLinearGrid)
    print(io, "CartesianLinearGrid(ncells=", grid.ncells, ", dimensions: ", grid.dimensions, ", dimensions_extra: ", grid.dimensions_extra,
        ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", grid::CartesianLinearGrid)
    println(io, "CartesianLinearGrid")
    println(io, "  ncells=", grid.ncells)
    _show_subdomains_dimensions_coordinates(io, grid)
    return nothing
end

available_spaces(grid::CartesianLinearGrid{2}) = (PB.ScalarSpace, PB.CellSpace)
available_spaces(grid::CartesianLinearGrid{3}) = (PB.ScalarSpace, PB.ColumnSpace, PB.CellSpace)

PB.has_internal_cartesian(grid::CartesianLinearGrid, ::Type{PB.CellSpace}) = true
PB.internal_size(space::Type{PB.CellSpace}, grid::CartesianLinearGrid) = (grid.ncells, )
PB.internal_size(space::Type{PB.ColumnSpace}, grid::CartesianLinearGrid{3}) = (grid.ncolumns, )
PB.cartesian_size(space::Type{PB.CellSpace}, grid::CartesianLinearGrid) = Tuple(d.size for d in grid.dimensions)

function PB.get_dimensions(grid::CartesianLinearGrid, space::Type{PB.CellSpace}; expand_cartesian=false)
    if expand_cartesian
        return grid.dimensions
    else
        return (PB.NamedDimension("cells", grid.ncells), )
    end
end

PB.get_dimensions(grid::CartesianLinearGrid{3}, space::Type{PB.ColumnSpace}; expand_cartesian=false) = 
    (PB.NamedDimension("columns", grid.ncolumns), )

_extra_dimensions(grid::CartesianLinearGrid) = grid.dimensions_extra


"""
    CartesianArrayGrid <: PB.AbstractMesh

nD grid with netcdf CF1.0 coordinates, using n-dimensional Arrays for PALEO internal representation of Variables

# Fields
- `ncells::Int64`: number of cells in Domain = product of cartesian dimension sizes
- `dimensions::Vector{PB.NamedDimensions}`: names and sizes of cartesian spatial dimensions (ordered list)
- `coordinates::Dict{String, Vector{String}}`: optional attached coordinates
"""
Base.@kwdef mutable struct CartesianArrayGrid{N} <: PB.AbstractMesh
    
    ncells::Int64                           = -1 
 
    dimensions::Vector{PB.NamedDimension}   = Vector{PB.NamedDimension}(undef, N) # netcdf dimension names
    dimensions_extra::Vector{PB.NamedDimension} = [PB.NamedDimension("bnds", 2)] # additional dimension names

    "index of longitude dimension (if any)"
    londim::Int                            = 0
    "index of latitude dimension (if any)"
    latdim::Int                             = 0
    "index of z dimension (if any)"
    zdim::Int                               = 0
    "index of surface in z coordinate (if any) (1 or length(z dim))" 
    zidxsurface::Int                       = 0
    "multiplier to use for display (eg -1.0 to convert depth to height)"
    display_mult::Vector{Float64}          = ones(N)

    subdomains::Dict{String, PB.AbstractSubdomain}   = Dict{String, PB.AbstractSubdomain}()

    coordinates::Dict{String, Vector{String}} = Dict{String, Vector{String}}()
end

function Base.show(io::IO, grid::CartesianArrayGrid)
    print(io, "CartesianArrayGrid(ncells=", grid.ncells, ", dimensions: ", grid.dimensions, ", dimensions_extra ", grid.dimensions_extra,
        ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", grid::CartesianArrayGrid)
    println(io, "CartesianArrayGrid")
    println(io, "  ncells=", grid.ncells)
    _show_subdomains_dimensions_coordinates(io, grid)
    return nothing
end

available_spaces(grid::CartesianArrayGrid) = (PB.ScalarSpace, PB.CellSpace)
PB.internal_size(space::Type{PB.CellSpace}, grid::CartesianArrayGrid) = Tuple(d.size for d in grid.dimensions)

PB.get_dimensions(grid::CartesianArrayGrid, space::Type{PB.CellSpace}; expand_cartesian=false) =
    grid.dimensions

_extra_dimensions(grid::CartesianArrayGrid) = grid.dimensions_extra

"""
    CartesianGrid(
        gridtype, dimensions::Vector{PB.NamedDimension};
        [londim] [, latdim] [,zdim=0], [,zidxsurface=0], [, ztoheight=1.0]) -> grid::CartesianLinearGrid

Create a CartesianLinearGrid or CartesianArrayGrid from dimensions and coordinates.
"""
function CartesianGrid(
    GridType::Type{<:Union{CartesianLinearGrid, CartesianArrayGrid}},
    dimensions::Vector{PB.NamedDimension};
    londim=findfirst(d-> d.name == "lon", dimensions),
    latdim=findfirst(d-> d.name == "lat", dimensions), 
    zdim=0,
    zidxsurface=0, 
    ztoheight=1.0
) 

    grid = GridType{length(dimensions)}(;
        dimensions=dimensions,
        londim=londim,
        latdim=latdim,
        zdim=zdim,
        zidxsurface=zidxsurface
    )

    if zdim > 0
        grid.display_mult[zdim] = ztoheight
    end

    sz = Tuple(d.size for d in dimensions)
    if GridType == CartesianLinearGrid
        # create linear index (with no points)
        grid.linear_index = Array{Union{Missing, Int32},length(grid.dimensions)}(undef, sz...)
        grid.ncells = 0
    else
        grid.ncells = prod(sz)
    end

    return grid
end
  
"""
    CartesianGrid(GridType, ncfilename::AbstractString, dimnames::Vector{<:AbstractString};
        equalspacededges=false, convert_edges_to_bounds=false) -> grid::GridType, coord_vars, coord_edges_bounds_vars

Read a netcdf file with CF1.0 coordinates, and create corresponding a CartesianLinearGrid or CartesianArrayGrid from dimensions `dimnames`.
`coord_vars` and optional `coord_edges_bounds_vars` are Vectors of `coorddimname => values` for creating coordinate variables.
"""
function CartesianGrid(
    GridType::Type{<:Union{CartesianLinearGrid, CartesianArrayGrid}},
    ncfilename::AbstractString, 
    dimnames::Vector{<:AbstractString};
    equalspacededges=false,
    convert_edges_to_bounds=false,
)   
    @info "CartesianGrid creating $GridType{$(length(dimnames))} from dimnames=$dimnames in netcdf file $(ncfilename)"

    grid = GridType{length(dimnames)}()

    coord_vars = []
    coord_edges_bounds_vars = []

    NCDatasets.Dataset(ncfilename) do ds
        # read dimensions and coordinates
        for i in eachindex(dimnames)
            dimname = dimnames[i]            
            dim = ds.dim[dimname]
            @info "  read dimension $dimname = $dim"

            grid.dimensions[i] = PB.NamedDimension(dimname, dim)
            v_cf = ds[dimname]
            coord_vals = Array(v_cf)
            @assert ndims(coord_vals) == 1
            coord_vals = [cv for cv in coord_vals] # narrow type to eliminate unnecessary Missing
            push!(coord_vars, dimname=>coord_vals)
            if haskey(v_cf.attrib, "edges")
                edgesname = v_cf.attrib["edges"]
                edgesdim = ds.dim[edgesname]            
                @info "  reading coordinate edges from $edgesname = $edgesdim"
                edges_vals = Array(ds[edgesname])              
                @assert size(edges_vals) == (edgesdim, )
                edges_vals = [ev for ev in edges_vals]
                if convert_edges_to_bounds
                    boundsname = dimname*"_bnds"
                    @info "    converting edges to bounds $boundsname"
                    coord_bounds = similar(coord_vals, eltype(coord_vals), 2, length(coord_vals))
                    coord_bounds[1, :] = edges_vals[1:end-1]
                    coord_bounds[2, :] = edges_vals[2:end]
                    push!(coord_edges_bounds_vars, boundsname => coord_bounds)
                else
                    push!(grid.dimensions_extra,  PB.NamedDimension(edgesname, edgesdim))
                    push!(coord_edges_bounds_vars, edgesname => edges_vals)
                end
            elseif equalspacededges
                ew = coord_vals[2] - coord_vals[1]
                if convert_edges_to_bounds
                    boundsname = dimname*"_bnds"
                    @info "  assuming equal spacing $ew to calculate coordinate bounds $boundsname"
                    coord_bounds = similar(coord_vals, eltype(coord_vals), 2, length(coord_vals))
                    coord_bounds[1, :] = coord_vals .- ew/2.0
                    coord_bounds[2, :] = coord_vals .+ ew/2.0
                    push!(coord_edges_bounds_vars, boundsname => coord_bounds)
                else
                    edgesname = dimname*"_edges"
                    @info "  assuming equal spacing $ew to calculate coordinate edges $edgesname"
                    edges_vals = [coord_vals .- ew/2.0; coord_vals[end] + ew/2.0 ]
                    push!(grid.dimensions_extra,  PB.NamedDimension(edgesname, dim + 1))
                    push!(coord_edges_bounds_vars, edgesname => edges_vals)
                end
            else
                @info "   no edges attribute and equalspacededges=false"
            end

            if v_cf.attrib["standard_name"] == "longitude"
                @info "  dim $i got standard_name=='longitude'"
                grid.londim = i
            elseif v_cf.attrib["standard_name"] == "latitude"
                @info "  dim $i got standard_name=='latitude'"
                grid.latdim = i
            elseif v_cf.attrib["standard_name"] == "depth"
                @info "  dim $i got standard_name=='depth', using this as z dimension"
                grid.zdim = i
                grid.display_mult[i] = -1.0
                if coord_vals[1] < coord_vals[end]       
                    grid.zidxsurface = 1 
                else             
                    grid.zidxsurface = length(coord_vals)
                end
                @info "  surface is index $(grid.zidxsurface)"
            end          
        end
    end

    grid_size = Tuple(nd.size for nd in grid.dimensions)
    if GridType == CartesianLinearGrid
        # create linear index (with no points)        
        grid.linear_index = Array{Union{Missing,Int32}, length(grid_size)}(undef, grid_size...)
        grid.ncells = 0
    else
        grid.ncells = prod(grid_size)
    end

    return grid, coord_vars, coord_edges_bounds_vars
end

"""
    set_linear_index(grid::CartesianLinearGrid{3}, v_i, v_j, v_k)

Set 3D grid linear index (given by v_i, v_j, v_k Vectors defining Cartesian [i,j,k] for each v_i[l], v_j[l], v_k[l])
"""
function set_linear_index(grid::CartesianLinearGrid{3}, v_i, v_j, v_k)

    grid.ncells = length(v_i)
    grid.ncolumns = 0

    grid_size = Tuple(nd.size for nd in grid.dimensions)
    grid.linear_index = Array{Union{Missing, Int32}, 3}(undef, grid_size...)
    grid.cartesian_index = Vector{CartesianIndex{3}}(undef, grid.ncells)
    fill!(grid.linear_index, missing)
    for l in eachindex(v_i)
        grid.linear_index[v_i[l], v_j[l], v_k[l]] = l
        grid.cartesian_index[l] = CartesianIndex(v_i[l], v_j[l], v_k[l])
        if v_k[l] == grid.zidxsurface
            grid.ncolumns += 1
        end
    end

    return nothing
end

"""
    set_linear_index(grid::CartesianLinearGrid{2}, v_i, v_j)

Set 2D grid linear index (given by v_i, v_j Vectors defining Cartesian [i,j] for each v_i[l], v_j[l])
"""
function set_linear_index(grid::CartesianLinearGrid{2}, v_i, v_j)
    
    grid.ncells = length(v_i)

    grid_size = Tuple(nd.size for nd in grid.dimensions)
    grid.linear_index = Array{Union{Missing,Int32}, 2}(undef, grid_size)
    grid.cartesian_index = Vector{CartesianIndex{2}}(undef, grid.ncells)
    fill!(grid.linear_index, missing)
    for l in eachindex(v_i)
        grid.linear_index[v_i[l], v_j[l]] = l
        grid.cartesian_index[l] = CartesianIndex(v_i[l], v_j[l])
    end

    return nothing
end


"""
    cartesian_to_internal(grid::CartesianLinearGrid, griddata::AbstractArray) -> lindata::Vector

 Convert Cartesian Array `griddata` to Vector on `grid.linear_index`
 """   
function cartesian_to_internal(grid::CartesianLinearGrid, griddata::AbstractArray)
    size(grid.linear_index) == size(griddata) || error("grid and data size mismatch")
 
    # recreate vector to change type back from eg Union{Missing, Float64} to Float64, where missing will have been added by internal_to_cartesian
    return [griddata[i] for i in grid.cartesian_index]
end

"""
    internal_to_cartesian(grid::CartesianLinearGrid, internaldata::AbstractArray [,missing_value=missing]) -> griddata::Array

 Convert Array `internaldata` (first index on `grid.linear_index`, up to 2 additional indices for additional data dimensions)
 to Cartesian Array `griddata` (with `missing_value` where no data)
 """  
function internal_to_cartesian(grid::CartesianLinearGrid, internaldata::AbstractArray{T, N}; missing_value=missing) where {T, N}
    grid.ncells == size(internaldata, 1) || error("grid and data size mismatch")

    if missing_value isa Missing
        gr_eltype = Union{Missing, T}
    else
        gr_eltype = T
    end

    if N == 1 # internaldata is a Vector; no extra data dimensions
        griddata = similar(grid.linear_index, gr_eltype)
        fill!(griddata, missing_value)
        griddata[grid.cartesian_index] .= internaldata
    elseif N == 2 # 1 extra data dimension``
        griddata = similar(grid.linear_index, gr_eltype, size(grid.linear_index)..., size(internaldata, 2))
        fill!(griddata, missing_value)
        griddata[grid.cartesian_index, :] .= internaldata
    elseif N == 3 # 2 extra data dimensions
        griddata = similar(grid.linear_index, gr_eltype, size(grid.linear_index)..., size(internaldata, 2), size(internaldata, 3))
        fill!(griddata, missing_value)
        griddata[grid.cartesian_index, :, :] .= internaldata
    else
        error("data with $(N-1) extra data dimensions not supported")
    end
    
    return griddata
end


function get_lon_idx(grid::CartesianLinearGrid, linear_idx::Integer)
    cartesian_idx = grid.cartesian_index[linear_idx]
    return cartesian_idx[grid.londim]
end

get_lon_idx(grid::CartesianArrayGrid{N}, cartesian_idx::CartesianIndex{N}) where {N} = cartesian_idx[grid.londim]

function get_lat_idx(grid::CartesianLinearGrid, linear_idx::Integer)
    cartesian_idx = grid.cartesian_index[linear_idx]
    return cartesian_idx[grid.latdim] 
end

get_lat_idx(grid::CartesianArrayGrid{N}, cartesian_idx::CartesianIndex{N}) where {N} = cartesian_idx[grid.latdim]


"""
    linear_indices_cartesian_ranges(grid::CartesianLinearGrid, rangestuple) -> lindices

Find linear indices corresponding to a region in a Cartesian grid

# Example

    lindices = linear_indices_cartesian_ranges(grid, (1:2, 10:20, :))
"""
#=
function linear_indices_cartesian_ranges(grid::CartesianLinearGrid{2}, rangestuple  )
    r1 = rangestuple[1]
    r2 = rangestuple[2]
    return [l for l in grid.linear_index[r1, r2] if !ismissing(l)]
end

function linear_indices_cartesian_ranges(grid::CartesianLinearGrid{3}, (r1, r2, r3)  )
    return [l for l in grid.linear_index[r1, r2, r3] if !ismissing(l)]
end
=#

"""
    cellrange_cartesiantile

Create a range of cells within a region of CartesianLinearGrid specified by `rangestuple` in a specified [`PB.Domain`](@ref)

`rangestuple` is a tuple of Cartesian index ranges eg (1:9, :, :) for a 3D grid.
"""
function cellrange_cartesiantile(domain, grid::CartesianLinearGrid{2}, rangestuple; operatorID=0)
    indices = Int[]
    r1, r2 = _expandrangetuple(grid.dimensions, rangestuple[1:2])
   
    for i=1:length(grid.cartesian_index)
        ci = grid.cartesian_index[i]
        if ci[1] in r1 && ci[2] in r2
            push!(indices, i)
        end
    end

    return PB.CellRange(domain=domain, indices=PB.replace_contiguous_range(indices), operatorID=operatorID)
end

"3D case return a CellRangeColumns"
function cellrange_cartesiantile(domain, grid::CartesianLinearGrid{3}, rangestuple; operatorID=0)

    grid.zdim == 3 || error("grid.zdim = $(grid.zdim) not supported (must be last dimension = 3)")
    # recreate (all) surface indices
    surfindices = [ci for ci in grid.cartesian_index if ci[3] == grid.zidxsurface]

    r1, r2, r3 = _expandrangetuple(grid.dimensions, rangestuple)
    indices = Int[]
    colindices = Vector{Pair{Int, Vector{Int}}}()
    # iterate through surface indices and generate columns 
    for is in 1:length(surfindices)
        sci = surfindices[is]
        # println("is ", is, " sci ", sci)
        if sci[1] in r1 && sci[2] in r2
            colinds = Int[]
            # iterate through column in order surface -> floor
            if grid.zidxsurface == 1
                krange = 1:size(grid.linear_index,3)
            else
                krange = reverse(1:size(grid.linear_index,3))
            end
            for k in krange
                li = grid.linear_index[sci[1], sci[2], k ]
                if k in r3 && !ismissing(li)
                    push!(colinds, li)
                end
            end
            if !isempty(colinds)
                colinds = PB.replace_contiguous_range(colinds)
                append!(indices, colinds)
                push!(colindices, is=>colinds)
            end
        end
    end

    return PB.CellRangeColumns(
        domain=domain, 
        indices=PB.replace_contiguous_range(indices),
        columns=colindices,
        operatorID=operatorID
    )
end

function cellrange_cartesiantile(domain, grid::CartesianArrayGrid{2}, rangestuple; operatorID=0)
    rt = _expandrangetuple(grid.dimensions, rangestuple[1:2])

    return PB.CellRange(domain=domain, indices=CartesianIndices(rt), operatorID=operatorID)
end

function cellrange_cartesiantile(domain, grid::CartesianArrayGrid{3}, rangestuple; operatorID=0)
    rt = _expandrangetuple(grid.dimensions, rangestuple[1:3])

    return PB.CellRange(domain=domain, indices=CartesianIndices(rt), operatorID=operatorID)
end

"""
    create_default_cellrange(domain, grid::Union{CartesianLinearGrid, CartesianArrayGrid} [; operatorID=0]) -> CellRangeColumns

Create a CellRange for entire `domain`. Return a [`PB.CellRangeColumns`](@ref) provided by [`cellrange_cartesiantile`](@ref)
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::Union{CartesianLinearGrid, CartesianArrayGrid}; operatorID=0)
    return cellrange_cartesiantile(domain, grid, (:, :, :), operatorID=operatorID)
end

"expand (1:2, :, :) replacing : with ranges"
function _expandrangetuple(nameddims::Vector{PB.NamedDimension},  rangestuple)
    ranges = []
    for i = 1:length(rangestuple)
        r = rangestuple[i]
        if r isa Colon
            push!(ranges, 1:nameddims[i].size)
        else
            push!(ranges, r)
        end
    end
    return Tuple(ranges)
end

"return partitioning of a 3D Domain into n tiles"
function cellrange_cartesiantile(domain::PB.AbstractDomain, grid::CartesianLinearGrid{3}, ntiles::Int; operatorID=0)
    # get entire Domain
    default_cellrange = create_default_cellrange(domain, grid)
    cumulative_cells = Int[] # total cells including column i
    cc = 0
    for (icol, cellscol) in default_cellrange.columns
        cc += length(cellscol)
        push!(cumulative_cells, cc)
    end
    total_cols = length(default_cellrange.columns)
    total_cells = last(cumulative_cells)
    @info "cellrange_cartesiantile total_cols=$total_cols, total_cells=$total_cells"

    cellranges = []
    firsttilecol = 1
    for t in 1:ntiles
        indices = Int[]
        colindices = Vector{Pair{Int,Vector{Int}}}()
        tilecells = 0
        for tc in firsttilecol:total_cols
            (icol, cellscol) = default_cellrange.columns[tc]
            append!(indices, cellscol)
            push!(colindices, icol=>cellscol)
            tilecells += length(cellscol)
            
            if tilecells >= total_cells/ntiles && t != ntiles
                firsttilecol = tc+1
                break
            end
        end
        @info "  tile $t columns $(first(colindices)[1]):$(last(colindices)[1]) total cells $(length(indices))"
        push!(
            cellranges, 
            PB.CellRangeColumns(
                domain=domain, 
                indices=PB.replace_contiguous_range(indices), 
                columns=colindices, 
                operatorID=operatorID
            )
        )
    end

    return cellranges
end

"derive CellRange for a 2D boundary Domain from a CellRangeColumns for the 3D interior Domain"
function cellrange_from_columns(boundarydomain, grid::CartesianLinearGrid{2}, interior_crcolumns::PB.CellRangeColumns)
    indices = [icol for (icol, cellscol) in interior_crcolumns.columns]
    return PB.CellRange(
        domain=boundarydomain, 
        indices=PB.replace_contiguous_range(indices),
        operatorID=interior_crcolumns.operatorID
    )
end

"""
    linear_kji_order(grid, v_i, v_j, v_k) -> perm::Vector

Find `perm` such that indices `v_i[perm], v_j[perm], v_k[perm]` correspond to 
traversal of the 3D grid in order k, j, i
"""
function linear_kji_order(grid, v_i, v_j, v_k)
    
    length(grid.dimensions) == 3 || error("linear_kji_order grid is not a 3D grid")

    # construct a 3D Array and fill with linear index values
    grid_size = Tuple(nd.size for nd in grid.dimensions)
    linear_index = Array{Union{Missing,Int32},length(grid_size)}(missing, grid_size...)
    for l in eachindex(v_i)
        linear_index[v_i[l], v_j[l], v_k[l]] = l
    end

    # visit 3D Array in k, j, i order and record linear index values
    perm = Vector{Int}(undef, length(v_i))
    ip = 0
    for i in 1:grid.dimensions[1].size
        for j in 1:grid.dimensions[2].size
            for k in 1:grid.dimensions[3].size
                if !ismissing(linear_index[i, j, k])
                    ip += 1
                    perm[ip] = linear_index[i, j, k]
                end
            end
        end
    end

    return (perm, sortperm(perm))
end



###########################
# Grid tiling 
###########################

"""
    get_tiled_cellranges(model::Model, tiles [; operatorID])  -> cellranges::Vector{Vector}

Partition Domains into tiles, assuming a single Cartesian gridded Domain + associated boundaries (eg ocean + oceansurface + ...)
and that all additional Domains are scalar Domains (eg atmosphere).

`tiles` should be a collection of tuples of ranges eg `tiles = [(1:9, :, :), (10:16, :, :), (17:25, :, :), (26:36, :, :)]`.

Returns `cellranges`, a Vector of Vectors of CellRanges providing one Vector of Cellranges per tile.
"""
function get_tiled_cellranges(model::PB.Model, tiles; operatorID=0)
    cellranges = []  # vector of vectors, 1 per tile

    for tr in tiles
        @info "creating CellRanges for tile $tr operatorID $operatorID"
        tcellranges = []
        for dom in model.domains            
            if !isnothing(dom.grid)
                push!(tcellranges, cellrange_cartesiantile(dom, dom.grid, tr, operatorID=operatorID))
                @info "  add CellRange for domain $(dom.name) length $(length(last(tcellranges).indices))"
            end
        end
        push!(cellranges, tcellranges)
    end
    # create cellranges for scalar Domains
    cellrangesscalar = []
    @info "creating CellRanges for scalar domains"
    for dom in model.domains            
        if isnothing(dom.grid)
            push!(cellrangesscalar, create_default_cellrange(dom, dom.grid, operatorID=operatorID))
            @info "  added CellRange for scalar domain $(dom.name)"
        end
    end
    @info "adding CellRanges for scalar domains to first tile"
    append!(first(cellranges), cellrangesscalar)

    return cellranges
end

"""
    get_tiled_cellranges(model::Model,  ntiles::Int, interior_domain_name [; operatorID])  -> cellranges::Vector{Vector}

Partition Domains into `ntiles` with approximately equal numbers of cells for `interior_domain_name`, 
assuming a single Cartesian gridded Domain `interior_domain_name` + associated boundaries (eg ocean + oceansurface + ...)
and that all additional Domains are scalar Domains (eg atmosphere).

Returns `cellranges`, a Vector of Vectors of CellRanges providing one Vector of Cellranges per tile.
"""
function get_tiled_cellranges(model::PB.Model, ntiles::Int, interior_domain_name::AbstractString; operatorID=0)
 
    interior_domain = PB.get_domain(model, interior_domain_name)
    interior_cellranges = cellrange_cartesiantile(interior_domain, interior_domain.grid, ntiles, operatorID=operatorID)
  
    io = IOBuffer()
    println(io, "get_tiled_cellranges:")

    cellranges = [] # vector of vectors, 1 per tile
    for it in 1:ntiles
        println(io, "  creating CellRanges for tile $it")      
        tcellranges = []
        push!(tcellranges, interior_cellranges[it])
        for dom in model.domains            
            if !isnothing(dom.grid) && dom != interior_domain
                push!(tcellranges, cellrange_from_columns(dom, dom.grid, interior_cellranges[it]))            
                println(io, "    add CellRange for domain $(dom.name) length $(length(last(tcellranges).indices))")
            end
        end
        push!(cellranges, tcellranges)
    end
    # create cellranges for scalar Domains
    cellrangesscalar = []
    println(io, "  creating CellRanges for scalar domains")
    for dom in model.domains            
        if isnothing(dom.grid)
            push!(cellrangesscalar, create_default_cellrange(dom, dom.grid, operatorID=operatorID))
            println(io, "    added CellRange for scalar domain $(dom.name)")
        end
    end
    println(io, "  adding CellRanges for scalar domains to first tile")
    
    append!(first(cellranges), cellrangesscalar)

    @info String(take!(io))

    return cellranges
end


end # module
