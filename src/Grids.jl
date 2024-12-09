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

[`PB.internal_size`](@ref), optionally [`PB.cartesian_size`](@ref)

[`PB.Grids.set_subdomain!`](@ref), [`PB.Grids.get_subdomain`](@ref)

[`PB.Grids.create_default_cellrange`](@ref), [`PB.Grids.get_region`](@ref)
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
function create_default_cellrange(domain::PB.AbstractDomain, grid::Union{PB.AbstractMesh, Nothing}) end

"""
    get_region(grid::Union{PB.AbstractMesh, Nothing}, values; selectargs...) -> values_subset, (dim_subset::NamedDimension, ...)

Return the subset of `values` given by `selectargs` (Grid-specific keywords eg cell=, column=, ...)
and corresponding dimensions (with attached coordinates).
"""
function get_region(grid::Union{PB.AbstractMesh, Nothing}, values) end

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
 

# generic handler when subdomainname present
function PB.internal_size(space::Type{<:PB.AbstractSpace}, grid::Union{PB.AbstractMesh, Nothing}, subdomainname::AbstractString)
    if isempty(subdomainname)
        if space === PB.ScalarSpace
            return (1, )
        else
            return PB.internal_size(space, grid)
        end
    else
        if space === PB.CellSpace
            subdomain = get_subdomain(grid, subdomainname)
            return (length(subdomain.indices),)
        else
            error("internal_size: cannot specify subdomain with space=$space (subdomainname=$subdomainname)")
        end
    end
end

# scalar space always supported
PB.internal_size(space::Type{PB.ScalarSpace}, grid::Union{PB.AbstractMesh, Nothing}) = (1, )
# ignore subdomain 
PB.internal_size(space::Type{PB.ScalarSpace}, grid::Union{PB.AbstractMesh, Nothing}, subdomainname::AbstractString) = (1, )

###################################
# Fallbacks for Domain with no grid
#####################################

# allow Vector variables length 1 in a 0D Domain without a grid
PB.internal_size(space::Type{PB.CellSpace}, grid::Nothing) = (1, )

# allow single dimension
PB.cartesian_size(grid::Nothing) = (1, )
cartesian_to_internal(grid::Nothing, griddata::AbstractArray) = griddata

get_subdomain(grid::Nothing, subdomainname::AbstractString) = error("get_subdomain: no subdomain $subdomainname")

"""
    create_default_cellrange(domain, grid::Nothing [; operatorID=0]) -> CellRange

Create a CellRange for entire `domain`. Fallback for a domain with no grid
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::Nothing; operatorID=0)
    return PB.CellRange(domain=domain, indices=1:PB.get_length(domain), operatorID=operatorID)
end

"""
    get_region(grid::Nothing, values) -> values[]

Fallback for Domain with no grid, assumed 1 cell
"""
function get_region(grid::Nothing, values)
    length(values) == 1 ||
        throw(ArgumentError("grid==Nothing and length(values) != 1"))
    return values[], ()
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

    subdomains::Dict{String, PB.AbstractSubdomain}   = Dict{String, PB.AbstractSubdomain}()
end

function Base.show(io::IO, grid::UnstructuredVectorGrid)
    print(io, "UnstructuredVectorGrid(ncells=", grid.ncells, 
        ", cellnames=", grid.cellnames, ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

PB.internal_size(space::Type{PB.CellSpace}, grid::UnstructuredVectorGrid) = (grid.ncells, )

# single dimension
PB.cartesian_size(grid::UnstructuredVectorGrid) = (grid.ncells, )
cartesian_to_internal(grid::UnstructuredVectorGrid, griddata::AbstractArray) = griddata

"""
    get_region(grid::UnstructuredVectorGrid, values; cell=Nothing) -> 
        values_subset, (dim_subset::NamedDimension, ...)

# Keywords for region selection:
- `cell::Union{Nothing, Int, Symbol}`: an Int for cell number (first cell = 1), or a Symbol to look up in `cellnames`
  `cell = Nothing` is also allowed for a single-cell Domain.
"""
function get_region(grid::UnstructuredVectorGrid, values; cell::Union{Nothing, Int, Symbol}=nothing)
    if isnothing(cell)
        grid.ncells == 1 || throw(ArgumentError("'cell' argument (an Int or Symbol) is required for an UnstructuredVectorGrid with > 1 cell"))
        idx = 1
    elseif cell isa Int
        idx = cell
    else
        idx = get(grid.cellnames, cell, nothing)
        !isnothing(idx) ||
            throw(ArgumentError("cell ':$cell' not present in  grid.cellnames=$(grid.cellnames)"))
    end

    return (
        values[idx],
        (),  # no dimensions (ie squeeze out a dimension length 1 for single cell)
    )
end


"""
    create_default_cellrange(domain, grid::UnstructuredVectorGrid [; operatorID=0]) -> CellRange

Create a CellRange for entire `domain` - use linear index.
"""
function create_default_cellrange(domain::PB.AbstractDomain, grid::UnstructuredVectorGrid; operatorID=0)
    return PB.CellRange(domain=domain, indices=1:grid.ncells, operatorID=operatorID)
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
- `z_coords::Vector{FixedCoord}`: z coordinates of cell mid point, lower surface, upper surface
- `columnnames::Vector{Symbol}:` optional column names
"""
Base.@kwdef mutable struct UnstructuredColumnGrid <: PB.AbstractMesh
    ncells::Int64    
    Icolumns::Vector{Vector{Int}}    

    z_coords::Vector{PB.FixedCoord} = PB.FixedCoord[]

    "Define optional column names"
    columnnames::Vector{Symbol} = Symbol[]

    subdomains::Dict{String, PB.AbstractSubdomain} = Dict{String, PB.AbstractSubdomain}()
end

function Base.show(io::IO, grid::UnstructuredColumnGrid)
    print(io, "UnstructuredColumnGrid(ncells=", grid.ncells, ", columns: ", length(grid.Icolumns),
        ", columnnames=", grid.columnnames, ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

PB.internal_size(space::Type{PB.CellSpace}, grid::UnstructuredColumnGrid) = (grid.ncells, )
PB.internal_size(space::Type{PB.ColumnSpace}, grid::UnstructuredColumnGrid) = (length(grid.Icolumns), )

# single dimension
PB.cartesian_size(grid::UnstructuredColumnGrid) = (grid.ncells, )
cartesian_to_internal(grid::UnstructuredColumnGrid, griddata::AbstractArray) = griddata

"""
    get_region(grid::UnstructuredColumnGrid, values; column, [cell=nothing]) -> 
        values_subset, (dim_subset::NamedDimension, ...)

# Keywords for region selection:
- `column::Union{Int, Symbol}`: (may be an Int, or a Symbol to look up in `columnames`)
- `cell::Int`: optional cell index within `column`, highest cell is cell 1
"""
function get_region(grid::UnstructuredColumnGrid, values; column, cell::Union{Nothing, Int}=nothing)

    if column isa Int
        column in 1:length(grid.Icolumns) ||
            throw(ArgumentError("column index $column out of range"))
        colidx = column
    else
        colidx = findfirst(isequal(column), grid.columnnames)
        !isnothing(colidx) || 
            throw(ArgumentError("columnname '$column' not present in  grid.columnnames=$(grid.columnnames)"))
    end

    if isnothing(cell)
        indices = grid.Icolumns[colidx]
        return (
            values[indices],
            (PB.NamedDimension("z", length(indices), PB.get_region(grid.z_coords, indices)), ),
        )
    else
        # squeeze out z dimension
        idx = grid.Icolumns[colidx][cell]
        return (
            values[idx],
            (),  # no dimensions (ie squeeze out a dimension length 1 for single cell)
        )
    end
    
end

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
- `dimnames::Vector{String}`: names of dimensions (ordered list)
- `dims::Vector{Int}`: sizes of dimensions (ordered list)
- `coords::Vector{Vector{Float64}}`: attached cell-centre coordinates for each dimension (ordered list)
- `coords_edges::Vector{Vector{Float64}}`: attached cell-edge coordinates for each dimension (ordered list)
"""
Base.@kwdef mutable struct CartesianLinearGrid{N} <: PB.AbstractMesh
    
    ncells::Int64                           = -1 
    ncolumns::Int64                         = -1

    dimnames::Vector{String}                = Vector{String}(undef, N)  # netcdf dimension names
    dims::Vector{Int}                       = Vector{Int}(undef, N)
    coords::Vector{Vector{Float64}}         = Vector{Vector{Float64}}()
    coords_edges::Vector{Vector{Float64}}   = Vector{Vector{Float64}}()

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

    "cartesian -> linear index mapping (initialized to `missing` ie linear index contains 0 cells)"
    linear_index::Array{Union{Missing,Int32}, N} = Array{Union{Missing,Int32},N}(undef, zeros(Int,N)...) 
    "linear -> cartesian index mapping (initialized to empty Vector ie linear index contains 0 cells)"
    cartesian_index::Vector{CartesianIndex{N}}  = Vector{CartesianIndex{N}}()
end

function Base.show(io::IO, grid::CartesianLinearGrid)
    dimtuple = NamedTuple{Tuple(Symbol.(grid.dimnames))}(Tuple(grid.dims))
    print(io, "CartesianLinearGrid(ncells=", grid.ncells, ", dimensions: ", dimtuple,
        ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

PB.internal_size(space::Type{PB.CellSpace}, grid::CartesianLinearGrid) = (grid.ncells, )
PB.internal_size(space::Type{PB.ColumnSpace}, grid::CartesianLinearGrid{3}) = (grid.ncolumns, )

PB.cartesian_size(grid::CartesianLinearGrid) = Tuple(grid.dims)


"""
    CartesianArrayGrid <: PB.AbstractMesh

nD grid with netcdf CF1.0 coordinates, using n-dimensional Arrays for PALEO internal representation of Variables

# Fields
- `ncells::Int64`: number of cells in Domain = `length(linear_index)` (may be subset of total points in `prod(dims)`)
- `dimnames::Vector{String}`: names of dimensions (ordered list)
- `dims::Vector{Int}`: sizes of dimensions (ordered list)
- `coords::Vector{Vector{Float64}}`: attached cell-centre coordinates for each dimension (ordered list)
- `coords_edges::Vector{Vector{Float64}}`: attached cell-edge coordinates for each dimension (ordered list)
"""
Base.@kwdef mutable struct CartesianArrayGrid{N} <: PB.AbstractMesh
    
    ncells::Int64                           = -1 
 
    dimnames::Vector{String}                = Vector{String}(undef, N)  # netcdf dimension names
    dims::Vector{Int}                       = Vector{Int}(undef, N)
    coords::Vector{Vector{Float64}}         = Vector{Vector{Float64}}()
    coords_edges::Vector{Vector{Float64}}   = Vector{Vector{Float64}}()

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
end

function Base.show(io::IO, grid::CartesianArrayGrid)
    dimtuple = NamedTuple{Tuple(Symbol.(grid.dimnames))}(Tuple(grid.dims))
    print(io, "CartesianArrayGrid(ncells=", grid.ncells, ", dimensions: ", dimtuple,
        ", subdomains: ", keys(grid.subdomains), ")")
    return nothing
end

PB.internal_size(space::Type{PB.CellSpace}, grid::CartesianArrayGrid) = Tuple(grid.dims)

PB.cartesian_size(grid::CartesianArrayGrid) = Tuple(grid.dims)


"""
    get_region(grid::Union{CartesianLinearGrid{2}, CartesianArrayGrid{2}} , internalvalues; [i=i_idx], [j=j_idx]) ->
        arrayvalues_subset, (dim_subset::NamedDimension, ...)

# Keywords for region selection:
- `i::Int`: optional, slice along first dimension
- `j::Int`: optional, slice along second dimension

`internalvalues` are transformed if needed from internal Field representation as a Vector length `ncells`, to
an Array (2D if neither i, j arguments present, 1D if i or j present, 0D ie one cell if both present)
"""
function get_region(
    grid::Union{CartesianLinearGrid{2}, CartesianArrayGrid{2}}, internalvalues; 
    i::Union{Integer, Colon}=Colon(), j::Union{Integer, Colon}=Colon()
)
    return _get_region(grid, internalvalues, [i, j])
end

"""
    get_region(grid::Union{CartesianLinearGrid{3}, CartesianArrayGrid{3}}, internalvalues; [i=i_idx], [j=j_idx]) ->
        arrayvalues_subset, (dim_subset::NamedDimension, ...)

# Keywords for region selection:
- `i::Int`: optional, slice along first dimension
- `j::Int`: optional, slice along second dimension
- `k::Int`: optional, slice along third dimension

`internalvalues` are transformed if needed from internal Field representation as a Vector length `ncells`, to
an Array (3D if neither i, j, k arguments present, 2D if one of i, j or k present, 1D if two present,
0D ie one cell if i, j, k all specified).
"""
function get_region(
    grid::Union{CartesianLinearGrid{3}, CartesianArrayGrid{3}}, internalvalues;
    i::Union{Integer, Colon}=Colon(), j::Union{Integer, Colon}=Colon(), k::Union{Integer, Colon}=Colon()
)
    return _get_region(grid, internalvalues, [i, j, k])
end

function _get_region(
    grid::Union{CartesianLinearGrid, CartesianArrayGrid}, internalvalues, indices
)
    if !isempty(grid.coords) && !isempty(grid.coords_edges)
        dims = [
            PB.NamedDimension(grid.dimnames[idx], grid.coords[idx], grid.coords_edges[idx])
            for (idx, ind) in enumerate(indices) if isa(ind, Colon)
        ]
    elseif !isempty(grid.coords)
        dims = [
            PB.NamedDimension(grid.dimnames[idx], grid.coords[idx])
            for (idx, ind) in enumerate(indices) if isa(ind, Colon)
        ]
    else
        dims = [
            PB.NamedDimension(grid.dimnames[idx])
            for (idx, ind) in enumerate(indices) if isa(ind, Colon)
        ]
    end

    values = internal_to_cartesian(grid, internalvalues)
    if !all(isequal(Colon()), indices)
        values = values[indices...]
    end

    return values, Tuple(dims)    
end

"""
    CartesianGrid(
        gridtype, dimnames::Vector{<:AbstractString}, dims, coords, coords_edges;
        [londim] [, latdim] [,zdim=0], [,zidxsurface=0], [, ztoheight=1.0]) -> grid::CartesianLinearGrid

Create a CartesianLinearGrid or CartesianArrayGrid from dimensions and coordinates.
"""
function CartesianGrid(
    GridType::Type{<:Union{CartesianLinearGrid, CartesianArrayGrid}},
    dimnames::Vector{<:AbstractString}, dims, coords=Vector{Vector{Float64}}(), coords_edges=Vector{Vector{Float64}}();
    londim=findfirst(isequal("lon"), dimnames),
    latdim=findfirst(isequal("lat"), dimnames), 
    zdim=0,
    zidxsurface=0, 
    ztoheight=1.0
) 

    grid = GridType{length(dimnames)}(
        dimnames=dimnames, 
        dims=dims, 
        coords=coords, 
        coords_edges=coords_edges,
        londim=londim,
        latdim=latdim,
        zdim=zdim,
        zidxsurface=zidxsurface
    )

    if zdim > 0
        grid.display_mult[zdim] = ztoheight
    end

    if GridType == CartesianLinearGrid
        # create linear index (with no points)
        grid.linear_index = Array{Union{Missing,Int32},length(grid.dims)}(undef, grid.dims...)
        grid.ncells = 0
    else
        grid.ncells = prod(grid.dims)
    end

    return grid
end
  
"""
    CartesianGrid(griddtype, ncfilename::AbstractString, dimnames::Vector{<:AbstractString}; equalspacededges=false) -> grid::CartesianLinearGrid

Read a netcdf file with CF1.0 coordinates, and create corresponding a CartesianLinearGrid or CartesianArrayGrid from dimensions `dimnames`.
"""
function CartesianGrid(
    GridType::Type{<:Union{CartesianLinearGrid, CartesianArrayGrid}},
    ncfilename::AbstractString, 
    dimnames::Vector{<:AbstractString};
    equalspacededges=false
)   
    @info "CartesianGrid creating $GridType{$(length(dimnames))} from dimnames=$dimnames in netcdf file $(ncfilename)"

    grid = GridType{length(dimnames)}()
    grid.coords = Vector{Vector{Float64}}(undef, length(dimnames))
    grid.coords_edges = Vector{Vector{Float64}}(undef, length(dimnames))

    NCDatasets.Dataset(ncfilename) do ds
        # read dimensions and coordinates
        for i in eachindex(dimnames)
            dimname = dimnames[i]
            grid.dimnames[i] = dimname
            grid.dims[i] = ds.dim[dimname]
            @info "  read dimension $(dimname) = $(grid.dims[i])"

            v_cf = ds[dimname]
            grid.coords[i] = Array(v_cf)
            if haskey(v_cf.attrib, "edges")
                edgesname = v_cf.attrib["edges"]
                @info "  reading coordinate edges from $edgesname"
                grid.coords_edges[i] = Array(ds[edgesname])
            elseif equalspacededges
                es = grid.coords[i][2] - grid.coords[i][1]
                @info "  assuming equal spacing $es to calculate coordinate edges"
                grid.coords_edges[i] = [grid.coords[i] .- es/2.0; grid.coords[i][end] + es/2.0 ]
            else
                error("   no edges attribute and equalspacededges=false")
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
                if grid.coords[i][1] < grid.coords[i][end]       
                    grid.zidxsurface = 1 
                else             
                    grid.zidxsurface = length(grid.coords[i])
                end
                @info "  surface is index $(grid.zidxsurface)"
            end          
        end
    end

    if GridType == CartesianLinearGrid
        # create linear index (with no points)
        grid.linear_index = Array{Union{Missing,Int32},length(grid.dims)}(undef, grid.dims...)
        grid.ncells = 0
    else
        grid.ncells = prod(grid.dims)
    end

    return grid
end

"""
    set_linear_index(grid::CartesianLinearGrid{3}, v_i, v_j, v_k)

Set 3D grid linear index (given by v_i, v_j, v_k Vectors defining Cartesian [i,j,k] for each v_i[l], v_j[l], v_k[l])
"""
function set_linear_index(grid::CartesianLinearGrid{3}, v_i, v_j, v_k)

    grid.ncells = length(v_i)
    grid.ncolumns = 0

    grid.linear_index = Array{Union{Missing,Int32}, 3}(undef, grid.dims...)
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

    grid.linear_index = Array{Union{Missing,Int32}, 2}(undef, grid.dims...)
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

    return griddata[grid.cartesian_index]
end

cartesian_to_internal(grid::CartesianArrayGrid, griddata::AbstractArray) = griddata

"""
    internal_to_cartesian(grid::CartesianLinearGrid, internaldata::AbstractVector [,missing_value=missing]) -> griddata::Array

 Convert Vector `internaldata` (on `grid.linear_index`) to Cartesian Array `griddata` (with `missing_value` where no data)
 """  
function internal_to_cartesian(grid::CartesianLinearGrid, internaldata::AbstractArray; missing_value=missing)
    grid.ncells == length(internaldata) || error("grid and data size mismatch")

    if missing_value isa Missing
        gr_eltype = Union{Missing, eltype(internaldata)}
    else
        gr_eltype = eltype(internaldata)
    end
    griddata = similar(grid.linear_index, gr_eltype)
    fill!(griddata, missing_value)
    griddata[grid.cartesian_index] .= internaldata
    
    return griddata
end

internal_to_cartesian(grid::CartesianArrayGrid, internaldata::AbstractArray; missing_value=missing) = internaldata

function get_lon(grid::CartesianLinearGrid, linear_idx::Integer)
    cartesian_idx = grid.cartesian_index[linear_idx]
    return _get_coord(grid.coords, grid.londim, "longitude", cartesian_idx)  
end

get_lon(grid::CartesianArrayGrid{N}, cartesian_idx::CartesianIndex{N}) where {N} =
    _get_coord(grid.coords, grid.londim, "longitude", cartesian_idx)

function get_lon_edges(grid::CartesianLinearGrid, linear_idx::Integer)
    cartesian_idx = grid.cartesian_index[linear_idx]
    return _get_coord_edges(grid.coords_edges, grid.londim, "longitude", cartesian_idx)
end

get_lon_edges(grid::CartesianArrayGrid{N}, cartesian_idx::CartesianIndex{N}) where {N} =
    _get_coord_edges(grid.coords_edges, grid.londim, "longitude", cartesian_idx)
    

function get_lat(grid::CartesianLinearGrid, linear_idx::Integer)
    cartesian_idx = grid.cartesian_index[linear_idx]
    return _get_coord(grid.coords, grid.latdim, "latitude", cartesian_idx)   
end

get_lat(grid::CartesianArrayGrid{N}, cartesian_idx::CartesianIndex{N}) where {N} =
    _get_coord(grid.coords, grid.latdim, "latitude", cartesian_idx)

function get_lat_edges(grid::CartesianLinearGrid, linear_idx::Integer)
    cartesian_idx = grid.cartesian_index[linear_idx]
    return _get_coord_edges(grid.coords_edges, grid.latdim, "latitude", cartesian_idx)   
end

get_lat_edges(grid::CartesianArrayGrid{N}, cartesian_idx::CartesianIndex{N}) where {N} =
    _get_coord_edges(grid.coords_edges, grid.latdim, "latitude", cartesian_idx)
    

function _get_coord(coords, cdim, dimname, cartesian_idx::CartesianIndex)
    cdim > 0 || error("grid has no $dimname dimension")
    return coords[cdim][cartesian_idx[cdim]]    
end

function _get_coord_edges(coords_edges, cdim, dimname, cartesian_idx::CartesianIndex)
    cdim > 0 || error("grid has no $dimname dimension")
    return (
        coords_edges[cdim][cartesian_idx[cdim]],
        coords_edges[cdim][cartesian_idx[cdim]+1],
    )
end

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
    r1, r2 = _expandrangetuple(grid.dims, rangestuple[1:2])
   
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

    r1, r2, r3 = _expandrangetuple(grid.dims, rangestuple)
    indices = Int[]
    colindices = Vector{Pair{Int,Vector{Int}}}()
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
    rt = _expandrangetuple(grid.dims, rangestuple[1:2])

    return PB.CellRange(domain=domain, indices=CartesianIndices(rt), operatorID=operatorID)
end

function cellrange_cartesiantile(domain, grid::CartesianArrayGrid{3}, rangestuple; operatorID=0)
    rt = _expandrangetuple(grid.dims, rangestuple[1:3])

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
function _expandrangetuple(dims,  rangestuple)
    ranges = []
    for i = 1:length(rangestuple)
        r = rangestuple[i]
        if r isa Colon
            push!(ranges, 1:dims[i])
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
    
    length(grid.dims) == 3 || error("linear_kji_order grid is not a 3D grid")

    # construct a 3D Array and fill with linear index values
    linear_index = Array{Union{Missing,Int32},length(grid.dims)}(missing, grid.dims...)
    for l in eachindex(v_i)
        linear_index[v_i[l], v_j[l], v_k[l]] = l
    end

    # visit 3D Array in k, j, i order and record linear index values
    perm = Vector{Int}(undef, length(v_i))
    ip = 0
    for i in 1:grid.dims[1]
        for j in 1:grid.dims[2]
            for k in 1:grid.dims[3]
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
  
    cellranges = [] # vector of vectors, 1 per tile
    for it in 1:ntiles
        println("creating CellRanges for tile $it")      
        tcellranges = []
        push!(tcellranges, interior_cellranges[it])
        for dom in model.domains            
            if !isnothing(dom.grid) && dom != interior_domain
                push!(tcellranges, cellrange_from_columns(dom, dom.grid, interior_cellranges[it]))            
                println("  add CellRange for domain $(dom.name) length $(length(last(tcellranges).indices))")
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


end # module
