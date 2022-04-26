module GridReactions

import PALEOboxes as PB

import Infiltrator # Julia debugger


"""
    ReactionUnstructuredVectorGrid

Create an [`PB.Grids.UnstructuredVectorGrid`](@ref) with `ncells` (from `ncells` Parameter).
"""
Base.@kwdef mutable struct ReactionUnstructuredVectorGrid{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParInt("ncells", 1,
            description="number of grid cells"),
    )
   
end
   
function PB.set_model_geometry(rj::ReactionUnstructuredVectorGrid, model::PB.Model)
    @info "set_model_geometry $(PB.fullname(rj))"

    grid = PB.Grids.UnstructuredVectorGrid(ncells=rj.pars.ncells.v)
    rj.domain.grid = grid
    
    @info "  set $(rj.domain.name) Domain size=$(grid.ncells) grid=$grid"
         
    return nothing
end

PB.register_methods!(rj::ReactionUnstructuredVectorGrid) = nothing


"""
    ReactionCartesianGrid

Create a [`PB.Grids.CartesianArrayGrid`](@ref) with `dims` and `dimnames``
"""
Base.@kwdef mutable struct ReactionCartesianGrid{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParStringVec("dimnames", ["lat", "lon", "z"],
            description="grid dimension names"),
        PB.ParIntVec("dims", [2, 3, 4],
            description="grid dimensions"),        
    )
   
end
   
function PB.set_model_geometry(rj::ReactionCartesianGrid, model::PB.Model)
    @info "set_model_geometry $(PB.fullname(rj))"

    grid = PB.Grids.CartesianGrid(
        PB.Grids.CartesianArrayGrid,
        rj.pars.dimnames.v, rj.pars.dims.v
    )
    
    rj.domain.grid = grid
    
    @info "  set $(rj.domain.name) Domain size=$(grid.ncells) grid=$grid"
         
    return nothing
end

PB.register_methods!(rj::ReactionCartesianGrid) = nothing


"""
    ReactionGrid2DNetCDF

Create a 2D [`PB.Grids.CartesianLinearGrid`](@ref) from grid information in a NetCDF file.
"""
Base.@kwdef mutable struct ReactionGrid2DNetCDF{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParType(PB.AbstractMesh, "grid_type", PB.Grids.CartesianArrayGrid,
            allowed_values=(PB.Grids.CartesianArrayGrid, PB.Grids.CartesianLinearGrid),
            description="Cartesian grid type to create"),

        PB.ParString("grid_file", "",
            description="netcdf file with 2D grid information"),

        PB.ParStringVec("coordinate_names", ["longitude", "latitude"],
            description="coordinate names to read from netcdf file"),

        PB.ParBool("equalspacededges", false,
            description="true to calculate cell edges assuming an equal-spaced grid"),

        PB.ParString("area_var", "",
            description="netcdf variable with cell area (if available)"),

        PB.ParDouble("planet_radius", 6.371229e6,
            description="radius to calculate cell area from spherical geometry (if area_var = \"\")"),
    )
   
end
   
function PB.set_model_geometry(rj::ReactionGrid2DNetCDF, model::PB.Model)
    @info "set_model_geometry $(PB.fullname(rj))"

    @info "  reading 2D grid from $(rj.pars.grid_file.v)"
    grid2D = PB.Grids.CartesianGrid(
        rj.pars.grid_type.v,
        rj.pars.grid_file.v, rj.pars.coordinate_names.v,
        equalspacededges=rj.pars.equalspacededges.v
    )

    if rj.pars.grid_type.v == PB.Grids.CartesianLinearGrid
        # define a linear index including every cell, in column-major order (first indices are consecutive in memory)
        v_i = vec([i for i=1:grid2D.dims[1], j=1:grid2D.dims[2]])
        v_j = vec([j for i=1:grid2D.dims[1], j=1:grid2D.dims[2]])
        PB.Grids.set_linear_index(grid2D, v_i, v_j)
    end
    
    rj.domain.grid = grid2D        
    
    @info "  set $(rj.domain.name) Domain size=$(grid2D.ncells) grid=$grid2D"
         
    return nothing
end

function PB.register_methods!(rj::ReactionGrid2DNetCDF)
    @info "register_methods! $(PB.fullname(rj))"

    grid_vars = [
        PB.VarPropStateIndep("Asurf",  "m^2",  "horizontal area of surface"),
        PB.VarPropScalarStateIndep("Asurf_total",  "m^2",  "total horizontal area of surface"),
    ]

    PB.add_method_setup!(
        rj,
        setup_grid_2DNetCDF,
        (PB.VarList_namedtuple(grid_vars),)
    )

    return nothing
end

function setup_grid_2DNetCDF(
    m::PB.ReactionMethod,
    (grid_vars, ),
    cellrange::PB.AbstractCellRange,
    deltat
)
    rj = m.reaction

    grid2D = rj.domain.grid

    length(cellrange.indices) == grid2D.ncells ||
        error("tiled cellrange not supported")

    if !isempty(rj.pars.area_var.v)
        NCDatasets.Dataset(rj.pars.grid_file.v) do ds
            area2D = Array(ds[rj.pars.area_var.v][:, :, 1])
            grid_vars.Asurf .= cartesian_to_internal(rj.domain.grid, area2D)
        end
    else        
        for i in cellrange.indices
            lon_edges = PB.Grids.get_lon_edges(grid2D, i)
            lat_edges = PB.Grids.get_lat_edges(grid2D, i)
            area = calc_spherical_area(rj.pars.planet_radius.v, lon_edges, lat_edges)
            grid_vars.Asurf[i] = area
        end
    end

    
    grid_vars.Asurf_total[] = sum(grid_vars.Asurf)

    return nothing
end

"area of cell on regular spherical grid"
function calc_spherical_area(r, (lond_l, lond_u), (latd_l, latd_u))
    Astrip = 2*pi*r^2*(sind(latd_u) - sind(latd_l)) # area of strip between latitudes
    A = (lond_u - lond_l)/360.0 * Astrip # area of cell
    return A
end


end # module
