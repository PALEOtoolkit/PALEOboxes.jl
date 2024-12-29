module GridReactions

import PALEOboxes as PB
using ..DocStrings

import Infiltrator # Julia debugger


"""
    ReactionUnstructuredVectorGrid

Create an [`PB.Grids.UnstructuredVectorGrid`](@ref) with `ncells` (from `ncells` Parameter).

# Parameters
$(PARS)
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

    grid = PB.Grids.UnstructuredVectorGrid(ncells=rj.pars.ncells[])
    rj.domain.grid = grid
    
    @info "  set $(rj.domain.name) Domain size=$(grid.ncells) grid=$grid"
         
    return nothing
end

PB.register_methods!(rj::ReactionUnstructuredVectorGrid) = nothing


"""
    ReactionCartesianGrid

Create a [`PB.Grids.CartesianArrayGrid`](@ref) with `dims` and `dimnames``

# Parameters
$(PARS)
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
        [
            PB.NamedDimension(name, sz) 
            for (name, sz) in PB.IteratorUtils.zipstrict(
                rj.pars.dimnames.v, rj.pars.dims.v; 
                errmsg="parameters 'dimnames' and 'dims' are not the same length"
            )
        ]
    )
    
    rj.domain.grid = grid
    
    @info "  set $(rj.domain.name) Domain size=$(grid.ncells) grid=$grid"
         
    return nothing
end

PB.register_methods!(rj::ReactionCartesianGrid) = nothing


"""
    ReactionGrid2DNetCDF

Create a 2D [`PB.Grids.CartesianLinearGrid`](@ref) from grid information in a NetCDF file.

# Parameters
$(PARS)

# Methods
$(METHODS_SETUP)
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
   
    coord_varnames_values = nothing
    coord_edges_varnames_values = nothing
end
   
function PB.set_model_geometry(rj::ReactionGrid2DNetCDF, model::PB.Model)
    @info "set_model_geometry $(PB.fullname(rj))"

    @info "  reading 2D grid from $(rj.pars.grid_file[])"
    grid2D, rj.coord_varnames_values, rj.coord_edges_varnames_values = PB.Grids.CartesianGrid(
        rj.pars.grid_type[],
        rj.pars.grid_file[], rj.pars.coordinate_names.v,
        equalspacededges=rj.pars.equalspacededges[]
    )

    if rj.pars.grid_type[] == PB.Grids.CartesianLinearGrid
        # define a linear index including every cell, in column-major order (first indices are consecutive in memory)
        v_i = vec([i for i=1:grid2D.dimensions[1].size, j=1:grid2D.dimensions[2].size])
        v_j = vec([j for i=1:grid2D.dimensions[1].size, j=1:grid2D.dimensions[2].size])
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

    coord_vars = [
        PB.VarPropScalarStateIndep(coord_name,  "",  "coordinate variable"; 
            attributes=(:field_data=>PB.ArrayScalarData, :data_dims=>(coord_name,), ))
        for (coord_name, _) in rj.coord_varnames_values
    ]

    coord_edges_vars = [
        PB.VarPropScalarStateIndep(coord_name,  "",  "coordinate edge variable";
            attributes=(:field_data=>PB.ArrayScalarData, :data_dims=>(coord_name,), ))
        for (coord_name, _) in rj.coord_edges_varnames_values
    ]

    PB.add_method_setup!(
        rj,
        setup_grid_2DNetCDF,
        (PB.VarList_namedtuple(grid_vars), PB.VarList_vector(coord_vars), PB.VarList_vector(coord_edges_vars))
    )

    return nothing
end

function setup_grid_2DNetCDF(
    m::PB.ReactionMethod,
    pars,
    (grid_vars, coord_vars, coord_edges_vars),
    cellrange::PB.AbstractCellRange,
    attribute_name
)
    rj = m.reaction

    attribute_name == :setup || return

    @info "$(PB.fullname(m)):"

    grid2D = rj.domain.grid

    length(cellrange.indices) == grid2D.ncells ||
        error("tiled cellrange not supported")

    for (coord_var, (coord_name, coord_values)) in PB.IteratorUtils.zipstrict(coord_vars, rj.coord_varnames_values)
        coord_var .= coord_values
    end

    for (coord_edge_var, (coord_edge_name, coord_edge_values)) in PB.IteratorUtils.zipstrict(coord_edges_vars, rj.coord_edges_varnames_values)
        coord_edge_var .= coord_edge_values
    end

    if !isempty(pars.area_var[])
        NCDatasets.Dataset(pars.grid_file[]) do ds
            area2D = Array(ds[pars.area_var[]][:, :, 1])
            grid_vars.Asurf .= PB.Grids.cartesian_to_internal(rj.domain.grid, area2D)
        end
    elseif !isempty(coord_edges_vars)
        lonedgevar, latedgevar = coord_edges_vars[[grid2D.londim, grid2D.latdim]]
        for idx in cellrange.indices
            lonidx = PB.Grids.get_lon_idx(grid2D, idx)
            latidx = PB.Grids.get_lat_idx(grid2D, idx)
            lon_edges = lonedgevar[lonidx:lonidx+1]
            lat_edges = latedgevar[latidx:latidx+1]
            area = calc_spherical_area(pars.planet_radius[], lon_edges, lat_edges)
            grid_vars.Asurf[idx] = area
        end
    else
        @warn "no area_var or coords edges specified, not calculating Asurf"
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
