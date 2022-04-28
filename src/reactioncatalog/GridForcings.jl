module GridForcings

import NCDatasets
import MAT
import Interpolations

import PALEOboxes as PB
using ..DocStrings

using Infiltrator # Julia debugger

"""
    ReactionForceGrid
 
Apply time-dependent, periodic forcing from a variable in a netcdf or Matlab file, optionally applying a scaling, a constant linear offset,
and a linear offset generated from a scalar model variable.

Reads records `tidx_start`:`tidx_end` (assumed to be the last dimension) for `data_var` from a gridded dataset in `netcdf_file`
or `matlab_file`, maps grid to linear using the Domain grid (which must match that of `data_var`).

If `tidx_end > tidx_start` (ie multiple records), reads a netcdf or Matlab variable named `time_var`, applies periodicity `cycle_time`, and uses this to linearly interpolate to model time.

Then applies forcing:

    F = scale*`data_var` + constant_offset + scale_offset_var*`scalar_offset_var`

Optionally (if `interp_vars` is non-empty), interpolate forcing from additional dimensions in the netcdf file,
given values supplied by additional Variable dependencies.
NB: netcdf dimensions order must be `grid_vars` x `interp_vars`` x `time_var`, where order within `interp_vars` also must match netcdf order.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForceGrid{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParString("netcdf_file", "",
            description="netcdf file with gridded time-series data"),
        PB.ParString("matlab_file", "",
            description="matlab file with gridded time-series data"),

        PB.ParString("data_var", "",
            description="variable name in data file"),
        PB.ParString("time_var", "time",
            description="time variable name in data file (empty to generate evenly spaced times from cycle_time)"),

        PB.ParInt("tidx_start", 1,
            description="first record in data file to use"),
        PB.ParInt("tidx_end", 1,
            description="last record in data file to use (set equal to tidx_start for constant forcing from single record)"),
        PB.ParBool("use_timeav", false,
            description="true to average records and provide constant forcing at time-averaged value"),
        PB.ParDouble("cycle_time", 1.0, units="yr",
            description="time periodicity to apply (0.0 to disable periodic)"),

        PB.ParStringVec("interp_vars", String[],
            description="optional interpolation variables for additional grid dimensions. "*
                "NB: netcdf dimensions order must be grid_vars x interp_vars x time_var, where order within interp_vars also must match netcdf order"),
        PB.ParBoolVec("interp_log", Bool[],
            description="true to interpolate interp_vars in log space"),

        PB.ParDouble("scale", 1.0, units="",
            description="scaling factor to apply"),
        PB.ParDouble("constant_offset", 0.0, units="",
            description="constant offset to apply"),
        PB.ParDouble("scale_offset_var", 0.0, units="",
            description="scaling for additional scalar offset from model variable (0.0 to disable)"),
    )

 
    data_time           = nothing
    data_var            = nothing

    time_interp         = nothing
    interp_interp       = nothing
    interp_fn           = nothing
end

function PB.register_methods!(rj::ReactionForceGrid)

    @info "register_methods! $(PB.fullname(rj))"

    vars = [
        PB.VarDepScalar("global.tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarProp("F", "", "interpolated forcing"),    
    ]
    if rj.pars.scale_offset_var.v != 0.0
        @info "  adding scalar offset from Variable 'scalar_offset'"
        push!(vars, PB.VarDepScalar("scalar_offset", "",  "scalar offset"))
    end
    PB.setfrozen!(rj.pars.scale_offset_var)

    length(rj.pars.interp_vars.v) == length(rj.pars.interp_log.v) ||
        error("length(interp_vars) != length(interp_log)")
    interp_vars = []
    for (vname, vlog) in zip(rj.pars.interp_vars.v, rj.pars.interp_log.v)
        @info "  adding interpolation Variable $vname log $vlog"
        push!(interp_vars, PB.VarDepScalar(vname, "",  "interpolation variable log $vlog"))
    end
    PB.setfrozen!(rj.pars.interp_vars, rj.pars.interp_log)

    PB.add_method_do!(
        rj,
        do_force_grid,
        (PB.VarList_namedtuple(vars), PB.VarList_tuple(interp_vars)),
        preparefn=prepare_do_force_grid,
    )

    return nothing
end

function prepare_do_force_grid(
    m::PB.ReactionMethod, 
    (vars, interp_vars),
)
    rj = m.reaction
   
    @info "prepare_do_force_grid: $(PB.fullname(rj))"

    if !isempty(rj.pars.netcdf_file.v) && isempty(rj.pars.matlab_file.v)
        @info "    reading variable '$(rj.pars.data_var.v)' from netcdf file '$(rj.pars.netcdf_file.v)'"

        NCDatasets.Dataset(rj.pars.netcdf_file.v) do ds
            _prepare_data(rj, ds)
        end
    elseif !isempty(rj.pars.matlab_file.v) && isempty(rj.pars.netcdf_file.v)
        @info "    reading variable '$(rj.pars.data_var.v)' from matlab file '$(rj.pars.matlab_file.v)'"

        ds = MAT.matread(rj.pars.matlab_file.v) # must return a Dict varname=>vardata
        _prepare_data(rj, ds)
    else
        error("    both netcdf_file $(rj.pars.netcdf_file.v) and matlab_file $(rj.pars.matlab_file.v) are specified")
    end
 
    return (vars, interp_vars, rj.time_interp, rj.interp_interp, rj.interp_fn, rj.data_var)
end

"""
    _prepare_data(rj::ReactionForceGrid, ds)

Read data variable from `ds[rj.pars.data_var.v]` and optional time values from `ds[rj.pars.time_var.v]`.  
Map data variable to PALEO internal array layout and store in `rj.data_var`, generate time interpolator in `rj.time_interp`.
"""
function _prepare_data(rj::ReactionForceGrid, ds)
        
    num_nc_time_recs = rj.pars.tidx_end.v - rj.pars.tidx_start.v + 1

    if rj.pars.use_timeav.v 
        num_time_recs = 1
    else
        num_time_recs = num_nc_time_recs
    end

    # NB: A PALEO Cartesian Grid may define an internal_size for model Variables that has different dimensions (eg a linear Vector)
    # to the n-D cartesian_size of the forcings read from the NetCDF file 
    ncartesiandims = length(PB.cartesian_size(rj.domain.grid))  # as read from NetCDF
    cartesiancolons = fill(Colon(), ncartesiandims)
    ninternaldims = length(PB.internal_size(rj.domain.grid)) # PALEO array layout
    internalcolons = fill(Colon(), ninternaldims)

    ninterpdims = length(rj.pars.interp_vars.v)
    interpcolons = fill(Colon(), ninterpdims)
    interpdims = [length(ds[vn]) for vn in rj.pars.interp_vars.v]

    # map to grid internal storage (ie mapping ncartesiandims -> ninternaldims)
    #           grid             interp           time
    data_ndims = ninternaldims + length(interpdims) + 1
    rj.data_var = Array{Float64, data_ndims}(undef, PB.internal_size(rj.domain.grid)..., interpdims..., num_time_recs)
    @info "  size(data_var) = $(size(rj.data_var))"
        
    # TODO - reorder indices (currently require  gridvars..., interpvars..., timevar)
    # read netcdf data, taking only the time records we need
    tmp_var  = Array(ds[rj.pars.data_var.v][cartesiancolons..., interpcolons..., rj.pars.tidx_start.v:rj.pars.tidx_end.v])

    # copy into data_var,  mapping ngriddims -> 1 linear index and creating time average if necessary
    if rj.pars.use_timeav.v
        @views rj.data_var[internalcolons..., interpcolons..., 1] .= 0.0
        for i in 1:num_nc_time_recs
            if length(interpdims) == 0
                @views rj.data_var[internalcolons..., 1] .+= PB.Grids.cartesian_to_internal(rj.domain.grid, tmp_var[cartesiancolons..., i])/num_nc_time_recs
            elseif length(interpdims) == 1
                for id1 in 1:interpdims[1]
                    @views rj.data_var[internalcolons..., id1, 1] .+= PB.Grids.cartesian_to_internal(rj.domain.grid, tmp_var[cartesiancolons..., id1, i])/num_nc_time_recs
                end
            elseif length(interpdims) == 2
                for id1 in 1:interpdims[1], id2 in 1:interpdims[2]
                    @views rj.data_var[internalcolons..., id1, id2, 1] .+= PB.Grids.cartesian_to_internal(rj.domain.grid, tmp_var[cartesiancolons..., id1, id2, i])/num_nc_time_recs
                end
            else
                error("  > 2 interpdims not supported")
            end
        end
    else
        for i in 1:num_nc_time_recs 
            if length(interpdims) == 0           
                @views rj.data_var[internalcolons..., i] .= PB.Grids.cartesian_to_internal(rj.domain.grid, tmp_var[cartesiancolons..., i])
            elseif length(interpdims) == 1
                for id1 in 1:interpdims[1]
                    @views rj.data_var[internalcolons..., id1, i] .= PB.Grids.cartesian_to_internal(rj.domain.grid, tmp_var[cartesiancolons..., id1, i])
                end
            elseif length(interpdims) == 2
                for id1 in 1:interpdims[1], id2 in 1:interpdims[2]
                    @views rj.data_var[internalcolons..., id1, id2, i] .= PB.Grids.cartesian_to_internal(rj.domain.grid, tmp_var[cartesiancolons..., id1, id2, i])
                end
            else
                error("  > 2 interpdims not supported")
            end
        end
    end

    # create time interpolator
    if num_time_recs > 1
        if !isempty(rj.pars.time_var.v)
            @info "  reading variable '$(rj.pars.time_var.v)' to generate time-dependent forcing from $num_time_recs time records"
            rj.data_time = Array(ds[rj.pars.time_var.v][rj.pars.tidx_start.v:rj.pars.tidx_end.v])
        else
            @info "  generating evenly spaced intervals to apply time-dependent forcing from $num_time_recs time records"
            rj.data_time = collect(range(0.5/num_time_recs, step=1.0/num_time_recs, length=num_time_recs))
        end
        # create time interpolation
        rj.time_interp = PB.LinInterp(rj.data_time, rj.pars.cycle_time.v)
    else
        @info "  generating constant time forcing from $num_nc_time_recs time record(s)"
        rj.time_interp = nothing
    end

    # create variable interpolator(s) (if any)
    rj.interp_interp = []
    rj.interp_fn = []
    for (vidx, (vname, vlog)) in enumerate(zip(rj.pars.interp_vars.v, rj.pars.interp_log.v))
        if vlog
            interp_fn = log
        else
            interp_fn = identity
        end
        push!(rj.interp_fn, interp_fn)
        vvalues = interp_fn.(Array(ds[vname]))
        @info "  interpolate $vname (log $vlog) values $vvalues"
        length(vvalues) == size(rj.data_var)[ninternaldims+vidx] ||
            error("  number of values != array size of dimension $(ninternaldims+vidx)")
        push!(rj.interp_interp, PB.LinInterp(vvalues, extrap_const=true))
    end
    rj.interp_fn = tuple(rj.interp_fn...)
    rj.interp_interp = tuple(rj.interp_interp...)

    return nothing
end

function do_force_grid(
    m::PB.ReactionMethod, 
    (vars, interp_vars, time_interp, interp_interp, interp_fn, data_var),
    cellrange::PB.AbstractCellRange,
    deltat,
)
    rj = m.reaction

    if !isnothing(time_interp)
        # time dependent forcing
        tforce = vars.tforce[]
    
        ((wt_lo, rec_lo), (wt_hi, rec_hi)) = PB.interp(time_interp, tforce)
    else
        # constant forcing
        wt_lo = 1.0
        rec_lo = 1
        wt_hi = 0.0
        rec_hi = 1
    end

    # apply scaling
    wt_lo *= rj.pars.scale.v
    wt_hi *= rj.pars.scale.v
    # apply scalar offset
    scalar_offset = rj.pars.constant_offset.v
    if rj.pars.scale_offset_var.v != 0.0
        scalar_offset += rj.pars.scale_offset_var.v * vars.scalar_offset[]
    end

    if length(interp_vars) == 0
        _do_interp_0(vars, data_var, cellrange, wt_lo, wt_hi, rec_lo, rec_hi, scalar_offset)
    elseif length(interp_vars) == 1
        _do_interp_1(vars, data_var, cellrange, wt_lo, wt_hi, rec_lo, rec_hi, scalar_offset, interp_interp, interp_fn, interp_vars)
    elseif length(interp_vars) == 2
        _do_interp_2(vars, data_var, cellrange, wt_lo, wt_hi, rec_lo, rec_hi, scalar_offset, interp_interp, interp_fn, interp_vars)
    else
        error("  > 2 interp_vars not supported")
    end

    return nothing
end

function _do_interp_0(vars, data_var, cellrange, wt_lo, wt_hi, rec_lo, rec_hi, scalar_offset)
    @inbounds for i in cellrange.indices
        vars.F[i] = wt_lo*data_var[i, rec_lo] + wt_hi*data_var[i, rec_hi] + scalar_offset
    end
end

function _do_interp_1(vars, data_var, cellrange, wt_lo, wt_hi, rec_lo, rec_hi, scalar_offset, interp_interp, interp_fn, interp_vars)
    ((wt1_lo, rec1_lo), (wt1_hi, rec1_hi)) = PB.interp(interp_interp[1], interp_fn[1](interp_vars[1][]))
    @inbounds for i in cellrange.indices
        vars.F[i] = (
            wt_lo*(wt1_lo*data_var[i, rec1_lo, rec_lo] + wt1_hi*data_var[i, rec1_hi, rec_lo])
          + wt_hi*(wt1_lo*data_var[i, rec1_lo, rec_hi] + wt1_hi*data_var[i, rec1_hi, rec_hi])
          + scalar_offset
        )
    end
end

function _do_interp_2(vars, data_var, cellrange, wt_lo, wt_hi, rec_lo, rec_hi, scalar_offset, interp_interp, interp_fn, interp_vars)
    ((wt1_lo, rec1_lo), (wt1_hi, rec1_hi)) = PB.interp(interp_interp[1], interp_fn[1](interp_vars[1][]))
    ((wt2_lo, rec2_lo), (wt2_hi, rec2_hi)) = PB.interp(interp_interp[2], interp_fn[2](interp_vars[2][]))
   
    @inbounds for i in cellrange.indices
        vars.F[i] = (
            wt_lo*(
                wt1_lo*(wt2_lo*data_var[i, rec1_lo, rec2_lo, rec_lo] + wt2_hi*data_var[i, rec1_lo, rec2_hi, rec_lo])
              + wt1_hi*(wt2_lo*data_var[i, rec1_hi, rec2_lo, rec_lo] + wt2_hi*data_var[i, rec1_hi, rec2_hi, rec_lo])
            )
          + wt_hi*(
                wt1_lo*(wt2_lo*data_var[i, rec1_lo, rec2_lo, rec_hi] + wt2_hi*data_var[i, rec1_lo, rec2_hi, rec_hi])
              + wt1_hi*(wt2_lo*data_var[i, rec1_hi, rec2_lo, rec_hi] + wt2_hi*data_var[i, rec1_hi, rec2_hi, rec_hi])
            )
          + scalar_offset              
        )

    end
end

"""
    insolMITgcmDIC(Timeyr,latdeg; albedo=0.6, solar=1360.0, parfrac=1.0) -> sfac

MITgcm DIC package insol function directly translated from fortran.
Similar to [Brock1981](@cite).

NB: there are three normalization constants here: `solar`, `albedo`, `parfrac` attempting to define
top-of-atmosphere flux (from astronomical formulae) -> a crude approx to ground level flux (taking into account clouds etc) -> photosynthetic PAR flux,
where the second two would seem to have been poorly defined or even conflated in the MITgcm literature.

    C !DESCRIPTION:
    C find light as function of date and latitude
    C based on paltridge and parson

# Arguments:
- `Timeyr`: yr, model time, NB: year assumed to start in winter
- `latdeg`: deg, latitudes
- `albedo`: planetary albedo (ie correct for top-of-atmosphere to ground-level, clouds etc)
- `solar`: W m-2 solar constant
- `parfrac`: photosynthetically active fraction

# Returns:
- `sfac`: daily average photosynthetically active solar radiation just below surface
"""
function insolMITgcmDIC(Timeyr,latdeg; albedo=0.6, solar=1360.0, parfrac=1.0)
  
    # C find day (****NOTE for year starting in winter*****)
    dayfrac= Timeyr - floor(Timeyr) # fraction of year
    yday = 2*π*dayfrac                    # convert to radians
    delta = (0.006918                
        -(0.399912 *cos(yday))           # cosine zenith angle
        +(0.070257 *sin(yday))           # (paltridge+platt)
        -(0.006758 *cos(2*yday))             
        +(0.000907 *sin(2*yday))
        -(0.002697 *cos(3*yday))
        +(0.001480 *sin(3*yday)) )

    # C latitude in radians
    lat=latdeg*2*π/360

    sun1 = -sin(delta)/cos(delta) * sin(lat)/cos(lat)
    sun1 = max(sun1,-0.999)
    sun1 = min(sun1, 0.999 )
    dayhrs = abs(acos(sun1))
    cosz = ( sin(delta)*sin(lat)+             # average zenith angle
         (cos(delta)*cos(lat).*sin(dayhrs)./dayhrs) )
    cosz = max(cosz, 5e-3)
    frac = dayhrs/π                           # fraction of daylight in day
    # C daily average photosynthetically active solar radiation just below surface
    fluxi = solar*(1-albedo)*cosz*frac*parfrac

    # C convert to sfac
    sfac = max(1e-5,fluxi)

    return sfac
end


"""
    ReactionForceInsolation
 
Calculate time and latitude dependent daily mean surface solar insolation.

# Parameters
$(PARS)

# Methods and Variables for default Parameters
$(METHODS_DO)
"""
Base.@kwdef mutable struct ReactionForceInsolation{P} <: PB.AbstractReaction
    base::PB.ReactionBase

    pars::P = PB.ParametersTuple(
        PB.ParDouble("albedo", 0.6,
            description="mean planetary albedo"),
        PB.ParDouble("parfrac", 1.0,
            description="fraction of radiation that is photosynthetically active"),
        PB.ParDoubleVec("latitude", Float64[], units="degrees N",
            description="if non-empty, override grid latitude and set explicitly for each surface cell"),
    )
end

function PB.register_methods!(rj::ReactionForceInsolation)

    vars = [
        PB.VarDepScalar("global.tforce", "yr",  "historical time at which to apply forcings, present = 0 yr"),
        PB.VarProp("insolation", "W m-2", "daily mean surface insolation"),
    ]

    PB.add_method_do!(
        rj,
        do_force_insolation, 
        (PB.VarList_namedtuple(vars),),
        p=rj.domain.grid,  # add as context so fully typed
    )
end


function PB.check_configuration(rj::ReactionForceInsolation, model::PB.Model)
    configok = true
    if !isempty(rj.pars.latitude.v) && !isnothing(rj.domain.grid)
        if rj.domain.grid.ncells != length(rj.pars.latitude.v)
            @warn "check_configuration $(PB.fullname(rj))  length(latitude) parameter $(length(rj.pars.latitude.v)) != grid.ncells $(grid.ncells)"
            configok = false
        end
    end
    return configok
end

function do_force_insolation(m::PB.ReactionMethod, (vars, ), cellrange::PB.AbstractCellRange, deltat)
    rj = m.reaction
    grid = m.p
    
    tforce = PB.value_ad(vars.tforce[])

    @inbounds for i in cellrange.indices
        if isempty(rj.pars.latitude.v)
            lat = PB.Grids.get_lat(grid, i)
        else
            lat = rj.pars.latitude.v[i]
        end
        vars.insolation[i] = insolMITgcmDIC(tforce,lat, albedo=rj.pars.albedo.v, parfrac=rj.pars.parfrac.v)
    end

    return nothing
end


end # module
