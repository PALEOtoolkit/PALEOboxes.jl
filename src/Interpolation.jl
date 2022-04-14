"""
    LinInterp(xvals, xperiod=0.0; require_sorted_input=true) -> LinInterp

Create a `LinInterp` struct for generic linear interpolation for (scalar) value `xval` into a set of records at `xvals`.

[`interp`](@ref) calculates weights and indices for the two records enclosing `xval`,
hence allows use for arbitrary record types.

# Args:
    - `xvals::Vector`: x values at which function is available.
    - `xperiod`: Periodicity of x values (0.0 for not periodic)
    - `require_sorted_input::Bool`: `true` to check that input is sorted in ascending order, 
     `false` to allow arbitrary order (which will then be sorted for internal use).
"""
struct LinInterp
    "supplied values at each record index"
    xvals::Vector{Float64}  
    "optional periodicity (0.0 means not periodic)"
    xperiod::Float64
    "modified values, mapped to `xperiod`, sorted in order, and padded at each end"
    _xvals::Vector{Float64}
    "record indices corresponding to _xvals"
    _recidx::Vector{Int}
    "out-of-range values set to constant"
    _extrap_const::Bool

    function LinInterp(xvals, xperiod=0.0; require_sorted_input=true, extrap_const=false)
        if require_sorted_input
            issorted(xvals) || error("`xvals` not sorted and `require_sorted_input==true`")
        end

        if xperiod > 0.0 && xperiod < Inf   
            # periodic - normalize, sort, and pad
            normxvals = mod.(xvals, Float64(xperiod)) # all normxvals now in range 0 - xperiod
            p = sortperm(normxvals) # if xvals contained eg out-of-range values, normxvals will now be out of ascending order
            _xvals = zeros(length(normxvals)+2)
            _xvals[2:end-1] = normxvals[p]
            _recidx = zeros(Int, length(normxvals)+2)
            _recidx[2:end-1] = p
            # pad (these values will be outside range 0 - xperiod)
            _xvals[1] = _xvals[end-1] - xperiod
            _xvals[end] = _xvals[2] + xperiod
            _recidx[1] = _recidx[end-1]
            _recidx[end] = _recidx[2]
            # check sorted (failure not expected here)
            issorted(_xvals) || error("programming error: normalized periodic xvalues not sorted xvals=$xvals, _xvals=$_xvals")
        elseif xperiod == 0.0
            p = sortperm(xvals)
            _xvals = xvals[p]
            _recidx = p       
        else
            error("invalid `xperiod` $xperiod")
        end

        return new(copy(xvals), xperiod, _xvals, _recidx, extrap_const)
    end
end

"""
    interp(lininterp::LinInterp, xval) -> ((wt_lo, recidx_lo), (wt_hi, recidx_hi))

Linear interpolation for (scalar) value `xval`. See [`LinInterp`](@ref).
"""
function interp(lininterp::LinInterp, xval)
        
    if lininterp.xperiod != 0.0
        # AD workaround (ForwardDiff generates NaN for derivative if on mod boundary ?)
        # as we know the solution is periodic.
        xval_val = mod(value_ad(xval), lininterp.xperiod) # sign of divisor ie +ve
        xval += xval_val - value_ad(xval)
    end

    # Return the index of the last value in a less than or equal to x
    idx = searchsortedlast(lininterp._xvals, xval)

    # exact equality to last record means we need the previous record (as we take a linear combination of idx and idx+1)
    if idx==length(lininterp._xvals) && xval == lininterp._xvals[end]
        idx -= 1
    end

    if idx == 0 
        if lininterp._extrap_const
            idx = 1
            xval = lininterp._xvals[1]
        else
            error("xval out of range")
        end
    end

    if idx == length(lininterp._xvals)
        if lininterp._extrap_const
            idx = length(lininterp._xvals) - 1
            xval = lininterp._xvals[end]
        else
            error("xval out of range")
        end
    end

    x_lo = lininterp._xvals[idx]
    x_hi = lininterp._xvals[idx+1]
    rec_lo = lininterp._recidx[idx]
    rec_hi = lininterp._recidx[idx+1]

    wt_lo = (x_hi - xval)/(x_hi - x_lo)
    wt_hi = (xval - x_lo)/(x_hi - x_lo)

    return ((wt_lo, rec_lo), (wt_hi, rec_hi))
end

"""
    interp(lininterp::LinInterp, xval, yvals) -> yval

Linear interpolation for (scalar) value `xval`. See [`LinInterp`](@ref).
"""
function interp(lininterp::LinInterp, xval, yvals)
    (wt_lo, rec_lo), (wt_hi, rec_hi) = interp(lininterp, xval)

    yval = wt_lo*yvals[rec_lo] + wt_hi*yvals[rec_hi]

    return yval
end
