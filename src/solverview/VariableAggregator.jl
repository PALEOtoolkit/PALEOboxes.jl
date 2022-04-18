
# import Infiltrator

"""
    VariableAggregator{T, F, C}

Aggregate multiple VariableDomains into a contiguous Vector for use by a numerical solver.
    
Access as a contiguous Vector using `copyto!(va::VariableAggregator, x)` and `copyto!(x, va::VariableAggregator)`
"""
mutable struct VariableAggregator{T, F <: Tuple, C <: Tuple}
    # Variables 
    vars::Vector{VariableDomain}

    # indices in contiguous Vector for each Variable
    indices::Vector{UnitRange{Int64}}

    # field and corresponding CellRange for each Variable
    fields::F
    cellranges::C
end

"""
    VariableAggregator(vars, cellranges, modeldata)

Create a [`VariableAggregator`](@ref) for collection of Variables `vars`,
with indices from corresponding `cellranges`, for `modeldata`.

`cellranges` may contain `nothing` entries to indicate whole Domain.
"""
function VariableAggregator(vars, cellranges, modeldata)
    length(vars) == length(cellranges) || 
        throw(ArgumentError("'vars' and 'cellranges' must be of same length"))

    fields = []
    indices = UnitRange{Int64}[]
    nextidx = 1
    for (v, cr) in zip(vars, cellranges)
        f = get_field(v, modeldata)

        dof = dof_field(f, cr)

        r = nextidx:(nextidx+dof-1)
    
        push!(fields, f)
        push!(indices, r)

        nextidx += dof
    end

    fields=Tuple(fields)
    cellranges=Tuple(cellranges)

    return VariableAggregator{eltype(modeldata), typeof(fields), typeof(cellranges)}(
        vars, indices, fields, cellranges
    )
end


"compact form"
function Base.show(io::IO, va::VariableAggregator)
    print(io, "VariableAggregator(length=$(length(va)), number of variables=$(length(va.vars)))")
    return nothing
end

"multiline form"
function Base.show(io::IO, ::MIME"text/plain", va::VariableAggregator)
    println(io, "VariableAggregator(length=$(length(va)), number of variables=$(length(va.vars))):")
    if length(va) > 0
        println(
            io,
            "    ",
            rpad("var", 6),
            rpad("indices", 14),
            rpad("name", 40),
            rpad("cellrange.indices", 14)
        )
        for i in eachindex(va.vars)
            println(
                io,
                "    ",
                rpad(i,6),
                rpad(va.indices[i], 14),
                rpad("$(fullname(va.vars[i]))", 40),
                rpad(isnothing(va.cellranges[i]) ? "-" : va.cellranges[i].indices, 14),
            )
        end
    end
    return nothing
end

Base.eltype(va::VariableAggregator{T, F, C}) where{T, F, C} = T

function Base.length(va::VariableAggregator)
    if isempty(va.indices)
        return 0
    else
        return last(last(va.indices))
    end
end


"""
    copyto!(dest::VariableAggregator, src::AbstractVector; sof=1) -> num_copied::Int

Set aggregated Variables `dest` = (contiguous) Vector `src`.

Optional `sof` sets first index in `src` to use.
"""
function Base.copyto!(dest::VariableAggregator, src::AbstractVector; sof::Int=1)
   
    function copy(si, field, cellrange)
        si += copyto!(field, cellrange, src, si)
        return si
    end

    fof = IteratorUtils.reduce_longtuple(copy, sof, dest.fields, dest.cellranges)

    return fof - sof
end

"""
    copyto!(dest::AbstractVector, va::VariableAggregator; dof=1) -> num_copied::Int

Set (contiguous) Vector `dest` = aggregated Variables `src`

Optional `dof` sets first index in `dest`
"""
function Base.copyto!(dest::AbstractVector, src::VariableAggregator; dof::Int=1)


    function copy(di, field, cellrange)
        di += copyto!(dest, di, field, cellrange)
        return di
    end

    fof = IteratorUtils.reduce_longtuple(copy, dof, src.fields, src.cellranges)

    return fof - dof
end

"""
    get_data(va::VariableAggregator) -> Vector

Allocate Vector and set to values of aggregated Variables `va`.
"""
function get_data(va::VariableAggregator)
    x = Vector{eltype(va)}(undef, length(va))
    copyto!(x, va)
    return x
end

"vay += a * vax"
function add_data!(vay::VariableAggregator{T, F, C}, a, vax::VariableAggregator{T, F, C}) where {T, F, C}
    length(vay) == length(vax) ||
        throw(ArgumentError("'vay' and 'vax' must be of same length"))
 
    adapt_add_field!(dest, cellrange, src) = add_field!(dest, a, cellrange, src)

    IteratorUtils.foreach_longtuple(adapt_add_field!, vay.fields, vay.cellranges, vax.fields)

    return nothing
end

"vay += a * vax"
function add_data!(vay::VariableAggregator{T, F, C}, a, vax::AbstractArray) where {T, F, C}
    length(vay) == length(vax) ||
        throw(ArgumentError("'vay' and 'vax' must be of same length"))
 
    function adapt_add_field_vec!(vidx, dest, cellrange)
        vidx += add_field_vec!(dest, a, cellrange, vax, vidx)
        return vidx
    end
    IteratorUtils.reduce_longtuple(adapt_add_field_vec!, 1, vay.fields, vay.cellranges)

    return nothing
end

"get Variables"
get_vars(va::VariableAggregator) = va.vars
   
num_vars(va::VariableAggregator) = length(va.vars)


"""
    VariableAggregatorNamed

Aggregate multiple VariableDomains into a NamedTuple, for easy access to 
host-dependent Variables by name.
"""

mutable struct VariableAggregatorNamed{V <: NamedTuple}
    # Variables 
    vars::Vector{VariableDomain}

    # Values for each Variable as NamedTuple
    values::V
end

Base.length(va::VariableAggregatorNamed) = length(va.vars)

"compact form"
function Base.show(io::IO, va::VariableAggregatorNamed)
    print(io, "VariableAggregatorNamed($([fullname(v) for v in va.vars]))")
    return nothing
end

function VariableAggregatorNamed(
    vars, modeldata::AbstractModelData;
    reallocate_hostdep_eltype=Float64,
)

    # If requested, change data type eg to remove AD type
    if !isnothing(reallocate_hostdep_eltype)
        for v in vars
            v_data = get_data(v, modeldata)
            if v_data isa AbstractArray && eltype(v_data) != reallocate_hostdep_eltype
                @info "VariableAggregatorNamed: reallocate $(fullname(v)) data $(eltype(v_data)) -> $reallocate_hostdep_eltype"
                set_data!(v, modeldata, similar(v_data, reallocate_hostdep_eltype))
            end
        end
    end

    # sort by Domain into Dict of Dicts
    vars_domains = Dict()
    for v in vars
        d = get!(vars_domains, Symbol(v.domain.name), Dict())
        d[Symbol(v.name)] = get_data(v, modeldata)
    end

    # convert to Dict of NamedTuples
    vars_domains_nt = Dict()
    for (dname, dvars) in vars_domains
        vars_domains_nt[dname] = NamedTuple{Tuple(keys(dvars))}(Tuple(values(dvars)))
    end
    
    vars_values = NamedTuple{Tuple(keys(vars_domains_nt))}(Tuple(values(vars_domains_nt)))

    return VariableAggregatorNamed(copy(vars), vars_values)
end

"""
    set_values!(van::VariableAggregatorNamed, domainname, varname, val; allow_missing=false)

Set a data array in `van` to `val`.
"""
function set_values!(
    van::VariableAggregatorNamed, domainname::AbstractString, varname::AbstractString, val; 
    allow_missing::Bool=false
)
    return set_values!(van, Symbol(domainname), Symbol(varname), val, allow_missing=allow_missing)
end

function set_values!(
    van::VariableAggregatorNamed, domainname::Symbol, varname::Symbol, val;
    allow_missing::Bool=false
)
    
    if hasfield(typeof(van.values), domainname)
        dvalues = getfield(van.values, domainname) 
        if hasfield(typeof(dvalues), varname)
            vvalues = getfield(dvalues, varname)
            vvalues .= val
        else
            allow_missing ||
                error("VariableAggregatorNamed has no Variable $domainame.$varname")
        end
    else
        allow_missing ||
            error("VariableAggregatorNamed has no Variable $domainame.$varname")
    end

    return nothing
end

"""
    set_tforce!(van::VariableAggregatorNamed, val; allow_missing::Bool=true)

Optimised non-allocating version of `set_values!(van, :global, :tforce, val)`
"""
function set_tforce!(
    van::VariableAggregatorNamed, val;
    allow_missing::Bool=true
)
    
    if hasfield(typeof(van.values), :global)
        dvalues = getfield(van.values, :global) 
        if hasfield(typeof(dvalues), :tforce)
            vvalues = getfield(dvalues, :tforce)
            vvalues .= val
        else
            allow_missing ||
                error("VariableAggregatorNamed has no Variable global.tforce")
        end
    else
        allow_missing ||
            error("VariableAggregatorNamed has no Variable global.tforce")
    end

    return nothing
end

