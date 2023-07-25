
import Infiltrator

#################################################################
# VariableAggregator
#################################################################

struct VariableAggregator{T, F <: Tuple, C <: Tuple, PF <: Tuple, PC <: AbstractVector}
    # modeldata used for data arrays (for diagnostic output only)
    modeldata::AbstractModelData  # not typed
    arrays_idx::Int

    # Variables 
    vars::Vector{VariableDomain}

    # indices in contiguous Vector for each Variable
    indices::Vector{UnitRange{Int64}}

    # field and corresponding CellRange for each Variable
    fields::F
    cellranges::C

    # optional packed fields and corresponding cellranges
    packed_vars::Vector{VariableDomain}
    packed_fields::PF
    packed_cellranges::PC
    packed_indices::UnitRange{Int64}
end

"""
    VariableAggregator(vars, cellranges, modeldata, arrays_idx [, packed_vars] [, packed_cellranges]) -> VariableAggregator

Aggregate multiple VariableDomains into a flattened list (a contiguous Vector).

Creates a `VariableAggregator` for collection of Variables `vars`, with indices from corresponding `cellranges`,
for data arrays in `modeldata` with `arrays_idx`.

`cellranges` may contain `nothing` entries to indicate whole Domain.

Optionally add `packed_vars` which should be CellSpace variables from a single Domain with cells defined by `packed_cellranges`

This is mostly useful for aggregating state Variables, derivatives, etc to implement an interface to a generic ODE/DAE etc solver.

Values may be copied to and from a Vector using [`copyto!`](@ref)

"""
function VariableAggregator(vars, cellranges, modeldata::AbstractModelData, arrays_idx::Int, packed_vars=[], packed_cellranges=[])

    IteratorUtils.check_lengths_equal(vars, cellranges; errmsg="'vars' and 'cellranges' must be of same length")

    fields = []
    indices = UnitRange{Int64}[]
    nextidx = 1
    for (v, cr) in IteratorUtils.zipstrict(vars, cellranges)
        f = get_field(v, modeldata, arrays_idx)

        dof = dof_field(f, cr)

        r = nextidx:(nextidx+dof-1)
    
        push!(fields, f)
        push!(indices, r)

        nextidx += dof
    end

    packed_fields = []
    packed_dof = 0
    for v in packed_vars
        f = get_field(v, modeldata, arrays_idx)    
        push!(packed_fields, f)
        for cr in packed_cellranges
            packed_dof += dof_field(f, cr)
        end
    end
    packed_indices = nextidx:(nextidx+packed_dof-1)

    fields = Tuple(fields)
    cellranges = Tuple(cellranges)
    packed_fields = Tuple(packed_fields)

    return VariableAggregator{eltype(modeldata, arrays_idx), typeof(fields), typeof(cellranges), typeof(packed_fields), typeof(packed_cellranges)}(
        modeldata, arrays_idx, copy(vars), indices, fields, cellranges, copy(packed_vars), packed_fields, packed_cellranges, packed_indices
    )
end


# compact form
function Base.show(io::IO, va::VariableAggregator)
    print(io, "VariableAggregator(modeldata=$(va.modeldata), arrays_idx=$(va.arrays_idx), length=$(length(va)), number of variables=$(length(va.vars)), number of packed variables=$(length(va.packed_vars)))")
    return nothing
end

# multiline form
function Base.show(io::IO, ::MIME"text/plain", va::VariableAggregator)
    println(io, "VariableAggregator(modeldata=$(va.modeldata), arrays_idx=$(va.arrays_idx), length=$(length(va)), number of variables=$(length(va.vars)), number of packed variables=$(length(va.packed_vars))):")
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
        
        if !isempty(va.packed_vars)
            println(
                io,
                "    ",
                rpad("var", 6),
                rpad("indices", 14),
                rpad("name", 40),
            )
            for i in eachindex(va.packed_vars)
                println(
                    io,
                    "    ",
                    rpad(i,6),
                    rpad(va.packed_indices, 14),
                    rpad("$(fullname(va.packed_vars[i]))", 40),
                )
            end
        end
    end
    return nothing
end

Base.eltype(va::VariableAggregator{T, F, C}) where{T, F, C} = T

Base.length(va::VariableAggregator) = last(last(va.packed_indices))

"""
    get_indices(va::VariableAggregator, varnamefull::AbstractString; allow_not_found=false) -> indices::UnitRange{Int64} 

Return indices in flattened Vector corresponding to Variable `varnamefull`
"""
function get_indices(va::VariableAggregator, varnamefull::AbstractString; allow_not_found=false)
    indices = nothing
    for (v, ind) in IteratorUtils.zipstrict(va.vars, va.indices)
        if fullname(v) == varnamefull
            indices = ind
            break
        end
    end

    !isnothing(indices) || allow_not_found || error("get_indices: Variable $varnamefull not found")
    return indices
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

    for cr in dest.packed_cellranges
        fof = IteratorUtils.reduce_longtuple_p(copy, fof, dest.packed_fields, cr)
    end

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

    for cr in src.packed_cellranges
        fof = IteratorUtils.reduce_longtuple_p(copy, fof, src.packed_fields, cr)
    end

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

#################################################################
# VariableAggregatorNamed
#################################################################


struct VariableAggregatorNamed{V <: NamedTuple, VV <: NamedTuple}
    # modeldata used for data arrays (for diagnostic output only)
    modeldata::AbstractModelData  # not typed
    arrays_idx::Int

    # Variables as a NamedTuple 
    vars::V

    # Values for each Variable as NamedTuple
    values::VV
end

function num_vars(va::VariableAggregatorNamed)
    nv = 0
    for vnt in va.vars
        nv += length(vnt)
    end
    return nv
end

function Base.length(va::VariableAggregatorNamed)
    l = 0
    for vnt in va.values
        for v in vnt
            l += length(v)
        end
    end
    return l
end    

# compact form
function Base.show(io::IO, van::VariableAggregatorNamed)
    print(io, "VariableAggregatorNamed(modeldata=$(van.modeldata), arrays_idx=$(van.arrays_idx), ")
    for (domainname, dvarsnt) in zip(keys(van.vars), van.vars)
        print(io, "$domainname.*$([k for k in keys(dvarsnt)]), ")
    end
    print(io, ")")
    return nothing
end

# multiline form
function Base.show(io::IO, ::MIME"text/plain", van::VariableAggregatorNamed)
    println(io, "VariableAggregatorNamed(modeldata=$(van.modeldata), arrays_idx=$(van.arrays_idx):")
    for (domainname, dvaluesnt) in zip(keys(van.vars), van.values)
        println(io, "  $domainname:")
        for (varname, varvalues) in zip(keys(dvaluesnt), dvaluesnt)
            println(io, "    $(rpad(varname, 20))$varvalues")
        end
    end
    return nothing
end

"""
    VariableAggregatorNamed(modeldata, arrays_idx=1) -> VariableAggregatorNamed
    VariableAggregatorNamed(vars, modeldata, arrays_idx=1) -> VariableAggregatorNamed
    VariableAggregatorNamed(va::VariableAggregator; ignore_cellranges=false) -> VariableAggregatorNamed

Aggregate VariableDomains into nested NamedTuples, with Domain and Variable names as keys and
data arrays (from `get_data`) as values.

Any `/` characters in Variable names are replaced with `__` (double underscore)

This provides direct access to Variables by name, and is mostly useful for testing or for small models.

# Fields
- `vars`: nested NamedTuples (domainname, varname) of VariableDomains
- `values`: nested NamedTuples (domainname, varname) of data arrays.
"""
function VariableAggregatorNamed(
    modeldata::AbstractModelData, arrays_idx=1,
)
    # get all domain Variables
    all_domvars = Vector{VariableDomain}()
    for dom in modeldata.model.domains
        append!(all_domvars, get_variables(dom))
    end

    return VariableAggregatorNamed(all_domvars, modeldata, arrays_idx)
end

function VariableAggregatorNamed(va::VariableAggregator, ignore_cellranges=false)
    all(va.cellranges .== nothing) || ignore_cellranges || error("VariableAggregatorNamed: VariableAggregator is using cellranges to select part of Domain")
    return VariableAggregatorNamed(va.vars, va.modeldata, va.arrays_idx)
end

function VariableAggregatorNamed(
    vars, modeldata::AbstractModelData, arrays_idx=1,
)

    # sort into Dicts of
    #    domain.name => (Dict of v.name=>v)
    #    domain.name => (Dict of v.name=>data)

    # Sort by name and then use OrderedDict to preserve alphabetical order (domainname, varname)
    domains_vars = OrderedCollections.OrderedDict()
    domains_values = OrderedCollections.OrderedDict()
    for v in sort(vars, by=fullname)
        sym_vname = Symbol(replace(v.name, "/"=>"__")) # Reaction-private Variables use <reactionname>/<varname>, but / isn't a legal Symbol name
        d_vars = get!(domains_vars, Symbol(v.domain.name), OrderedCollections.OrderedDict())
        d_vars[sym_vname] = v
        d_values = get!(domains_values, Symbol(v.domain.name), OrderedCollections.OrderedDict())
        d_values[sym_vname] = get_data(v, modeldata, arrays_idx)
    end

    # convert to nested NamedTuples, via a Dict of NamedTuples
    vars_dictnt = OrderedCollections.OrderedDict()
    for (dname, d_vars) in domains_vars
        vars_dictnt[dname] = NamedTuple{Tuple(keys(d_vars))}(Tuple(values(d_vars)))
    end
    vars_nt = NamedTuple{Tuple(keys(vars_dictnt))}(Tuple(values(vars_dictnt)))

    vars_values_dictnt = OrderedCollections.OrderedDict()
    for (dname, d_values) in domains_values
        vars_values_dictnt[dname] = NamedTuple{Tuple(keys(d_values))}(Tuple(values(d_values)))
    end
    vars_values_nt = NamedTuple{Tuple(keys(vars_values_dictnt))}(Tuple(values(vars_values_dictnt)))

    return VariableAggregatorNamed(modeldata, arrays_idx, vars_nt, vars_values_nt)
end


"""
    set_values!(van::VariableAggregatorNamed, domainname_sym::Symbol, varname_sym::Symbol, val; allow_not_found=false)
    set_values!(van::VariableAggregatorNamed, domainname::AbstractString, varname::AbstractString, val; allow_not_found=false)
    set_values!(van::VariableAggregatorNamed, ::Val{domainname_sym}, ::Val{varname_sym}, val; allow_not_found=false)    

Set a data array in `var`: ie set `van.values.domainname.varname .= val`.

If Variable `domainname.varname` not present, errors if `allow_not_found==false`, or has no effect if `allow_not_found==true`.

As an optimisation where `domainname` and `varname` are known at compile time, they can be passed as Val(domainname_sym), Val(varname_sym) to maintain
type stability and avoid allocations.
"""
function set_values!(
    van::VariableAggregatorNamed, domainname::Symbol, varname::Symbol, val;
    allow_missing::Bool=false, # deprecated
    allow_not_found::Bool=allow_missing,    
)
    # check Variable exists
    if hasfield(typeof(van.values), domainname)
        if hasfield(typeof(getfield(van.values, domainname)), varname)
            getfield(getfield(van.values, domainname), varname) .= val
        else
            allow_not_found ||
                error("VariableAggregatorNamed has no Variable $domainame.$varname")
        end
    else
        allow_not_found ||
            error("VariableAggregatorNamed has no Variable $domainame.$varname")
    end

    return nothing
end

set_values!(
    van::VariableAggregatorNamed, domainname::AbstractString, varname::AbstractString, val; 
    allow_missing::Bool=false, # deprecated
    allow_not_found::Bool=allow_missing,
) = set_values!(van, Symbol(domainname), Symbol(varname), val; allow_not_found)

# type stable version
function set_values!(
    van::VariableAggregatorNamed, domainname::Val{DN}, varname::Val{VN}, val;
    allow_missing::Bool=false, # deprecated
    allow_not_found::Bool=allow_missing,
) where {DN, VN}
    # check Variable exists
    if hasfield(typeof(van.values), DN)
        if hasfield(typeof(getfield(van.values, DN)), VN)
            getfield(getfield(van.values, DN), VN) .= val
        else
            allow_not_found ||
                error("VariableAggregatorNamed has no Variable $DN.$VN")
        end
    else
        allow_not_found ||
            error("VariableAggregatorNamed has no Variable $DN.$VN")
    end

    return nothing
end