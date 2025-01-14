############################################################################################
# ParameterAggregator
###########################################################################################

"""
    ParameterAggregator(parfullnames::Vector{String}, model; eltype=Float64) -> ParameterAggregator

Represent a subset of model parameters given by `parfullnames` as a flattened Vector

`parfulnames` is a Vector of form `["domainname.reactionname.parname", ...]` defining a subset of
model parameters (NB: must be of type `ParDouble` or `ParDoubleVec` ie scalar or vector of Float64).

`norm_values` can be used to specify normalisation of the flattened parameter vector (defaults to 1.0).

The parameters can then be set from and copied to a flattened Vector using:

    copyto!(pa::ParameterAggregator, newvalues::Vector)  # set from newvalues .* norm_values
    copyto!(currentvalues::Vector, pa::ParameterAggregator) # copy to currentvalues, dividing by norm_values
    get_currentvalues(pa::ParameterAggregator) -> currentvalues::Vector

The subset of parameters are then defined by the `p` parameter Vector used by SciML solvers, and 
combined with the full set (from the yaml file) to eg solve an ODE to enable sensitivity studies.

`eltype` can be eg a Dual number to support ForwardDiff automatic differentiation for parameter Jacobians.
"""
mutable struct ParameterAggregator{T}

    # Parameter full names 
    parfullnames::Vector{String}

    # replacement Parameters to be used (order matches parfullnames)
    replacement_parameters::Vector{Union{Parameter{T, Nothing}, VecParameter{T, Nothing}}}

    # indices in flattened p Vector for each replacement parameter
    indices::Vector{UnitRange{Int64}}

    # normalization (as flattened vector) for parameters
    norm_values::Vector{Float64}

    # reactions with replacement parameter values (each entry is a Dict of :par_name => index in replacement_parameters)
    reactpars::Dict{AbstractReaction, Dict{Symbol, Int}}

    # replacement ParametersTuple (merging replacement parameters with all reaction parameters),
    # for those reactions that need parameter replacement
    reactpartuples::Dict{AbstractReaction, NamedTuple}
end

# compact form
function Base.show(io::IO, pa::ParameterAggregator)
    print(io, "ParameterAggregator(parfullnames='", pa.parfullnames,"', indices='", pa.indices,"')")
end
# multiline form
function Base.show(io::IO, ::MIME"text/plain", pa::ParameterAggregator)
    println(io, typeof(pa))
    Printf.@printf(io, "%40s%20s\n", "parfullname", "indices")
    for (pfn, i) in IteratorUtils.zipstrict(pa.parfullnames, pa.indices)
        Printf.@printf(io, "%40s%20s\n", pfn, string(i))
    end
end




function ParameterAggregator(model::AbstractModel, parfullnames::Vector{String}; eltype=Float64)

    reactpars = Dict{AbstractReaction, Dict{Symbol, Int}}()
    replacement_parameters = Vector{Union{Parameter{eltype, Nothing}, VecParameter{eltype, Nothing}}}()
    indices = UnitRange{Int}[]
    nextidx = 1

    # iterate through parfullnames and assemble lists of replacement parameters and corresponding indices in flattened vector
    for (pidx, domreactpar) in enumerate(parfullnames)
        domainname, reactionname, parname = split(domreactpar, ".")

        react = get_reaction(model, domainname, reactionname; allow_not_found=false)
        p = get_parameter(react, parname)

        if p isa Parameter{Float64, Nothing}
            replace_p = Parameter{eltype, Nothing}(
                p.name, p.description, p.units, eltype(p.v), eltype(p.default_value), eltype[], false, p.external
            )
        elseif p isa VecParameter{Float64, Nothing}
            replace_p = VecParameter{eltype, Nothing}(
                p.name, p.description, p.units, Vector{eltype}(p.v), Vector{eltype}(p.default_value), eltype[], false, p.external
            )
        else
            error("parameter $domreactpar $p is not a ParDouble or ParDoubleVec")
        end
        rparsindices = get!(reactpars, react, Dict{Symbol, Int}())
        rparsindices[Symbol(parname)] = pidx

        push!(replacement_parameters, replace_p)
        endidx = nextidx + length(replace_p) - 1
        push!(indices, nextidx:endidx)
        
        nextidx = endidx + 1

    end

    # generate new ParametersTuple for those reactions that need parameter replacement
    reactpartuples = Dict{AbstractReaction, NamedTuple}()
    for (react, rparsindices) in reactpars
        newparstuple = (haskey(rparsindices, k) ? replacement_parameters[rparsindices[k]] : v for (k, v) in pairs(react.pars))
        reactpartuples[react] = NamedTuple{keys(react.pars)}(newparstuple)
    end
  
    norm_values = ones(indices[end][end])

    return ParameterAggregator{eltype}(
        parfullnames,
        replacement_parameters,
        indices,
        norm_values,
        reactpars,
        reactpartuples,
    )
end

Base.copy(pa::ParameterAggregator{old_eltype}) where {old_eltype} = copy_new_eltype(old_eltype, pa)

function copy_new_eltype(new_eltype, pa::ParameterAggregator{old_eltype}) where {old_eltype}

    replacement_parameters = Vector{Union{Parameter{new_eltype, Nothing}, VecParameter{new_eltype, Nothing}}}()
    for p in pa.replacement_parameters
        if p isa Parameter{old_eltype, Nothing}
            replace_p = Parameter{new_eltype, Nothing}(
                p.name, p.description, p.units, new_eltype(p.v), new_eltype(p.default_value), new_eltype[], false, p.external
            )
        elseif p isa VecParameter{old_eltype, Nothing}
            replace_p = VecParameter{new_eltype, Nothing}(
                p.name, p.description, p.units, Vector{new_eltype}(p.v), Vector{new_eltype}(p.default_value), new_eltype[], false, p.external
            )
        else
            error("parameter $p is not a scalar or vector parameter with eltype $old_eltype")
        end
        push!(replacement_parameters, replace_p)
    end

    # generate new ParametersTuple for those reactions that need parameter replacement
    reactpartuples = Dict{AbstractReaction, NamedTuple}()
    for (react, rparsindices) in pa.reactpars
        newparstuple = (haskey(rparsindices, k) ? replacement_parameters[rparsindices[k]] : v for (k, v) in pairs(react.pars))
        # @Infiltrator.infiltrate
        reactpartuples[react] = NamedTuple{keys(react.pars)}(newparstuple)
    end

    pa_net = ParameterAggregator{new_eltype}(
        pa.parfullnames,
        replacement_parameters,
        pa.indices,
        pa.norm_values,
        pa.reactpars,
        reactpartuples,
    )

    return pa_net
end

# for use by solver: test whether `reaction` has modified parameters
has_modified_parameters(pa::ParameterAggregator, reaction::AbstractReaction) = haskey(pa.reactpartuples, reaction)

# for use by solver: retrieve modified parameters for `reaction`, or return `nothing` if no modified parameters
get_parameters(pa::ParameterAggregator, reaction::AbstractReaction) = get(pa.reactpartuples, reaction, nothing)

function Base.copyto!(pa::ParameterAggregator, newvalues::Vector)

    lastidx = pa.indices[end][end]
    lastidx == length(newvalues) || 
        error("ParameterAggregator length $lastidx != length(newvalues) $(length(newvalues))")

    for (p, indices) in IteratorUtils.zipstrict(pa.replacement_parameters, pa.indices)
        if p isa Parameter
            # p.v = only(view(newvalues, indices))
            setvalue!(p, only(view(newvalues, indices)) * only(view(pa.norm_values, indices)))
        elseif p isa VecParameter
            p.v .= view(newvalues, indices) .* view(pa.norm_values, indices)
        else
            error("invalid Parameter type $p")
        end
    end

    return pa
end

function Base.copyto!(currentvalues::Vector, pa::ParameterAggregator)

    lastidx = pa.indices[end][end]
    lastidx == length(currentvalues) || 
        error("ParameterAggregator length $lastidx != length(currentvalues) $(length(currentvalues))")

    for (p, indices) in IteratorUtils.zipstrict(pa.replacement_parameters, pa.indices)
        currentvalues[indices] .= p.v ./ view(pa.norm_values, indices)
    end

    return currentvalues
end

function get_currentvalues(pa::ParameterAggregator{T}) where T
    currentvalues = Vector{T}(undef, pa.indices[end][end])
    return Base.copyto!(currentvalues, pa)
end
