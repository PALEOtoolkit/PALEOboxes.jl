




"""
    ReactionMethod

Defines a callback function `methodfn` with Variables `vars`, to be called from the Model framework
either during setup or as part of the main loop.

# Fields
- `methodfn::Function`: callback from Model framework
- `reaction::AbstractReaction`: the Reaction that created this ReactionMethod.
- `name::String`: a descriptive name, eg generated from the name of methodfn
- `vars::Tuple{Vararg{AbstractVarList}}`: [`AbstractVarList`](@ref)s of [`VariableReaction`](@ref)s.
  Corresponding Variable accessors (views on Arrays) will be provided to the `methodfn` callback.
- `p::P`: optional context field (of arbitrary type) to store data needed by methodfn.
- `operatorID::Vector{Int64}`: operator ID (to allow operator splitting)
- `domain::Domain`: the Domain that `methodfn` operates on (ie will be called with `cellrange` for this `domain`).
- `preparefn::Function`: optional function to eg define buffers before main loop starts.

# methodfn
The `methodfn` callback is:

    methodfn(m::ReactionMethod, vardata::Tuple, cellrange::AbstractCellRange, modelctxt)

With arguments:
- `m::ReactionMethod`: context is available as `m.reaction::AbstractReaction` (the Reaction that defined the `ReactionMethod`),
  and `m.p` (an arbitrary extra context field supplied when `ReactionMethod` created).
- `vardata`: A Tuple of collections of views on Domain data arrays, defined by `vars`
- `cellrange::AbstractCellRange`: range of cells to calculate.
- `modelctxt`:
    - for a setup method, `:initial_value` or `:norm_value` defining the type of setup requested
    - for a main loop method `deltat` providing timestep information eg for rate throttling.

# preparefn
An optional `preparefn` callback can be supplied eg to allocate buffers that require knowledge of
the data types of `vardata` or to cache expensive calculations:

    preparefn(m::ReactionMethod, vardata::Tuple) -> vardata

This is called after model arrays are allocated, and prior to setup.
"""
mutable struct ReactionMethod{M, R, P, V} <: AbstractReactionMethod
    methodfn::M 
    reaction::R
    name::String
    vars::V
    p::P
    operatorID::Vector{Int64}
    domain::Domain
    preparefn::Union{Nothing, Function} # NB: not concretely typed as not performance-critical

    function ReactionMethod(
        methodfn::M,
        reaction::R,
        name,
        vars::V,
        p::P, 
        operatorID::Vector{Int64}, 
        domain; 
        preparefn = nothing
    ) where {M <: Function, R <: AbstractReaction, P, V <: Tuple{Vararg{AbstractVarList}}}        
        
        copy_vars = deepcopy(vars)

        newmethod = new{M, R, P, V}(
            methodfn, 
            reaction,
            name,
            copy_vars,
            p,
            operatorID, 
            domain, 
            preparefn
        )

        for va in copy_vars
            for v in get_variables(va)
                v.method = newmethod
                set_attribute!(v, :operatorID, newmethod.operatorID, allow_create=true)
            end
        end

        return newmethod
    end
end

"Get VariableReactions from `vars` as a Tuple"
function get_variables_tuple(method::AbstractReactionMethod)
    vars_tuple = []
    for va in method.vars
        push!(vars_tuple, get_variables(va))
    end
    return Tuple(vars_tuple)
end

"""
    get_variables(method::AbstractReactionMethod; filterfn = v -> true) -> Vector

Get VariableReactions from `method.vars` as a flat Vector, optionally restricting to those that match `filterfn`
"""
function get_variables(
    method::AbstractReactionMethod;
    filterfn = v -> true
)
    vars = VariableReaction[]
    for va in method.vars
        append!(vars, filter(filterfn, get_variables(va)))
    end
    return vars
end

"Get a single VariableReaction by localname"
function get_variable(
    method::AbstractReactionMethod, localname::AbstractString; 
    allow_not_found=false
)
    matchvars = get_variables(method; filterfn = v -> v.localname==localname)
    length(matchvars) <= 1 || error("duplicate variable localname", localname)

    !isempty(matchvars) || allow_not_found || error("method ", fullname(method), " no variable localname=", localname)

    return isempty(matchvars) ? nothing : matchvars[1]
end

fullname(method::AbstractReactionMethod) = fullname(method.reaction)*"."*method.name

is_method_setup(method::AbstractReactionMethod) = (method in method.reaction.methods_setup)
is_method_initialize(method::AbstractReactionMethod) = (method in method.reaction.methods_initialize)
is_method_do(method::AbstractReactionMethod) = (method in method.reaction.methods_do)

get_rate_stoichiometry(m::ReactionMethod) = []

###########################################
# Pretty printing
############################################

"compact form"
function Base.show(io::IO, method::ReactionMethod)
    print(
        io, 
        "ReactionMethod(fullname='", fullname(method), 
            "', methodfn=", string(method.methodfn), 
            ", domain='", method.domain.name ,
            "', operatorID=", method.operatorID,
            ")",
    )
end
